struct Interpolated{X, Y, I}
    x0::X
    xf::X
    y0::Y
    yf::Y
    interpolator::I
    clamp::Bool

    function Interpolated(x0::X, xf::X, y0::Y, yf::Y, f=identity; min_num_derivs=Val(0), clamp::Bool=true) where {X, Y}
        interpolator = make_interpolator(f, min_num_derivs)
        new{X, Y, typeof(interpolator)}(x0, xf, y0, yf, interpolator, clamp)
    end
end

function Interpolated(x0, xf, y0, yf, f=identity; min_num_derivs=Val(0))
    Interpolated(promote(x0, xf)..., promote(y0, yf)..., f, min_num_derivs=min_num_derivs)
end

# Call overload with derivatives
function (trajectory::Interpolated)(x, ::Val{num_derivs}) where num_derivs
    x0, xf = trajectory.x0, trajectory.xf
    if !trajectory.clamp && (x < x0 || x > xf)
        throw_trajectory_domain_error(x, x0, xf)
    end

    # 1. Compute θ(x) = clamp((x - x0) / (xf - x0), 0, 1) and dθ/dx.
    Δx = xf - x0
    θ = (x - x0) / Δx
    if θ <= zero(θ)
        θ = zero(θ)
        dθdx = zero(θ)
    elseif θ >= one(θ)
        θ = one(θ)
        dθdx = zero(θ)
    else
        dθdx = one(x - x0) / Δx
    end

    # 2. Warp θ using interpolator to obtain α(θ(x)) and its partial derivatives w.r.t θ.
    α, α_θ_derivs = let α_vals = trajectory.interpolator(θ, Val(num_derivs))
        α_vals[1], Base.tail(α_vals)
    end

    # 3. Compute derivatives of α w.r.t x using higher-order chain rule,
    # where higher order derivatives of θ w.r.t. x are zero.
    α_x_derivs = let α_θ_derivs=α_θ_derivs, dθdx=dθdx
        ntuple(i -> α_θ_derivs[i] * dθdx^i, Val(num_derivs))
    end

    # 3. Interpolate between y0 and yf using α.
    y0, yf = trajectory.y0, trajectory.yf
    y, dydα = interpolate(y0, yf, α)

    # 4. Compute derivatives of y w.r.t. x using higher-order chain rule,
    # where higher order derivatives of y w.r.t. α are zero.
    y_x_derivs = ntuple(i -> dydα * α_x_derivs[i], Val(num_derivs))

    y, y_x_derivs...
end

# Single-argument call overload just returns the value:
(trajectory::Interpolated)(x) = trajectory(x, Val(0))[1]

# Interpolation
function interpolate(y0, yf, α)
    Δy = yf - y0
    y = y0 + α * Δy
    dydα = Δy
    y, dydα
end

function interpolate(y0::Rotation{3}, yf::Rotation{3}, α)
    interpolate(QuatRotation(y0), QuatRotation(yf), α)
end

function interpolate(y0::QuatRotation, yf::QuatRotation, α)
    Δy = AngleAxis(y0 \ yf)
    angle = rotation_angle(Δy)
    axis = rotation_axis(Δy)
    y = y0 * QuatRotation(AngleAxis(α * angle, axis...))
    dydα = axis * angle # Lie derivative
    y, dydα
end

# `make_interpolator` overloads
# function returned by `make_interpolator` has signature (θ, n) -> (f(θ), (f′(θ), ...,f^{(n)}(θ)))
function make_interpolator(f, ::Val)
    function (θ, ::Val{num_derivs}) where {num_derivs}
        # Default to calling the `nthderiv` function
        f(θ), ntuple(n -> nthderiv(f, θ, n), Val(num_derivs))...
    end
end

function make_interpolator(p::Polynomial, ::Val{min_num_derivs}) where min_num_derivs
    PolyDerivEvaluator(p, Val(min_num_derivs))
end

# `nthderiv` overloads
nthderiv(::typeof(identity), θ, n::Int) = n === 1 ? one(θ) : zero(θ)