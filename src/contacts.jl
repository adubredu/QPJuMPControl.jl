"""
z_up_transform(origin::Point3D, zaxis::FreeVector3D, from::CartesianFrame3D)
Return `Transform3D` from `from` to the frame in which both `origin` and `zaxis`
is expressed, such that z-axis of `from` is `zaxis`, and origin of `from` is
`origin`.
"""
function z_up_transform(origin::Point3D, zaxis::FreeVector3D,       from::CartesianFrame3D)
    @framecheck origin.frame zaxis.frame
    to = origin.frame
    rotation = Rotations.rotation_between(SVector(0., 0., 1.), zaxis.v)
    translation = origin.v
    Transform3D(from, to, rotation, translation)
end

function forcebasis(μ::Float64, num_basis_vectors::Val{N}) where N
    Δθ = 2 * π / N
    basisvectors = ntuple(num_basis_vectors) do i
        θ = (i - 1) * Δθ
        normalize(SVector(μ * cos(θ), μ * sin(θ), 1.0))
    end
    hcat(basisvectors...)
end


struct ContactPoint{N}
    normal_aligned_frame::CartesianFrame3D
    ρ::Vector{VariableRef}
    force_local::FreeVector3D{Vector{VariableRef}}
    wrench_world::Wrench{VariableRef}
    position::Any 
    normal::Any 
    μ::Any 
    weight::Ref{Float64}
    maxnormalforce::Ref{Float64} 

    function ContactPoint{N}(position::Point3D, normal::FreeVector3D, μ::Float64, state::MechanismState, model::JuMP.Model) where N 
        
        #frames 
        normal_aligned_frame = CartesianFrame3D()
        worldframe = root_frame(state.mechanism)

        # variables 
        ρ = @variable(model, [1:N])
        force_local_vars = @variable(model, [1:3])
        force_local = FreeVector3D(normal_aligned_frame, force_local_vars) 
        wrench_angular = @variable(model, [1:3])
        wrench_linear = @variable(model, [1:3])
        wrench_world = Wrench(worldframe, wrench_angular, wrench_linear)

        ret = new{N}(normal_aligned_frame, ρ, force_local, wrench_world, position, normal, μ, Ref(0.0), Ref(0.0))

        basis = forcebasis(μ, Val(N))
        maxnormalforce = ret.maxnormalforce[]
        maxρ = (maxnormalforce / (N * sqrt(μ*μ + 1))) * ones(N)
        state_param = state 
        toroot = transform_to_root(state_param, position.frame) * z_up_transform(position, normal, normal_aligned_frame)
        hat = RBD.Spatial.hat 

        @constraint(model, force_local.v .== basis * ρ)
        @constraint(model, ρ .>= zeros(N))
        @constraint(model, ρ .<= maxρ)
        @constraint(model, linear(wrench_world) .== rotation(toroot)*force_local.v)
        @constraint(model, angular(wrench_world) .== hat(translation(toroot)) * linear(wrench_world))

        return ret 
    end
end 

disable!(point::ContactPoint) = point.maxnormalforce[] = 0 
isenabled(point::ContactPoint) = point.maxnormalforce[] > 0

function objectiveterm(point::ContactPoint, model::JuMP.Model)
    weight = point.weight[]
    f = point.force_local
    return weight * (f ⋅ f)
end

function JuMP.value(m::JuMP.Model, wrench::Wrench{VariableRef})
    τ = angular(wrench)
    f = linear(wrench)
    @inbounds τval = SVector(value.(τ)...)
    @inbounds fval = SVector(value.(f)...)
    return Wrench(wrench.frame, τval, fval)
end

    