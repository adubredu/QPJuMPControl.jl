using QPJuMPControl.Trajectories
using Test
using LinearAlgebra
using Random
using Rotations
using RigidBodyDynamics
using StaticArrays
using StaticUnivariatePolynomials
using ForwardDiff

using QPJuMPControl.Trajectories: subfunctions, breaks
using StaticUnivariatePolynomials: exponential_integral

const SUP = StaticUnivariatePolynomials

@testset "fit_cubic" begin
    rng = MersenneTwister(15)
    for i = 1 : 10
        x0 = rand(rng)
        xf = rand(rng)
        y0 = rand(rng)
        yd0 = rand(rng)
        yf = rand(rng)
        ydf = rand(rng)

        p = fit_cubic(x0=x0, xf=xf, y0=y0, yd0=yd0, yf=yf, ydf=ydf)
        pd = SUP.derivative(p)
        pdd = SUP.derivative(pd)

        atol = 1e-6
        @test p(x0) ≈ y0 atol=atol
        @test pd(x0) ≈ yd0 atol=atol
        @test p(xf) ≈ yf atol=atol
        @test pd(xf) ≈ ydf atol=atol

        @test length(p.coeffs) == 4
    end
end

@testset "fit_quintic" begin
    rng = MersenneTwister(15)
    for i = 1 : 10
        x0 = rand(rng)
        xf = rand(rng)
        y0 = rand(rng)
        yd0 = rand(rng)
        ydd0 = rand(rng)
        yf = rand(rng)
        ydf = rand(rng)
        yddf = rand(rng)

        p = fit_quintic(x0=x0, xf=xf, y0=y0, yd0=yd0, ydd0=ydd0, yf=yf, ydf=ydf, yddf=yddf)
        pd = SUP.derivative(p)
        pdd = SUP.derivative(pd)

        atol = 1e-6
        @test p(x0) ≈ y0 atol=atol
        @test pd(x0) ≈ yd0 atol=atol
        @test pdd(x0) ≈ ydd0 atol=atol
        @test p(xf) ≈ yf atol=atol
        @test pd(xf) ≈ ydf atol=atol
        @test pdd(xf) ≈ yddf atol=atol

        @test length(p.coeffs) == 6
    end
end

@testset "Interpolated, identity interpolator function" begin
    let y0 = 1, yf = 2
        traj = Interpolated(0.0, 1.0, y0, yf)

        y = traj(0.5)
        @test y === 1.5

        y, yd, ydd = traj(0.5, Val(2))
        @test y === 1.5
        @test yd === 1.0
        @test ydd === 0.0
    end

    let f = CartesianFrame3D(), y0 = Point3D(f, 1, 2, 3), yf = Point3D(f, 0, 2, 4)
        traj = Interpolated(0.0, 1.0, y0, yf)
        y, yd, ydd = traj(0.5, Val(2))
        @test y ≈ Point3D(f, 0.5, 2.0, 3.5) atol=1e-15
        @test yd ≈ FreeVector3D(f, -1.0, 0.0, 1.0) atol=1e-15
        @test ydd ≈ FreeVector3D(f, 0.0, 0.0, 0.0) atol=1e-15
    end

    let angle = π / 2, axis = SVector(1.0, 0.0, 0.0), y0 = one(QuatRotation), yf = QuatRotation(AngleAxis(angle, axis...))
        traj = Interpolated(0.0, 1.0, y0, yf)

        y = traj(0.0)
        @test y ≈ y0 atol=1e-15

        y = traj(1.0)
        @test y ≈ yf atol=1e-15

        y, yd, ydd = traj(0.5, Val(2))
        @test y ≈ QuatRotation(AngleAxis(0.5 * angle, axis...)) atol=1e-15
        @test normalize(yd) ≈ axis atol=1e-15
        @test norm(yd) ≈ abs(angle) atol=1e-15
        @test ydd ≈ SVector(0.0, 0.0, 0.0) atol=1e-15
    end
end

@testset "Interpolated, polynomial interpolator function" begin
    interpolator = fit_quintic(x0=0.0, xf=1.0, y0=0.0, yd0=0.0, ydd0=0.0, yf=1.0, ydf=0.0, yddf=0.0)
    x0 = -1.0
    xf = 2.0
    y0 = 2.0
    yf = 3.0
    traj = Interpolated(x0, xf, y0, yf, interpolator, min_num_derivs=Val(2))

    y, yd, ydd = traj(x0, Val(2))
    @test y ≈ y0 atol=1e-10
    @test yd ≈ 0.0 atol=1e-10
    @test ydd ≈ 0.0 atol=1e-10

    y, yd, ydd = traj(xf, Val(2))
    @test y ≈ yf atol=1e-10
    @test yd ≈ 0.0 atol=1e-10
    @test ydd ≈ 0.0 atol=1e-10

    for x in range(x0; stop=xf, length=10)
        _, yd, ydd = traj(x, Val(2))
        yd_forward = ForwardDiff.derivative(traj, x)
        ydd_forward = ForwardDiff.derivative(x -> ForwardDiff.derivative(traj, x), x)
        @test yd ≈ yd_forward atol=1e-4
        @test ydd ≈ ydd_forward atol=1e-4
    end
end

function test_piecewise(traj::Piecewise)
    n = length(subfunctions(traj))
    for (i, t) in enumerate(breaks(traj))
        subtraj = subfunctions(traj)[min(i, n)]
        @test traj(t) == subtraj(t - breaks(traj)[min(i, n)])
        if i < n
            tnext = breaks(traj)[i + 1]
            tmid = (t + tnext) / 2
            @test traj(tmid) == subtraj(tmid - t)
            @test traj(tmid, Val(2)) == subtraj(tmid - t, Val(2))
        end
    end
    t0 = first(breaks(traj))
    tf = last(breaks(traj))
    if traj.clamp
        @test traj(t0 - 1) == traj(t0)
        @test traj(tf + 1) == traj(tf)
    else
        @test_throws DomainError traj(t0 - 1)
        @test_throws DomainError traj(tf + 1)
    end
end

@testset "Piecewise constant" begin
    n = 5
    subtrajectories = [Constant(i) for i = 1 : n]
    breaks = [i^2 - 1 for i = 1 : n + 1]
    test_piecewise(Piecewise(subtrajectories, breaks; clamp=true))
    test_piecewise(Piecewise(subtrajectories, breaks; clamp=false))
end

@testset "Piecewise interpolated" begin
    rng = MersenneTwister(1)
    n = 5
    breaks = [i^2 - 1 for i = 1 : n + 1]
    subtrajectories = map(diff(breaks)) do Δt
        Interpolated(0.0, Float64(Δt), rand(rng), rand(rng))
    end
    test_piecewise(Piecewise(subtrajectories, breaks; clamp=true))
    test_piecewise(Piecewise(subtrajectories, breaks; clamp=false))
end

@testset "PointTrajectory, FreeVectorTrajectory" begin
    rng = MersenneTwister(2)
    interpolator = fit_cubic(x0=0.0, xf=1.0,
        y0=rand(rng), yd0=rand(rng), yf=rand(rng), ydf=rand(rng))
    x0 = -1.0
    xf = 2.0
    y0 = rand(rng, SVector{3})
    yf = rand(rng, SVector{3})
    interpolated = Interpolated(x0, xf, y0, yf, interpolator, min_num_derivs=Val(2))
    frame = CartesianFrame3D()
    t = 1.0
    let
        traj = PointTrajectory(frame, interpolated)
        @test traj(t) === Point3D(frame, interpolated(t))
        point, vel, accel = traj(t, Val(2))
        @test vel === FreeVector3D(frame, interpolated(t, Val(2))[2])
        @test accel === FreeVector3D(frame, interpolated(t, Val(2))[3])
        allocs = let traj=traj, t=t
            @allocated traj(t, Val(2))
        end
        @test allocs == 0
    end
    let
        traj = FreeVectorTrajectory(frame, interpolated)
        @test traj(t) === FreeVector3D(frame, interpolated(t))
        x, xd, xdd = traj(t, Val(2))
        @test xd === FreeVector3D(frame, interpolated(t, Val(2))[2])
        @test xdd === FreeVector3D(frame, interpolated(t, Val(2))[3])
        allocs = let traj=traj, t=t
            @allocated traj(t, Val(2))
        end
        @test allocs == 0
    end
end

@testset "SE(3)" begin
    angular = let angle = π / 2, axis = SVector(1.0, 0.0, 0.0), y0 = one(QuatRotation), yf = QuatRotation(AngleAxis(angle, axis...))
        Interpolated(0.0, 1.0, y0, yf)
    end
    linear = Interpolated(0.0, 1.0, SVector(0.0, 1.0, 2.0), SVector(2.0, 3.0, 4.0))
    body = CartesianFrame3D()
    base = CartesianFrame3D()
    frame = body
    traj = SE3Trajectory(body, base, angular, linear)
    H, T, Td = traj(0.5, Val(2))
    @test H.from == frame
    @test H.to == base
    @test T.body == body
    @test T.base == base
    @test T.frame == frame
    @test Td.body == body
    @test Td.base == base
    @test Td.frame == frame
end
 