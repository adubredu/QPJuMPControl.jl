@testset "SE3PDController" begin
    val = Valkyrie()
    mechanism = val.mechanism
    state = MechanismState(mechanism)
    rng = MersenneTwister(1)
    base = rand(rng, bodies(mechanism))
    body = rand(rng, bodies(mechanism))
    angular = let angle = π / 2, axis = SVector(1.0, 0.0, 0.0), y0 = one(QuatRotation), yf = QuatRotation(AngleAxis(angle, axis...))
        QPC.Trajectories.Interpolated(0.0, 1.0, y0, yf)
    end
    linear = QPC.Trajectories.Interpolated(0.0, 1.0, SVector(0.0, 1.0, 2.0), SVector(2.0, 3.0, 4.0))
    frame = default_frame(body)
    gains = SE3PDGains(FramePDGains(frame, PDGains(100.0, 20.0)), FramePDGains(frame, PDGains(1000.0, 200.0)))
    trajectory = QPC.Trajectories.SE3Trajectory(default_frame(body), default_frame(base), angular, linear)
    weight = Diagonal(vcat(zeros(3), fill(10.0, 3)))
    controller = SE3PDController(BodyID(base), BodyID(body), trajectory, weight, gains)
    controller(0.5, state) 
end

@testset "parameterized contacts" begin
    # Construct a mechanism consisting of a single body which can
    # rotate about its origin
    world = RigidBody{Float64}("world")
    mechanism = Mechanism(world, gravity=SVector(0, 0, -9.81))
    frame = CartesianFrame3D("body")
    inertia = SpatialInertia(frame, SDiagonal(1., 1, 1), SVector(0., 0, 0), 10.0)
    body = RigidBody(inertia)
    joint = Joint("rx", Revolute(SVector(1., 0, 0)))
    attach!(mechanism, world, body, joint)

    contactmodel = SoftContactModel(hunt_crossley_hertz(k = 500e3), ViscoelasticCoulombModel(0.8, 20e3, 100.))
    add_contact_point!(body, Contact.ContactPoint(Point3D(default_frame(body), 0., 0, 0), contactmodel))

    # Make the contact normal a parameter, so that it updates its representation
    # in body frame to match a fixed orientation in world frame
    state = MechanismState(mechanism)
    model = Model(OSQP.Optimizer)
    set_silent(model)
    position = Point3D(default_frame(body), 0., 0, 0)
    μ = 1.0 
    normal = transform(state, FreeVector3D(root_frame(state.mechanism), 0., 0, 1), default_frame(body))
    controller_contact = QPJuMPControl.ContactPoint{4}(position, normal, μ, state, model)
    controller_contact.maxnormalforce[] = 1e3
    @test QPJuMPControl.isenabled(controller_contact)

    # Constrain the contact force to be non-zero for testing
    @constraint(model, controller_contact.force_local.v .== [0.0, 0.0, 1.0])

    for θ in range(-π, stop=π, length=2)
        set_configuration!(state, [θ])
        optimize!(model)
        # Sanity check our constraint
        @test value.(controller_contact.force_local.v) ≈ [0.0, 0.0, 1.0] atol=1e-3

        # No matter how we rotate the robot, the contact-aligned frame will still
        # be aligned to the contact normal, which is fixed in world frame. So the
        # linear component of the contact wrench in world frame will always be
        # along [0, 0, 1].
        @test value.(linear(controller_contact.wrench_world)) ≈ [0.0, 0.0, 1.0] atol=1e-3
    end
end