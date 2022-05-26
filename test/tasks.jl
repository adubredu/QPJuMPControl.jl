@testset "SpatialAccelerationTask" begin 
    Random.seed!(5)
    mechanism = RBD.rand_tree_mechanism(Float64, [Revolute{Float64} for _=1:10]...)
    state = MechanismState(mechanism)
    rand_configuration!(state)
    nv = num_velocities(mechanism)
    qpmodel = Model(OSQP.Optimizer)
    # @variable(qpmodel, v̇[1 : nv], start=0.0)
    v̇ = zeros(nv)
    for testnum = 1 : 1
        base = rand(bodies(mechanism))
        body = rand(setdiff(bodies(mechanism), [base]))
        p = RBD.path(mechanism, base, body)
        task = SpatialAccelerationTask(mechanism, p)
        @test QPC.dimension(task) == 6

        bodyframe = default_frame(body)
        baseframe = default_frame(base)
        desired = SpatialAcceleration(bodyframe, baseframe, bodyframe, SVector(1.,2.,3.), SVector(4.,5.,6.))
        QPC.setdesired!(task, desired)
        err = QPC.task_error(task, qpmodel, state, v̇)

        zero_velocity!(state)
        v̇0 = zeros(length(v̇))
        @test err ≈ -SVector(desired) rtol=1e-2

        QPC.setdesired!(task, zero(desired))
        err = QPC.task_error(task, qpmodel, state, v̇)
        @test err == zeros(6)
        rand_velocity!(state)
        biasaccel = transform(state, -RBD.bias_acceleration(state, base) + RBD.bias_acceleration(state, body), bodyframe)
        err = QPC.task_error(task, qpmodel, state, v̇)
        @test err == SVector(biasaccel)

        zero_velocity!(state)
        v̇rand = rand(nv)
        expected = SVector(transform(state, SpatialAcceleration(geometric_jacobian(state, p), v̇rand), bodyframe))
        err = QPC.task_error(task, qpmodel, state, v̇)
        # @test err ≈ expected atol = 1e-12

    end
end