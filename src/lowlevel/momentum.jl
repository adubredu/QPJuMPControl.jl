mutable struct MomentumBasedController{N, S<:MechanismState}
    state::S 
    result::DynamicsResult{Float64, Float64}
    floatingjoint::Union{Nothing, Joint{Float64}}
    floating_joint_velocity_range::UnitRange{Int}
    centroidalframe::CartesianFrame3D
    momentum_matix::MomentumMatrix{Matrix{Float64}}
    contacts::Dict{RigidBody{Float64}, Vector{ContactPoint{N}}}
    contactwrenches::Dict{BodyID, Wrench{Float64}}
    qpmodel::JuMP.Model
    v̇::Vector{VariableRef}
    objective::QuadExpr
    initialized::Bool 

    function MomentumBasedController{N}(
        mechanism::Mechanism{Float64}, optimizer; floatingjoint=nothing)   where {N}
        state = MechanismState(mechanism)
        result = DynamicsResult(mechanism)
        floating_joint_velocity_range = floatingjoint === nothing ? (1 : 0) : velocity_range(state, floatingjoint)
        worldframe = root_frame(mechanism)
        centroidalframe = CartesianFrame3D("centroidal")
        nv = num_velocities(state)
        momentum_matix = MomentumMatrix(worldframe, zeros(3, nv), zeros(3, nv)) 
        contacts = Dict{RigidBody{Float64}, Vector{ContactPoint{N}}}()
        contactwrenches = Dict{BodyID, Wrench{Float64}}()
        qpmodel = JuMP.Model(optimizer)
        # set_silent(qpmodel)
        @variable(qpmodel, v̇[1:nv])
        objective = @expression(qpmodel, 0 * v̇[1]^2) #add_to_expression!(ey, 1.0, y[2], y[2])

        new{N, typeof(state)}(state, result, floatingjoint, floating_joint_velocity_range, centroidalframe, momentum_matix, contacts, contactwrenches, qpmodel, v̇, objective, false)
    end
end



Base.show(io::IO, controller::MomentumBasedController{N, S}) where {N, S} = print(io, "MomentumBasedController{$N, $S}(…)")

centroidal_frame(controller::MomentumBasedController) = controller.centroidalframe


function (controller::MomentumBasedController)(τ::AbstractVector, t::Number, x::Union{<:Vector, <:MechanismState})
    if !controller.initialized 
        initialize!(controller)
        controller.initialized = true
    end

    qpmodel = controller.qpmodel
    state = controller.state 
    objective = controller.objective
    result = controller.result 
    contacts = controller.contacts 
    contactwrenches = controller.contactwrenches
    worldframe = root_frame(state.mechanism)

    copyto!(state, x)
    @objective(qpmodel, Min, objective)
    optimize!(qpmodel)
    
    vals = value.(controller.v̇)
    @inbounds for i in eachindex(vals)
        result.v̇[i] = vals[i]
        # @show vals[i]
    end
    empty!(contactwrenches)
    for body in keys(contacts)
        contactwrench = RBD.zero(Wrench{Float64}, worldframe)
        for data in contacts[body]
            # wl = value.(wrench_linear)
            # wa = value.(wrench_angular)
            wl = value.(linear(data.wrench_world))
            wa = value.(angular(data.wrench_world))
            # @show wa, wl
            wrench_world = Wrench(worldframe, wa, wl)
            contactwrench += wrench_world 
        end
        # @show typeof(contactwrench)
        contactwrenches[BodyID(body)] = contactwrench 
    end

    inverse_dynamics!(τ, result.jointwrenches, result.accelerations, state, result.v̇, contactwrenches)

    zero_floating_joint_torques!(τ, state, controller.floating_joint_velocity_range)

    return τ
end

function zero_floating_joint_torques!(τ::AbstractVector, state::MechanismState, velocity_range::UnitRange)
    for i in velocity_range
        τ[i] = 0
    end
end

function addtask!(controller::MomentumBasedController, task::AbstractMotionTask)
    model = controller.qpmodel 
    err = task_error(task, model, controller.state, controller.v̇)
    taskdim = dimension(task)
    @constraint(model, err.==zeros(taskdim))
    nothing 
end

function addtask!(controller::MomentumBasedController, task::AbstractMotionTask, weight::Number)
    e = add_task_error_slack_variables!(controller, task)
    add_to_expression!(controller.objective, weight, e ⋅ e)
    return e 
end

function addtask!(controller::MomentumBasedController, task::AbstractMotionTask, weight::AbstractMatrix)
    e = add_task_error_slack_variables!(controller, task)
    add_to_expression!(controller.objective, transpose(e), weight, e)
    return e 
end

function add_task_error_slack_variables!(controller::MomentumBasedController, task::AbstractMotionTask)
    model = controller.qpmodel 
    state = controller.state 
    v̇ = controller.v̇
    e = @variable(model, [1:dimension(task)])
    @constraint(model, e .== task_error(task, model, state, v̇))
    return e 
end

function regularize!(controller::MomentumBasedController, joint::Joint, weight)
    v̇joint = controller.v̇[velocity_range(controller.state, joint)]
    add_to_expression!(controller.objective, weight, v̇joint ⋅ v̇joint)
end

function addcontact!(controller::MomentumBasedController{N}, body::RigidBody{Float64}, point::ContactPoint{N}) where N 
    push!(get!(Vector{ContactPoint{N}}, controller.contacts, body), point)
    objterm = objectiveterm(point, controller.qpmodel)
    add_to_expression!(controller.objective, objterm)
    return point
end

function addcontact!(controller::MomentumBasedController{N}, body::RigidBody{Float64}, position::Point3D, normal::FreeVector3D,  μ::Float64) where N
    addcontact!(controller, body, ContactPoint{N}(position, normal, μ, controller.state, controller.qpmodel))
end

@noinline function initialize!(controller::MomentumBasedController)
setobjective!(controller)
    mechanism = controller.state.mechanism
    if controller.floatingjoint !== nothing
        add_wrench_balance_constraint!(controller, controller.floatingjoint)
    end
end

function setobjective!(controller::MomentumBasedController)
    @objective(controller.qpmodel, Min, controller.objective)
end 

function add_wrench_balance_constraint!(controller::MomentumBasedController{N}, joint::Joint) where N 
    qpmodel = controller.qpmodel
    state = controller.state
    mechanism = state.mechanism
    fg = mass(mechanism) * mechanism.gravitational_acceleration
    v̇ = controller.v̇
    A = momentum_matrix!(controller.momentum_matix, state)
    Ȧv = momentum_rate_bias(state)
    Wg = Wrench(center_of_mass(state) × fg, fg)
    floatingbody = successor(joint, mechanism)
    qjoint = configuration(state, joint)
    H = transform_to_root(state, floatingbody)
    S = transform(motion_subspace(joint, qjoint), H)

    torque = angular(Wg)
    force = linear(Wg)
    for contactvec in values(controller.contacts)
        for contact in contactvec 
            torque = torque + angular(contact.wrench_world)
            forcce = force + linear(contact.wrench_world)
        end 
    end
    nv = num_velocities(joint)
    @constraint(qpmodel, transpose(angular(S)) * (angular(A) * v̇ + angular(Ȧv) - torque) + transpose(linear(S)) * (linear(A) * v̇ + linear(Ȧv) - force) .== zeros(nv))

    nothing 
end