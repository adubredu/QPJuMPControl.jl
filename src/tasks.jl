abstract type AbstractMotionTask end

struct SpatialAccelerationTask <: AbstractMotionTask 
    path::TreePath{RigidBody{Float64}, Joint{Float64}}
    jacobian::GeometricJacobian{Matrix{Float64}}
    desired::Base.RefValue{SpatialAcceleration{Float64}}

    function SpatialAccelerationTask(
            mechanism::Mechanism,
            path::TreePath{RigidBody{Float64}, Joint{Float64}};
            frame::CartesianFrame3D = default_frame(target(path)))
        nv = num_velocities(mechanism)
        bodyframe = default_frame(target(path))
        baseframe = default_frame(source(path))
        jacobian = GeometricJacobian(bodyframe, baseframe, frame, zeros(3, nv), zeros(3, nv))
        desired = Ref{SpatialAcceleration{Float64}}(zero(SpatialAcceleration{Float64}, bodyframe, baseframe, frame))
        new(path, jacobian, desired)
    end
end

dimension(task::SpatialAccelerationTask) = 6 

function setdesired!(task::SpatialAccelerationTask, desired::SpatialAcceleration)
    @framecheck task.desired[].body desired.body 
    @framecheck task.desired[].base desired.base 
    @framecheck task.desired[].frame desired.frame 
    task.desired[] = desired 
    nothing 
end

function task_error(task::SpatialAccelerationTask, qpmodel, state::MechanismState, v̇::AbstractVector)
    world_to_desired = inv(transform_to_root(state, task.desired[].frame))
    J = geometric_jacobian!(task.jacobian, state, task.path, world_to_desired)

    bias = -bias_acceleration(state, source(task.path)) + bias_acceleration(state, target(task.path))
    J̇v = transform(state, bias, task.desired[].frame)

    desired = task.desired[]

    return [
        angular(J) * v̇ + angular(J̇v) - angular(desired);
        linear(J) * v̇ + linear(J̇v) - linear(desired)
    ]
end

struct AngularAccelerationTask <: AbstractMotionTask
    path::TreePath{RigidBody{Float64}, Joint{Float64}}
    jacobian::GeometricJacobian{Matrix{Float64}}
    desired::Base.RefValue{FreeVector3D{SVector{3, Float64}}}

    function AngularAccelerationTask(
            mechanism::Mechanism,
            path::TreePath{RigidBody{Float64}, Joint{Float64}};
            frame::CartesianFrame3D = default_frame(target(path)))
        nv = num_velocities(mechanism)
        bodyframe = default_frame(target(path))
        baseframe = default_frame(source(path))
        jacobian = GeometricJacobian(bodyframe, baseframe, frame, zeros(3, nv), zeros(3, nv))
        desired = Ref(FreeVector3D(frame, 0.0,0.0,0.0))
        new(path, jacobian, desired)
    end
end

dimension(task::AngularAccelerationTask) = 3

function setdesired!(task::AngularAccelerationTask, desired::FreeVector3D)
    @framecheck task.desired[].frame desired.frame 
    task.desired[] = desired 
    nothing 
end

function task_error(task::AngularAccelerationTask, qpmodel, state::MechanismState, v̇::AbstractVector)
    world_to_desired = inv(transform_to_root(state, task.desired[].frame))
    J = geometric_jacobian!(task.jacobian, state, task.path, world_to_desired)

    bias = -bias_acceleration(state, source(task.path)) + bias_acceleration(state, target(task.path))
    J̇v = transform(state, bias, task.desired[].frame)

    desired = task.desired[].v 
    return angular(J) * v̇ + angular(J̇v) - desired 
end


struct LinearAccelerationTask <: AbstractMotionTask 
    path::TreePath{RigidBody{Float64}, Joint{Float64}}
    jacobian::GeometricJacobian{Matrix{Float64}}
    desired::Base.RefValue{FreeVector3D{SVector{3, Float64}}}

    function LinearAccelerationTask(
            mechanism::Mechanism,
            path::TreePath{RigidBody{Float64}, Joint{Float64}};
            frame::CartesianFrame3D = default_frame(target(path)))
        nv = num_velocities(mechanism)
        bodyframe = default_frame(target(path))
        baseframe = default_frame(source(path))
        jacobian = GeometricJacobian(bodyframe, baseframe, frame, zeros(3, nv), zeros(3, nv))
        desired = Ref(FreeVector3D(frame, 0.0,0.0,0.0))
        new(path, jacobian, desired)
    end
end

dimension(task::LinearAccelerationTask) = 3

function setdesired!(task::LinearAccelerationTask, desired::FreeVector3D)
    @framecheck task.desired[].frame desired.frame 
    task.desired[] = desired 
    nothing
end

function task_error(task::LinearAccelerationTask, qpmodel, state::MechanismState, v̇::AbstractVector)
    world_to_desired = inv(transform_to_root(state, task.desired[].frame))
    J = geometric_jacobian!(task.jacobian, state, task.path, world_to_desired)

    bias = -bias_acceleration(state, source(task.path)) + bias_acceleration(state, target(task.path))
    J̇v = transform(state, bias, task.desired[].frame)

    desired = task.desired[].v 
    return linear(J) * v̇ + linear(J̇v) - desired 
end


struct PointAccelerationTask <: AbstractMotionTask
    path::TreePath{RigidBody{Float64}, Joint{Float64}}
    jacobian::PointJacobian{Matrix{Float64}}
    point::Point3D{SVector{3, Float64}}
    desired::Base.RefValue{FreeVector3D{SVector{3, Float64}}}

    function PointAccelerationTask(
            mechanism::Mechanism,
            path::TreePath{RigidBody{Float64}, Joint{Float64}},
            point::Point3D)
        nv = num_velocities(mechanism)
        bodyframe = default_frame(target(path))
        baseframe = default_frame(source(path))
        @framecheck point.frame bodyframe 
        jacobian = PointJacobian(baseframe, zeros(3, nv))
        desired = Ref(FreeVector3D(baseframe, 0.0, 0.0, 0.0))
        new(path, jacobian, point, desired)
    end
end

dimension(task::PointAccelerationTask) = 3

function setdesired!(task::PointAccelerationTask, desired::FreeVector3D)
    @framecheck task.desired[].frame desired.frame 
    task.desired[] = desired 
    nothing
end

function task_error(task::PointAccelerationTask, qpmodel, state::MechanismState, v̇::AbstractVector)
    frame = task.desired[].frame 
    point_in_task_frame = transform(state, task.point, frame)

    J = point_jacobian!(task.jacobian, state, task.path, point_in_task_frame)

    bias = -bias_acceleration(state, source(task.path)) + bias_acceleration(state, target(task.path))
    J̇v = transform(state, bias, frame)

    desired =  task.desired[].v 

    T = transform(state, relative_twist(state, target(task.path), source(task.path)), frame)
    @framecheck T.frame frame 
    ω = angular(T)
    ṗ = point_velocity(T, point_in_task_frame)
    return ω × ṗ.v + J.J * v̇ + (angular(J̇v) × point_in_task_frame.v + linear(J̇v)) - desired 
end


struct JointAccelerationTask{JT<:JointType{Float64}} <: AbstractMotionTask 
    joint::Joint{Float64, JT}
    desired::Vector{Float64}

    function JointAccelerationTask(joint::Joint{Float64, JT}) where {JT<:JointType{Float64}}
        new{JT}(joint, zeros(num_velocities(joint)))
    end
end

dimension(task::JointAccelerationTask) = length(task.desired)
setdesired!(task::JointAccelerationTask, desired) = set_velocity!(task.desired, task.joint, desired)

function task_error(task::JointAccelerationTask, qpmodel, state::MechanismState, v̇::AbstractVector)
    desired = task.desired 
    # @show desired
    v̇joint = v̇[velocity_range(state, task.joint)]
    return v̇joint - desired 
end


struct MomentumRateTask <: AbstractMotionTask
    momentum_matrix::MomentumMatrix{Matrix{Float64}}
    desired::Base.RefValue{Wrench{Float64}}

    function MomentumRateTask(mechanism::Mechanism, centroidalframe::CartesianFrame3D)
        nv = num_velocities(mechanism)
        momentum_matrix = MomentumMatrix(centroidalframe, zeros(3, nv), zeros(3, nv))
        desired = Ref(zero(Wrench{Float64}, centroidalframe))
        new(momentum_matrix, desired)
    end
end

function momentum_rate_task_params(task, qpmodel, state, v̇)
    centroidalframe = task.momentum_matrix.frame 
    com = center_of_mass(state)
    centroidal_to_world = Transform3D(centroidalframe, com.frame, com.v)
    world_to_centroidal = inv(centroidal_to_world)
    A = momentum_matrix!(task.momentum_matrix, state, world_to_centroidal)
    
    Ȧv = transform(momentum_rate_bias(state), world_to_centroidal)

    return A, Ȧv 
end

dimension(task::MomentumRateTask) = 6

function setdesired!(task::MomentumRateTask, desired::Wrench)
    @framecheck task.momentum_matrix.frame desired.frame 
    task.desired[] = desired 
end

function task_error(task::MomentumRateTask, qpmodel, state::MechanismState, v̇::AbstractVector)
    A, Ȧv = momentum_rate_task_params(task, qpmodel, state, v̇)
    desired = task.desired[]
    return [angular(A) * v̇ + angular(Ȧv) - angular(desired);
            linear(A) * v̇ + linear(Ȧv) - linear(desired)]
end


struct LinearMomentumRateTask <: AbstractMotionTask
    momentum_matrix::MomentumMatrix{Matrix{Float64}}
    desired::Base.RefValue{FreeVector3D{SVector{3, Float64}}}

    function LinearMomentumRateTask(mechanism::Mechanism, centroidalframe::CartesianFrame3D = CartesianFrame3D())
        nv = num_velocities(mechanism)
        momentum_matrix = MomentumMatrix(centroidalframe, zeros(3,nv), zeros(3, nv))
        desired = Ref(FreeVector3D(centroidalframe, 0.0, 0.0, 0.0))
        new(momentum_matrix, desired)
    end
end

dimension(task::LinearMomentumRateTask) = 3 

function setdesired!(task::LinearMomentumRateTask, desired::FreeVector3D)
    @framecheck task.momentum_matrix.frame desired.frame 
    task.desired[] = desired 
end

function task_error(task::LinearMomentumRateTask, qpmodel, state::MechanismState, v̇::AbstractVector)
    A, Ȧ = momentum_rate_task_params(task, qpmodel, state, v̇)
    desired = task.desired[]
    return linear(A) * v̇ + linear(Ȧv) - desired
end
