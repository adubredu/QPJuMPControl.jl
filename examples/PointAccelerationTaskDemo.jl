using Revise
using RigidBodyDynamics
using RigidBodyDynamics.PDControl
using RigidBodyDynamics.Graphs: target
using JuMP
using OSQP
using QPJuMPControl
using StaticArrays
using MeshCatMechanisms
using MeshCat: PointCloud, setobject!, LineSegments, LineBasicMaterial
using GeometryBasics: Point

# Load URDF
urdf = joinpath(dirname(pathof(RigidBodyDynamics)), "..", "test", "urdf", "Acrobot.urdf")
mechanism = parse_urdf(Float64, urdf)

mvis = MechanismVisualizer(mechanism, URDFVisuals(urdf))
open(mvis)

optimizer = OSQP.Optimizer

const center = Point3D(root_frame(mechanism), 0.0, 0.25, 2.2)
const radius = 0.5
const velocity = 1.0

# Create reference position, velocity, and acceleration functions
p_reference(t) = center + FreeVector3D(root_frame(mechanism), 
    radius * cos(velocity * t), 0.0, radius * sin(velocity * t))
ṗ_reference(t) = FreeVector3D(root_frame(mechanism), 
    -velocity * radius * sin(velocity * t), 0.0, velocity * radius * cos(velocity * t))
p̈_reference(t) = FreeVector3D(root_frame(mechanism),
    -velocity^2 * radius * cos(velocity * t), 0.0, -velocity^2 * radius * sin(velocity * t))


geometry = PointCloud([Point(p_reference(t).v) for t in range(0, stop=2π, length=100)])
setobject!(mvis.visualizer[:circle], LineSegments(geometry, LineBasicMaterial()));

# Create the low-level momentum-based controller
lowlevel = MomentumBasedController{4}(mechanism, optimizer)

# Add a point acceleration task to the low-level controller. The point 
# we'll control is located at the tip of the lower arm:
body = last(bodies(mechanism))
point = Point3D(default_frame(body), 0., 0, -2.05)
task = PointAccelerationTask(mechanism,
    path(mechanism, root_body(mechanism), body),
    point)
err = QPJuMPControl.task_error(task, lowlevel.qpmodel, MechanismState(mechanism), lowlevel.v̇)
@show err 

##=
# Add the task to the controller. This will create a hard constraint
# in the controller to force it to produce the desired point acceleration.
# addtask!(lowlevel, task)
# If we wanted to just add a penalty to the objective function instead,
# we could do:
    # addtask!(lowlevel, task, 1.0) 
# for some real value `cost`. 

# Also add a small regularization term to avoid unbounded joint 
# accelerations
for joint in joints(mechanism)
    regularize!(lowlevel, joint, 1e-6)
end

# Create the high-level controller. The high-level controller does the following
# at each time step:
#  1. Compute the reference position, velocity, and acceleration of the target point
#  2. Compute the desired acceleration of the point using a simple PD controller
#  3. Set the desired acceleration of the PointAccelerationTask in the low-level
#     controller
#  4. Run the low-level controller to produce the commanded torques
highlevel = let lowlevel = lowlevel, task = task, state = MechanismState(mechanism)
    function (τ, t, x)
        copyto!(state, x)
        
        # Reference position, velocity, and acceleration
        pref = p_reference(t)
        ṗref = ṗ_reference(t)
        p̈ref = p̈_reference(t)
        
        # Compute the current position and velocity of the point
        H = transform_to_root(state, target(task.path))
        T = twist_wrt_world(state, target(task.path))
        p = H * task.point
        ṗ = point_velocity(T, H * task.point)
        
        # Compute the desired acceleration using a PD law:
        p̈des = pd(PDGains(1.0, 1.0), p, pref, ṗ, ṗref) + p̈ref
        
        # Set the desired acceleration in the low-level controller
        setdesired!(task, p̈des)
        
        # Run the low-level controller to produce commanded torques  
        addtask!(lowlevel, task) 
        println(lowlevel.qpmodel)

        lowlevel(τ, t, x)
    end
end;

# Now we can run our controller from some initial state using the `simulate`
# function from RigidBodyDynamics. 
state = MechanismState(mechanism)
set_configuration!(state, [-π, 0.1])
ts, qs, vs = RigidBodyDynamics.simulate(state, 20.0, highlevel; Δt=1e-2);


# Show the location of our target point
setelement!(mvis, point)

# Animate the resulting trajectory in the visualizer
setanimation!(mvis, ts, qs)
# =#