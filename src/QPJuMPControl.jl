module QPJuMPControl

using LinearAlgebra
using RigidBodyDynamics
using RigidBodyDynamics.Graphs
using RigidBodyDynamics.Contact 
using RigidBodyDynamics.PDControl
using StaticArrays
using JuMP
using OSQP 
using Rotations 
using StaticUnivariatePolynomials

const RBD = RigidBodyDynamics

include("contacts.jl")
include("tasks.jl")
include("exceptions.jl")
include(joinpath("lowlevel", "momentum.jl"))
include(joinpath("lowlevel", "se3pdcontroller.jl"))
include(joinpath("highlevel", "standing.jl"))
include(joinpath("trajectories", "trajectories.jl")) 

# Task-related
export 
    AbstractMotionTask,
    SpatialAccelerationTask,
    AngularAccelerationTask,
    LinearAccelerationTask,
    LinearMomentumRateTask,
    JointAccelerationTask,
    MomentumRateTask,
    PointAccelerationTask, 
    setdesired!

# Low-level 
export 
    MomentumBasedController,
    SE3PDController,
    addtask!,
    addcontact!,
    disable!,
    regularize!,
    centroidal_frame 

# High-level 
export 
    StandingController

export
    SE3Trajectory

end
