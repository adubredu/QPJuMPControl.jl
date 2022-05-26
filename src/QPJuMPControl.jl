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

const RBD = RigidBodyDynamics

include("tasks.jl")

export 
    SpatialAccelerationTask

end
