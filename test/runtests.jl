using QPJuMPControl
using Test
using RigidBodyDynamics
using JuMP 
using OSQP 
using Random 
using LinearAlgebra
using StaticArrays
using ValkyrieRobot
using ValkyrieRobot.BipedControlUtil
using Rotations 
using RigidBodyDynamics.PDControl
using RigidBodyDynamics.Contact
using StaticUnivariatePolynomials

const QPC = QPJuMPControl 
const RBD = RigidBodyDynamics

# include("tasks.jl")
include("trajectories.jl")
# include("controller.jl")

