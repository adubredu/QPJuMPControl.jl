using QPJuMPControl
using Test
using RigidBodyDynamics
using JuMP 
using OSQP 
using Random 
using StaticArrays

const QPC = QPJuMPControl 
const RBD = RigidBodyDynamics

include("tasks.jl")

