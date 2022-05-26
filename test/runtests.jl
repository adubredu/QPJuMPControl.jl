using QPJuMPControl
using Test
using RigidBodyDynamics
using JuMP 
using OSQP 
using Random 
using LinearAlgebra
using StaticArrays

const QPC = QPJuMPControl 
const RBD = RigidBodyDynamics

include("tasks.jl")
include("trajectories.jl")

