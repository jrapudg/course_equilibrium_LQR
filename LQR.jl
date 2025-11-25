using Pkg
Pkg.develop(path="dependencies/PinnZoo")

Pkg.add("MeshCatMechanisms")
Pkg.add("RigidBodyDynamics")

using PinnZoo
using MeshCat
using MeshCatMechanisms
import RigidBodyDynamics: MechanismState, joints, parse_urdf, QuaternionFloating
import ForwardDiff
using LinearAlgebra

Pkg.add(Pkg.PackageSpec(name = "ForwardDiff", version = v"0.10.39"))

include("dynamics_pineapple.jl")

model = Pineapple();

robot = parse_urdf(model.urdf_path, floating=true);
vis = MechanismVisualizer(robot, URDFVisuals(model.urdf_path, package_path=[dirname(model.urdf_path)]));
render(vis);

let
    state = MechanismState(robot);
    # Retrieve configuration and velocity names
    config_names = [Symbol(joints(robot)[id].name) for id in state.q_index_to_joint_id];
    vel_names = [Symbol(joints(robot)[id].name) for id in state.v_index_to_joint_id];
    # Adjust for floating base
    for joint in joints(robot)
        if typeof(joint.joint_type) <: QuaternionFloating
            config_names[state.qranges[joint]] = [:q_w, :q_x, :q_y, :q_z, :x, :y, :z];
            vel_names[state.vranges[joint]] = [:ang_v_x, :ang_v_y, :ang_v_z, :lin_v_x, :lin_v_y, :lin_v_z];
        end;
    end;
    # Define torque names and generate conversions  
    torque_names = [name for name in vel_names if name in model.orders[:nominal].torque_names];
    model.orders[:rigidBodyDynamics] = StateOrder(config_names, vel_names, torque_names);
    generate_conversions(model.orders, model.conversions);
end;

x = init_state(model);

set_configuration!(vis, change_order(model, x[1:model.nq], :nominal, :rigidBodyDynamics))

foot_locs = kinematics(model, x)

wheel_radius = 0.075;

x_guess = init_state(model);
x_guess[3] = 0.1;

# indexing  
nλ = length(model.kinematics_bodies)*3;
const idx_x = 1:model.nx;
const idx_u = (model.nx+1):(model.nx+model.nu);
const idx_λ = (model.nx+model.nu+1):(model.nx+model.nu+nλ);
nc = model.nv + kinematics_size(model, Val(:None)) + model.nv + 2;
const idx_μ = (model.nx+model.nu+nλ+1):(model.nx+model.nu+nλ+nc);
const idx_kin_pos =  [1:3; 7 .+ (1:3)]; # Indices of contact positions in kinematics function with quaternions

function pineapple_cost(y::Vector)
    T = eltype(y)
    # cost function 
    @assert length(y) == model.nx+model.nu+nλ
    x = y[idx_x]
    u = y[idx_u]
    λ = y[idx_λ]
    
    # return cost 
    err = state_error(model, x, x_guess)
    Q = diagm([1e2*ones(T, 3); 1e3*ones(T, 3); 1e1ones(T, model.nu); ones(T, model.nv)])
    return 0.5*err'*Q*err + 0.5*1e-3*norm(u)^2+ 0.5*1e-3*norm(λ)^2
end

function pineapple_constraint(y::Vector)::Vector
    # constraint function 
    @assert length(y) == model.nx+model.nu+nλ
    x = y[idx_x]
    u = y[idx_u]
    λ = y[idx_λ]

    λ_full = zeros(eltype(y), kinematics_size(model))
    λ_full[idx_kin_pos] = λ
    
    # return constraint
    B = B_func(model)
    return [
        forward_dynamics(model, x, B*u + kinematics_jacobianTvp(model, x, λ_full)) # dynamics equal to zero (no change in state)
        kinematics(model, x)[1:2] - kinematics(model, x_guess)[1:2] # does not move left wheel position in x and y
        kinematics(model, x)[7 .+ (1:2)] - kinematics(model, x_guess)[7 .+ (1:2)] # does not move right wheel position in x and y
        kinematics(model, x)[3] - wheel_radius # make the left wheel to be on the floor
        kinematics(model, x)[7 + 3] - wheel_radius # make the right wheel to be on the floor
        x[model.nq + 1:end] # velocity equal to zero
        x[4:7]'*x[4:7] - 1 # make the attitude of the body to be a unit quaternion
        x[3]-x_guess[3] # make the height of the body to be the same as the initial height
    ]
end

function pineapple_kkt(z::Vector)::Vector
    @assert length(z) == model.nx+model.nu+nλ+nc
    x = z[idx_x]
    u = z[idx_u]
    λ = z[idx_λ]
    μ = z[idx_μ]

    y = [x;u;λ]
    
    # return the KKT conditions 
    return [
        ForwardDiff.gradient(pineapple_cost,y) + ForwardDiff.jacobian(pineapple_constraint, y)'*μ;
        pineapple_constraint(y)
           ]       
end

function pineapple_kkt_jac(z::Vector)::Matrix
    @assert length(z) == model.nx+model.nu+nλ+nc
    x = z[idx_x]
    u = z[idx_u]
    λ = z[idx_λ]

    y = [x;u;λ]
    A = ForwardDiff.jacobian(pineapple_constraint, y)
    H = ForwardDiff.hessian(pineapple_cost, y)
    ρ = 1e-4
    return [(H + ρ*I) A'; A -ρ*I] 
end

function pineapple_merit(z)
    # merit function for the quadruped problem 
    @assert length(z) == model.nx+model.nu+nλ+nc
    r = pineapple_kkt(z)
    return norm(r[1:(model.nx+model.nu+nλ)]) + 1e1*norm(r[(model.nx+model.nu+nλ+1):end])
end

z0 = [x_guess; zeros(model.nu); zeros(nλ); zeros(nc)];

Z = newtons_method(z0, pineapple_kkt, pineapple_kkt_jac, pineapple_merit; tol = 1e-9, verbose = true, max_iters = 200);

x0, u0, λ0 = Z[end][idx_x], Z[end][idx_u], Z[end][idx_λ];

set_configuration!(vis, change_order(model, x0[1:model.nq], :nominal, :rigidBodyDynamics))