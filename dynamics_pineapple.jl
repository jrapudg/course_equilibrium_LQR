###########################
# Contact Frame Utilities #
###########################
function get_wheel_contour(rot_normal)
    return atan(rot_normal[1], rot_normal[3])
end

function wheels_velocity_contribution(model, x, normals)
    # Define wheel parameters (better to do globally or as part of the model)
    ## Wheel radius
    # if model.nu == 6
    #     wheel_r = 0.07;  # Wheel radius for 6 DoF model
    # elseif model.nu == 8
    #     wheel_r = 0.075;  # Wheel radius for 8 DoF model
    # end 

    wheel_r = 0.075;  # Wheel radius
    ## Floor normals (assuming the wheels are on a flat surface)
    lw_normal = normals[1:3]
    rw_normal = normals[4:6]

    # println("Normals: $normals")
    # println("Left wheel normal: $lw_normal")
    # println("Right wheel normal: $rw_normal")
    ## Define the wheel axis of rotation in the local frame
    wheel_local_rot_axis = [0, 1, 0]  # Right wheel axis of rotation

    # Get the wheel rotations
    wheels_rotations = kinematics_rotation(model, x)
    #@show wheels_rotations
    # wheels_rotations = [0.7071067811865475 -0.0 -0.7071067811865476; 0.0 1.0 -0.0; 0.7071067811865476 0.0 0.7071067811865475; 0.7071067811865475 -0.0 -0.7071067811865476; 0.0 1.0 -0.0; 0.7071067811865476 0.0 0.7071067811865475]

    # @show size(wheels_rotations)
    wheels_rotations = reshape(vec(wheels_rotations), 6, 3) 
    # @show size(wheels_rotations)
    ## Rotation matrices from wheel to world frame
    lw_rotation = wheels_rotations[1:3, :]  # Left wheel rotation matrix
    rw_rotation = wheels_rotations[4:6, :]  # Right wheel rotation matrix

    # Floor normals in wheel frame
    lw_rot_normal = lw_rotation' * lw_normal
    rw_rot_normal = rw_rotation' * rw_normal

    # @show size(lw_rotation) size(rw_rotation)
    # @show size(lw_rot_normal) size(rw_rot_normal)

    # Compute vector from wheel center to contact point in the world frame
    ## Get the contact contour for each wheel
    lw_contour = get_wheel_contour(lw_rot_normal)
    rw_contour = get_wheel_contour(rw_rot_normal)

    ## Get the local contact vectors
    local_lw_radius_contact = [-wheel_r * sin(lw_contour), 0, -wheel_r * cos(lw_contour)]
    local_rw_radius_contact = [-wheel_r * sin(rw_contour), 0, -wheel_r * cos(rw_contour)]
    
    # Rotate the local contact vectors to the world frame
    ## Get the wheel centers in the world frame
    ## Compute the contact points in the world frame
    lw_radius_contact = lw_rotation * local_lw_radius_contact;
    rw_radius_contact = rw_rotation * local_rw_radius_contact;
    
    # Build wheel angular velocities Jacobian
    nj = Int(model.nu/2) - 1
    # println("nj: $nj")
    Jiw = [zeros(2,model.nq) zeros(2, 6)  zeros(2, nj)  [1, 0]  zeros(2, nj) [0, 1]];
    
    # Add rotation axis in world frame, which is assumed to be the z-axis
    rot_axis = wheels_rotations*wheel_local_rot_axis; # Rotation axis for the wheel contact plane
    P = [rot_axis[1:3]  zeros(3);
         zeros(3)   rot_axis[4:6]]
    
    # Compute the local angular velocities magnitude
    ω_local = Jiw*x;
    # Compute the angular velocities with rotation axis
    ω_ic = P*ω_local;

    # Compute the contact linear velocities
    ω_lic_cross = skew(ω_ic[1:3]);
    ω_ric_cross = skew(ω_ic[4:6]);
    vic_right = ω_ric_cross * rw_radius_contact;
    vic_left = ω_lic_cross * lw_radius_contact;

    # J_dot is [vic_right; vic_left] since this depends on wheel accelerations 
    # after derivativing it with respect to time
    # That is: [ω̇_IW] × r_WC(σ) in World frame
    
    # Lets compute J that depends on v in acceleration level constraint
    # Time derivative of the velocity contribution due to rotation:
    # [ω_IW]^2 × r_WC(σ)+ [ω_IW]× * R_IW * d/dσ (W_r_WC(σ))
    lw_radius_contact_local_dot = [-wheel_r * cos(lw_contour), 0, wheel_r * sin(lw_contour)];
    rw_radius_contact_local_dot = [-wheel_r * cos(rw_contour), 0, wheel_r * sin(rw_contour)];

    # # Compute the contact contour detivatives
    lw_contour_dot_normal = 1/(lw_rot_normal[1]^2 + lw_rot_normal[3]^2) * 
                            [lw_rot_normal[3], 0, -lw_rot_normal[1]]; # Derivative of the contact contour for the left wheel
    rw_contour_dot_normal = 1/(rw_rot_normal[1]^2 + rw_rot_normal[3]^2) * 
                            [rw_rot_normal[3], 0, -rw_rot_normal[1]]; # Derivative of the contact contour for the right wheel
    
    l_normal_dot = skew(wheel_local_rot_axis*ω_local[1])*lw_rot_normal;
    r_normal_dot = skew(wheel_local_rot_axis*ω_local[2])*rw_rot_normal;

    lw_contour_dot = lw_contour_dot_normal' * l_normal_dot;
    rw_contour_dot = rw_contour_dot_normal' * r_normal_dot;

    aic_left = ω_lic_cross * ω_lic_cross * lw_radius_contact + ω_lic_cross * lw_rotation * lw_radius_contact_local_dot * lw_contour_dot;
    aic_right = ω_ric_cross * ω_ric_cross * rw_radius_contact + ω_ric_cross * rw_rotation * rw_radius_contact_local_dot * rw_contour_dot;

    Jic = [vic_left; vic_right; aic_left; aic_right];
    return Jic
end

function jacobian_wheels_velocity_contribution(model, x, normals = [[0,0,1]; [0,0,1]])
    ForwardDiff.jacobian(x -> wheels_velocity_contribution(model, x, normals), x)
end

"""
    ẋ, λ = pinned_continuous_dynamics(model::Pineapple, x, u; K_pd = zeros(model.nu, model.nx))

Only supported for models where model.kinematics_ori is :Quaternion
"""
function pinned_continuous_dynamics(model::Pineapple, x, u;  K_pd = zeros(model.nu, model.nx))
    @show model.kinematics_ori
    @assert model.kinematics_ori == :Quaternion "model needs to have kinematics orientation with quaternions, Pineapple(kinematics_ori = :Quaterion)"
    # Get dynamics matrices
    M = M_func(model, x)
    C = C_func(model, x)
    B = B_func(model)

    # Map v into q̇ using velocity kinematics
    v = x[model.nq + 1:end]
    E = velocity_kinematics(model, x)
    q̇ = E*v

    # 
    pos_indices = [1:3; 7 .+ (1:3)]
    # Get constraint derivatives and split into d/dq and d/dv terms
    JJ_wheel_center = kinematics_velocity_jacobian(model, x)[pos_indices, :] # This is d/dx(J(q)E(q)v)
    Jdot_wheel_center  = JJ_wheel_center[:, 1:model.nq] * E
    J_wheel_center  = JJ_wheel_center[:, model.nq + 1:end]

    JJ_wheel_contact = jacobian_wheels_velocity_contribution(model, x)
    Jdot_wheel_contact = JJ_wheel_contact[7:end, model.nq+1:end]
    J_wheel_contact = JJ_wheel_contact[1:6, model.nq+1:end]
    
    Jdot = Jdot_wheel_center + Jdot_wheel_contact
    J = J_wheel_center + J_wheel_contact

    JA = [J[[1,3],:]; J[[4,6],:]]  # X-Z Jacobian
    JAdot = [Jdot[[1,3],:]; Jdot[[4,6],:]]  # X-Z Jacobian derivative
    
    JF = J[[2,5],:]  # Y Jacobian
    # JFdot = Jdot[[2,5],:]  # Y Jacobian derivative
    JFv = JF*v
    CF = -model.μ*[zeros(2,1) [tanh(JFv[1]), 0] zeros(2,1) [0, tanh(JFv[2])]]  # Y Coriolis + gravity vectorm
    # Combine all jacobians

    n_c = size(JA, 1)
    # Build and solve KKT system for v̇ and λ
    # M vdot - J' lambda = B u - C
    # J vdot = - Jdot v   
    kkt_matrix = [M (JA + CF'JF)'; JA zeros(n_c, n_c)]
    kkt_rhs_top = B*u - C - B*K_pd*x  
    kkt_rhs_bottom = -JAdot*v
    kkt_rhs = [kkt_rhs_top; kkt_rhs_bottom]
    res = kkt_matrix \ kkt_rhs

    # Extract v̇ and λ
    v̇ = res[1:model.nv]
    λ = res[model.nv + 1:end]
    return [q̇; v̇], λ
end

function linesearch(z::Vector, Δz::Vector, merit_fx::Function;
                    max_ls_iters = 10)::Float64 # optional argument with a default
    
    # return maximum α≤1 such that merit_fx(z + α*Δz) < merit_fx(z)
    # with a backtracking linesearch (α = α/2 after each iteration)

    # NOTE: DO NOT USE A WHILE LOOP 
    ϕ0 = merit_fx(z)
    α = 1.0 
    for i = 1:max_ls_iters
        if merit_fx(z + α*Δz) < ϕ0 
            return α
        else
            α = α/2
        end
    end
    return 1.0 
    # error("linesearch failed")
end

function newtons_method(z0::Vector, res_fx::Function, res_jac_fx::Function, merit_fx::Function;
                        tol = 1e-10, max_iters = 50, verbose = false)::Vector{Vector{Float64}}
    # - z0, initial guess 
    # - res_fx, residual function 
    # - res_jac_fx, Jacobian of residual function wrt z 
    # - merit_fx, merit function for use in linesearch 
    
    # optional arguments 
    # - tol, tolerance for convergence. Return when norm(residual)<tol 
    # - max iter, max # of iterations 
    # - verbose, bool telling the function to output information at each iteration
    
    # return a vector of vectors containing the iterates 
    # the last vector in this vector of vectors should be the approx. solution 
    
    # return the history of guesses as a vector
    Z = [zeros(length(z0)) for i = 1:max_iters]
    Z[1] = z0 
    
    
    for i = 1:(max_iters - 1)
        
        # evaluate current residual 
        r = res_fx(Z[i])
        norm_r = norm(r)
        if verbose 
            print("iter: $i    |r|: $norm_r   ")
        end
        
        # check convergence with norm of residual < tol 
        if norm_r < tol
            return Z[1:i]
        end
        
        # caculate Newton step (don't forget the negative sign)
        J = res_jac_fx(Z[i])
        Δz = -J\r 
        
        # linesearch and update z 
        α = linesearch(Z[i], Δz, merit_fx)
        Z[i+1] = Z[i] + α*Δz
        if verbose
            print("α: $α \n")
        end
        
    end
    error("Newton's method did not converge")
end