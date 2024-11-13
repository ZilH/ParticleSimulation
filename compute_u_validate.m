function u = compute_u_validate(particle_x, particle_y, params)
    % Extract parameters
    N = params.num_particles;
    num_particles = params.num_particles;
    Vc = params.Vc;
    R = params.R;
    d = params.d;
    epsilon = params.epsilon;
    isLubrication = params.isLubrication;
    
%     D1 = pdist([particle_x,particle_y]);
%     Dsq = squareform(D1);

    F_c_x = zeros(N, 1);
    F_c_y = zeros(N, 1);
    % Include any other parameters or forces needed
    
%     for i = 1:num_particles-1
%         for j = i+1:num_particles
%             x_ij = Dsq(i, j);
%             h_ij = x_ij - d;
% 
%             % Check for overlap
%             if h_ij < 2 * 10^-6 * d
%                 % Compute distance components
%                 dx = particle_x(j) - particle_x(i);
%                 dy = particle_y(j) - particle_y(i);
% 
%                 % Avoid division by zero
%                 if x_ij == 0
%                     n_ij_x = 0;
%                     n_ij_y = 0;
%                 else
%                     % Unit vector from particle i to particle j
%                     n_ij_x = dx / x_ij;
%                     n_ij_y = dy / x_ij;
%                 end
% 
%                 % Contact force magnitude
% %                 F_contact = F0;
%                 F_contact = 0;
% 
%                 % Apply forces to particle i
%                 F_c_x(i) = F_c_x(i) + F_contact * n_ij_x;
%                 F_c_y(i) = F_c_y(i) + F_contact * n_ij_y;
% 
%                 % Apply equal and opposite forces to particle j
%                 F_c_x(j) = F_c_x(j) - F_contact * n_ij_x;
%                 F_c_y(j) = F_c_y(j) - F_contact * n_ij_y;
%             end
%         end
%     end


    % Initialize A and b
    A = zeros(2*N, 2*N);
    b = zeros(2*N, 1);

    % Loop over particles to set up drag and contact forces
    for i = 1:N
        idx_i_x = 2*i - 1;
        idx_i_y = 2*i;

        % Background flow velocity at particle i
        u_inf_i_x = Vc * particle_y(i);
        u_inf_i_y = 0;

        % Hydrodynamic drag coefficients
        K = 1;  % Adjust if necessary

        % Diagonal terms (drag)
        A(idx_i_x, idx_i_x) = K;
        A(idx_i_y, idx_i_y) = K;

        % Right-hand side (drag against background flow and contact forces)
        b(idx_i_x) = K * u_inf_i_x - F_c_x(i);
        b(idx_i_y) = K * u_inf_i_y - F_c_y(i);
    end
    
    if isLubrication
    % Loop over unique particle pairs for lubrication forces
        for i = 1:N-1
            for j = i+1:N
                dx = particle_x(j) - particle_x(i);
                dy = particle_y(j) - particle_y(i);
                x_ij = sqrt(dx^2 + dy^2);

                h_ij = x_ij - d;

                if h_ij > 0 && h_ij <= epsilon  % Lubrication cutoff
                    % Lubrication coefficient
                    L_ij = 1 / (8 * h_ij);

                    % Unit vector
                    n_ij_x = dx / x_ij;
                    n_ij_y = dy / x_ij;

                    idx_i_x = 2*i - 1;
                    idx_i_y = 2*i;
                    idx_j_x = 2*j - 1;
                    idx_j_y = 2*j;

                    % Lubrication terms for particle i
                    A(idx_i_x, idx_i_x) = A(idx_i_x, idx_i_x) + L_ij * n_ij_x^2;
                    A(idx_i_x, idx_i_y) = A(idx_i_x, idx_i_y) + L_ij * n_ij_x * n_ij_y;
                    A(idx_i_y, idx_i_x) = A(idx_i_y, idx_i_x) + L_ij * n_ij_x * n_ij_y;
                    A(idx_i_y, idx_i_y) = A(idx_i_y, idx_i_y) + L_ij * n_ij_y^2;

                    % Lubrication terms for particle j
                    A(idx_j_x, idx_j_x) = A(idx_j_x, idx_j_x) + L_ij * n_ij_x^2;
                    A(idx_j_x, idx_j_y) = A(idx_j_x, idx_j_y) + L_ij * n_ij_x * n_ij_y;
                    A(idx_j_y, idx_j_x) = A(idx_j_y, idx_j_x) + L_ij * n_ij_x * n_ij_y;
                    A(idx_j_y, idx_j_y) = A(idx_j_y, idx_j_y) + L_ij * n_ij_y^2;

                    % Coupling terms
                    A(idx_i_x, idx_j_x) = A(idx_i_x, idx_j_x) - L_ij * n_ij_x^2;
                    A(idx_i_x, idx_j_y) = A(idx_i_x, idx_j_y) - L_ij * n_ij_x * n_ij_y;
                    A(idx_i_y, idx_j_x) = A(idx_i_y, idx_j_x) - L_ij * n_ij_x * n_ij_y;
                    A(idx_i_y, idx_j_y) = A(idx_i_y, idx_j_y) - L_ij * n_ij_y^2;

                    % Symmetric entries
                    A(idx_j_x, idx_i_x) = A(idx_j_x, idx_i_x) - L_ij * n_ij_x^2;
                    A(idx_j_x, idx_i_y) = A(idx_j_x, idx_i_y) - L_ij * n_ij_x * n_ij_y;
                    A(idx_j_y, idx_i_x) = A(idx_j_y, idx_i_x) - L_ij * n_ij_x * n_ij_y;
                    A(idx_j_y, idx_i_y) = A(idx_j_y, idx_i_y) - L_ij * n_ij_y^2;
                end
            end
        end
    end

    % Solve for velocities
    A = sparse(A);
    u = A \ b;
end
