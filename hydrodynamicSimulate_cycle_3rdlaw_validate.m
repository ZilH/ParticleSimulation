clear
close all

% Define simulation constants
N_cycles = 1;

makeVideo = false;
continuousVideo = true;
markersize = 45;

isLubrication = true;

%% Physical parameters

% Modify every experiment
Q = 1000 * 10^-3 ; %mL/min

% Fix parameters
rho_mix = 1180; %kg/m^3
T = 10.11; %s
D = 6.35; % mm, tube diameter
dp = 0.6; % mm, Particle diameter
A = pi/4 * D^2; % mm^2
mu_mix = 2.38; % Pa*S dynamic viscosity from Snook 2016

v_avg = (Q * 16.67) / A * 10^-3; % m/s
F_p_drag = 3 * pi * mu_mix * (dp * 10^-3) * v_avg;  % N

% dt = 0.0001;     % [s] Time step for integration (adjust as needed)
dt = 1 * 10^-6;
% dt = dt / ((dp * 10^-3) / v_avg);

% half_cycle_iters = 250;
% iters = 2 * half_cycle_iters * N_cycles + 1;
iters = 500000000;

%% Define simulation constants

% Set the parameters -- all non dimensional
phi = 0.8;
num_particles = 2;
% num_particles = 100;
d = 1;   
epsilon = d/2;

% mu = 1;        % Dynamic viscosity of the fluid
% hydro_coeff = 3 * pi * mu * d;  % Hydrodynamic drag coefficient
F0 = 0;        % Magnitude of the contact force (adjust as needed)
% Vc_magnitude = 0.01; % Flow between 2 plates


gamma0 = 5.0; % Strain Amplitude. Also the slope of the velocity profile
is_preshar = false;
is_upbotwallHard = true;
is_upbotwallBounce = false;
is_shearTimesConcent = false;
periodicHorizontal = false;

% Define the grid size and sigma
is_concentration_grad = false;
grid_size = [100, 100]; % Example grid size
sigma = 1; % Example standard deviation for the Gaussian

%%
% box_height = sqrt(num_particles * pi / (8 * phi)) * d;
% box_width = box_height * 2;

aspect_ratio = 10.58;   % b/d
box_height = d * aspect_ratio;
box_width = (num_particles * pi * d) / (4 * aspect_ratio * phi);

% Generate random positions for the particles
particle_x = rand(num_particles, 1) * box_width;
particle_y = rand(num_particles, 1) * (box_height - d) + d/2;

particle_x = [0; -2];
particle_y = [0; 0.1];

delta_x = zeros(iters / 1000,1);
delta_y = zeros(iters / 1000,1);
displace_x = zeros(iters / 1000,1);
displace_y = zeros(iters / 1000,1);

R = box_height / 2;

%% Particle motion
% active = zeros(iters,1);
% plug_areas = zeros(iters,1);
% msd_x = zeros(iters,1);
% msd_y = zeros(iters,1);

% msd_x_cycle = zeros(N_cycles,1);
% msd_y_cycle = zeros(N_cycles,1);


% Set up video writer
if makeVideo
    vidName = sprintf('Test_2DSimulation_Q%03d_3rdLaw_nondim_lub_validate',Q*1000);
    vidObj = VideoWriter(vidName);
    vidObj.FrameRate = 15;  % Set the frame rate
    open(vidObj);
end

% if is_preshar
%     for i = 1:num_particles
%         particle_x(i) = particle_x(i) + Vc * (1 - (particle_y(i) - R)^2/R^2);
%     end
%     init_x = particle_x;
%     init_y = particle_y;
% end
Vc_magnitude = 1 / R;

for it = 1:iters
%     if mod(floor((it - 1) / half_cycle_iters), 2) == 0
%         Vc = Vc_magnitude;  % Positive for iterations
%     else
%         Vc = -Vc_magnitude;  % Negative for iterations
%     end
    Vc = Vc_magnitude;
    
    particle_x_forward = 0 * particle_x;
    
    for i = 1:num_particles
        particle_x_forward(i) = particle_x(i) + Vc * (1 - (particle_y(i) - R)^2/R^2) * dt;
    end
    
    % Check Collision before shear
%     D1 = pdist([particle_x,particle_y]);
%     Dsq = squareform(D1);
%     ind = find(Dsq < d & Dsq ~= 0);
%     % Convert linear indices to row and column indices
%     [nrows, ncols] = size(Dsq);
%     [row_ind_0, col_ind] = ind2sub([nrows, ncols], ind);
%     
%     % Check Collision after shear
%     D1 = pdist([particle_x_forward,particle_y]);
%     Dsq = squareform(D1);
%     ind = find(Dsq < d & Dsq ~= 0);
%     % Convert linear indices to row and column indices
%     [nrows, ncols] = size(Dsq);
%     [row_ind_1, col_ind] = ind2sub([nrows, ncols], ind);
%     
%     
%     row_ind = unique([row_ind_0;row_ind_1]);
%     active(it) = length(row_ind);
    
%     figure
%     set(gcf, 'Position',  [100, 100, 1200, 400])
%     axis equal
%     scatter(particle_x, particle_y, markersize);
%     xlim([0 2*box_width]);
%     ylim([0 box_height]);
%     xlabel('X position (units)');
%     ylabel('Y position (units)');
%     title('Original Location');
%     hold on
%     scatter(particle_x_forward, particle_y, markersize);
%     scatter(particle_x_forward(row_ind, :), particle_y(row_ind, :), markersize);
    
    if makeVideo
        if continuousVideo
            condition = true;
        else
            condition = mod(it, 2*half_cycle_iters) == 0 || it == 1;
        end

            
        if condition
            scatter(particle_x, particle_y, markersize);
            hold on
            viscircles([particle_x, particle_y],ones(num_particles,1).*0.5,'Color','b');

            set(gcf, 'Position',  [100, 100, 1200, 400])
            axis equal
            % scatter(particle_x, particle_y, markersize);
    %         xlim([0 2*box_width]);
            xlim([-15 15]);
            ylim([0 box_height]);
            xlabel('X position (units)');
            ylabel('Y position (units)');
%             title(sprintf('2D Simulation, %c=%.1f, Cycle #%02d',947,gamma0,it));
            title(sprintf('2D Simulation, Vc=%.1f, Steps #%02d',Vc_magnitude, it));

            % Plot the active particles
%             scatter(particle_x(row_ind, :), particle_y(row_ind, :), markersize,'r');
            hold off

            currFrame = getframe(gcf);
            writeVideo(vidObj, currFrame);
        end
        
%         clf(f)
    end
    
 %% Contact Force
 
%     D1 = pdist([particle_x,particle_y]);
%     Dsq = squareform(D1);
%     
%     
    F_c_x = zeros(num_particles, 1);
    F_c_y = zeros(num_particles, 1);
%     % Loop over each unique pair of particles to avoid double counting
%     for i = 1:num_particles-1
%         for j = i+1:num_particles
%             x_ij = Dsq(i, j);
%             h_ij = x_ij - d;
% 
%             % Check for overlap
%             if h_ij < 0
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
%                 F_contact = F0;
% %                 F_contact = abs(R-particle_y(i)) * 6 / R;
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
% 
%     total_F_c_x(it) = sum(abs(F_c_x));
%     total_F_c_y(it) = sum(abs(F_c_y));
    
%     fprintf('Total contact force in x-direction: %e\n', total_F_c_x);
%     fprintf('Total contact force in y-direction: %e\n', total_F_c_y);

%% Lubrication Force
    if isLubrication
        N = num_particles;
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
            K = 1;

            % Diagonal terms (drag)
            A(idx_i_x, idx_i_x) = K;
            A(idx_i_y, idx_i_y) = K;

            % Right-hand side (drag against background flow and contact forces)
            b(idx_i_x) = K * u_inf_i_x - F_c_x(i);
            b(idx_i_y) = K * u_inf_i_y - F_c_y(i);
        end
        
        % Loop over unique particle pairs for lubrication forces
        for i = 1:N-1
            for j = i+1:N
                dx = particle_x(j) - particle_x(i);
                dy = particle_y(j) - particle_y(i);
                x_ij = sqrt(dx^2 + dy^2);

                h_ij = x_ij - d;

                if h_ij > 0 && h_ij <= epsilon % epsilon is a small cutoff value
                    % Lubrication coefficient 
                    L_ij = 1 / (8 * h_ij);

                    % Unit vector
                    n_ij_x = dx / x_ij;
                    n_ij_y = dy / x_ij;

                    idx_i_x = 2*i - 1;
                    idx_i_y = 2*i;
                    idx_j_x = 2*j - 1;
                    idx_j_y = 2*j;

                    % Lubrication terms
                    % For particle i
                    A(idx_i_x, idx_i_x) = A(idx_i_x, idx_i_x) + L_ij * n_ij_x^2;
                    A(idx_i_x, idx_i_y) = A(idx_i_x, idx_i_y) + L_ij * n_ij_x * n_ij_y;
                    A(idx_i_y, idx_i_x) = A(idx_i_y, idx_i_x) + L_ij * n_ij_x * n_ij_y;
                    A(idx_i_y, idx_i_y) = A(idx_i_y, idx_i_y) + L_ij * n_ij_y^2;

                    % For particle j
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
        u = A \ b;

        % Extract velocities
        for i = 1:N
            idx_i_x = 2*i - 1;
            idx_i_y = 2*i;

            u_i_x = u(idx_i_x);
            u_i_y = u(idx_i_y);

            % Update particle positions
            particle_x(i) = particle_x(i) + u_i_x * dt;
            particle_y(i) = particle_y(i) + u_i_y * dt;
            
            if mod(it, 1000) == 0
                displace_x(it / 1000) = u_i_x * dt;
                displace_y(it / 1000) = u_i_y * dt;
            end

%             displace_x(it) = u_i_x * dt;
%             displace_y(it) = u_i_y * dt;

            % Apply boundary conditions...
%             if is_upbotwallHard
%                 if particle_y(i) >= box_height - d/2
%                     particle_y(i) = box_height - d/2;
%                 elseif particle_y(i) < d/2
%                     particle_y(i) = d/2;
%                 end
%             end
%             if particle_y(i) < 0
%                 particle_y(i) = 0;
%             end
        end
               
    else
        for i = 1:num_particles
            % Compute background flow velocity at particle i using parabolic profile
            %         u_inf_i_x = Vc * (1 - ((particle_y(i) - R) / R)^2);
            u_inf_i_x = Vc * particle_y(i) * (2 * R - particle_y(i)) / R^2;
            u_inf_i_y = 0; % No flow in y-direction
            
            % Calculate particle velocity from force balance
            u_i_x = u_inf_i_x - F_c_x(i);
            u_i_y = u_inf_i_y - F_c_y(i);
            
            
            % Update particle positions
            particle_x(i) = particle_x(i) + u_i_x * dt;
            particle_y(i) = particle_y(i) + u_i_y * dt;
            
            %         msd_x(it) = msd_x(it) + (u_i_x * dt)^2;
            %         msd_y(it) = msd_y(it) + (u_i_y * dt)^2;
            
            % No penetration boundary condition:
            if is_upbotwallHard
                if particle_y(i) >= box_height - d/2
                    particle_y(i) = box_height - d/2;
                elseif particle_y(i) < d/2
                    particle_y(i) = d/2;
                end
            end
            
            % Bounding boundary condition :
            %         if is_upbotwallBounce
            %             if particle_y(row_ind(idx)) >= box_height - d/2
            %                 particle_y(row_ind(idx)) = box_height - d/2 - epsilon * rand;
            %             elseif particle_y(row_ind(idx)) < d/2
            %                 particle_y(row_ind(idx)) = d/2 + epsilon * rand;
            %             end
            %         end
        end        
    end
    
    if mod(it, 1000) == 0
        delta_x(it / 1000) = particle_x(1) - particle_x(2);
        delta_y(it / 1000) = particle_y(1) - particle_y(2);
    end
end
% Close video writer
if makeVideo
    close(vidObj);
end

plot(delta_x, delta_y)