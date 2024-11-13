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

dt = 10^-1;
iters = 6000;

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

record_interval = 1 / dt;
delta_x = zeros(iters / record_interval,1);
delta_y = zeros(iters / record_interval,1);
displace_x = zeros(iters / record_interval,1);
displace_y = zeros(iters / record_interval,1);

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
    vidName = sprintf('Test_2DSimulation_Q%03d_3rdLaw_nondim_lub_validate_RK4',Q*1000);
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

% Prepare parameters for compute_u
params.num_particles = num_particles;
params.Vc = Vc_magnitude;
params.R = R;
params.d = d;
params.epsilon = epsilon;
params.isLubrication = isLubrication;

for it = 1:iters
    %     if mod(floor((it - 1) / half_cycle_iters), 2) == 0
    %         Vc = Vc_magnitude;  % Positive for iterations
    %     else
    %         Vc = -Vc_magnitude;  % Negative for iterations
    %     end
    Vc = Vc_magnitude;
    
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
    
    
    N = num_particles;
    
    x_n = particle_x;
    y_n = particle_y;
    
 %% Eulerian
 
%     u = compute_u_validate(x_n, y_n, params);
%     
%     % Extract velocities
%     for i = 1:N
%         idx_i_x = 2*i - 1;
%         idx_i_y = 2*i;
%         
%         u_i_x = u(idx_i_x);
%         u_i_y = u(idx_i_y);
%         
%         % Update particle positions
%         particle_x(i) = particle_x(i) + u_i_x * dt;
%         particle_y(i) = particle_y(i) + u_i_y * dt;
%         
%         if mod(it, record_interval) == 0
%             displace_x(it / record_interval) = u_i_x * dt;
%             displace_y(it / record_interval) = u_i_y * dt;
%         end
%     end
 
 %% RK4
    
% RK4 steps
    u1 = compute_u_validate(x_n, y_n, params);
    k1_x = u1(1:2:end);
    k1_y = u1(2:2:end);

    x_k2 = x_n + (dt/2) * k1_x;
    y_k2 = y_n + (dt/2) * k1_y;
    u2 = compute_u_validate(x_k2, y_k2, params);
    k2_x = u2(1:2:end);
    k2_y = u2(2:2:end);

    x_k3 = x_n + (dt/2) * k2_x;
    y_k3 = y_n + (dt/2) * k2_y;
    u3 = compute_u_validate(x_k3, y_k3, params);
    k3_x = u3(1:2:end);
    k3_y = u3(2:2:end);

    x_k4 = x_n + dt * k3_x;
    y_k4 = y_n + dt * k3_y;
    u4 = compute_u_validate(x_k4, y_k4, params);
    k4_x = u4(1:2:end);
    k4_y = u4(2:2:end);

    particle_x = x_n + (dt/6) * (k1_x + 2*k2_x + 2*k3_x + k4_x);
    particle_y = y_n + (dt/6) * (k1_y + 2*k2_y + 2*k3_y + k4_y);

    
    if mod(it, record_interval) == 0
        delta_x(it / record_interval) = particle_x(1) - particle_x(2);
        delta_y(it / record_interval) = particle_y(1) - particle_y(2);
    end
end


% Close video writer
if makeVideo
    close(vidObj);
end

plot(delta_x, delta_y)