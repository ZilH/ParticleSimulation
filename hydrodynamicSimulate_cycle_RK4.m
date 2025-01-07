function hydrodynamicSimulate_cycle_RK4(Q, runNum)
% clear
% close all

% Define simulation constants
tic

N_cycles = 100;

makeVideo = true;
continuousVideo = false;
markersize = 45;

isLubrication = true;

isLoadTestXY = true;
savedata = true;

%% Physical parameters

% Modify every experiment
% Q = 7000 * 10^-3 ; %mL/min

% Fix parameters
rho_mix = 1180; %kg/m^3
T = 10.11; %s
D = 6.35; % mm, tube diameter
dp = 0.6; % mm, Particle diameter
A = pi/4 * D^2; % mm^2
mu_mix = 2.38; % Pa*S dynamic viscosity from Snook 2016

v_avg = (Q * 16.67) / A * 10^-3; % m/s
F_p_drag = 3 * pi * mu_mix * (dp * 10^-3) * v_avg;  % N

dt_exp = 0.01;     % [s] Time step for integration (adjust as needed)
dt = dt_exp / ((dp * 10^-3) / v_avg);

half_cycle_iters = round(5 / dt_exp);
iters = 2 * half_cycle_iters * N_cycles + 1;

%% Define simulation constants

% Set the parameters -- all non dimensional
phi = 0.8;
num_particles = 200;
% num_particles = 100;
d = 1;
epsilon = d/2;

% mu = 1;        % Dynamic viscosity of the fluid
% hydro_coeff = 3 * pi * mu * d;  % Hydrodynamic drag coefficient
F0 = 1;        % Magnitude of the contact force (adjust as needed)
Vc_magnitude = 1.5; % Flow between 2 plates


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
% particle_y = rand(num_particles, 1) * box_height;

if isLoadTestXY
    load('pack_xy_jiter.mat');
end

R = box_height / 2;

% Plot the particles as circles
init_x = particle_x;
init_y = particle_y;

%% Particle motion
active = zeros(iters,1);
plug_areas = zeros(iters,1);
msd_x = zeros(iters,1);
msd_y = zeros(iters,1);
n_contact = zeros(iters,1);
n_lubricate = zeros(iters,1);

msd_x_cycle = zeros(N_cycles,1);
msd_y_cycle = zeros(N_cycles,1);
particle_x_cycle = cell(N_cycles,1);
particle_y_cycle = cell(N_cycles,1);
n_contact_cycle = zeros(N_cycles,1);
n_lubricate_cycle = zeros(N_cycles,1);


% Set up video writer
if makeVideo
    if isLubrication
        vidName = sprintf('2DSimulation_Q%03d_3rdLaw_nondim_lub',Q*1000);
    else
        vidName = sprintf('2DSimulation_Q%03d_3rdLaw_nondim',Q*1000);
    end
    vidObj = VideoWriter(vidName);
    vidObj.FrameRate = 10;  % Set the frame rate
    open(vidObj);
end

% if is_preshar
%     for i = 1:num_particles
%         particle_x(i) = particle_x(i) + Vc * (1 - (particle_y(i) - R)^2/R^2);
%     end
%     init_x = particle_x;
%     init_y = particle_y;
% end

% Prepare parameters for compute_u
params.num_particles = num_particles;
params.R = R;
params.d = d;
params.epsilon = epsilon;
params.isLubrication = isLubrication;


for it = 1:iters
    if mod(it, 2*half_cycle_iters) == 1  % Start of a new cycle
        % Record the initial positions for the current cycle
        init_x_cycle = particle_x;
        init_y_cycle = particle_y;
    end
    
    if mod(floor((it - 1) / half_cycle_iters), 2) == 0
        Vc = Vc_magnitude;  % Positive for iterations
    else
        Vc = -Vc_magnitude;  % Negative for iterations
    end
    
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
            xlim([-box_width 3 * box_width]);
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
    params.Vc = Vc;
    x_n = particle_x;
    y_n = particle_y;
    
    %% Eulerian
    %     u = compute_u(x_n, y_n, params);
    %
    %     N = num_particles;
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
    %         % Apply boundary conditions...
    %         if is_upbotwallHard
    %             if particle_y(i) >= box_height - d/2
    %                 particle_y(i) = box_height - d/2;
    %             elseif particle_y(i) < d/2
    %                 particle_y(i) = d/2;
    %             end
    %         end
    %     end
    
    %% RK 4
    % RK4 steps
    [u1, num_contact_k1, num_lubrication_k1] = compute_u(x_n, y_n, params);
    k1_x = u1(1:2:end);
    k1_y = u1(2:2:end);
    
    x_k2 = x_n + (dt/2) * k1_x;
    y_k2 = y_n + (dt/2) * k1_y;
    [u2, num_contact_k2, num_lubrication_k2] = compute_u(x_k2, y_k2, params);
    k2_x = u2(1:2:end);
    k2_y = u2(2:2:end);
    
    x_k3 = x_n + (dt/2) * k2_x;
    y_k3 = y_n + (dt/2) * k2_y;
    [u3, num_contact_k3, num_lubrication_k3] = compute_u(x_k3, y_k3, params);
    k3_x = u3(1:2:end);
    k3_y = u3(2:2:end);
    
    x_k4 = x_n + dt * k3_x;
    y_k4 = y_n + dt * k3_y;
    [u4, num_contact_k4, num_lubrication_k4] = compute_u(x_k4, y_k4, params);
    k4_x = u4(1:2:end);
    k4_y = u4(2:2:end);
    particle_x = x_n + (dt/6) * (k1_x + 2*k2_x + 2*k3_x + k4_x);
    particle_y = y_n + (dt/6) * (k1_y + 2*k2_y + 2*k3_y + k4_y);
    
    % Vectorized boundary enforcement
    particle_y = min(max(particle_y, d/2), box_height - d/2);
    
    total_contact_pairs = (num_contact_k1 + num_contact_k2 + num_contact_k3 + num_contact_k4) / 4;
    total_lubrication_pairs = (num_lubrication_k1 + num_lubrication_k2 + num_lubrication_k3 + num_lubrication_k4) / 4;

    %     msd_x(it) = msd_x(it) + (particle_x - init_x)^2;
    %     msd_y(it) = msd_y(it) + (particle_y - init_y)^2;
    %
    %     msd_x(it) = msd_x(it) / num_particles;
    %     msd_y(it) = msd_y(it) / num_particles;
    
    if mod(it, 2*half_cycle_iters) == 0 || it == 1
        [~, ~, ~, ~, p_area] = cal_grad_concentration(particle_x, particle_y, grid_size, sigma, Vc, box_height / 2, is_shearTimesConcent);
        plug_areas(it) = p_area;
    end
    
    % At the end of the cycle (right before starting a new cycle), calculate MSD
    if mod(it, 2*half_cycle_iters) == 0  % End of a cycle
        msd_x_cycle(it / (2*half_cycle_iters)) = mean((particle_x - init_x_cycle).^2);
        msd_y_cycle(it / (2*half_cycle_iters)) = mean((particle_y - init_y_cycle).^2);
        particle_x_cycle{it/(2*half_cycle_iters)} = particle_x;
        particle_y_cycle{it/(2*half_cycle_iters)} = particle_y;
        n_contact_cycle(it/(2*half_cycle_iters)) = total_contact_pairs;
        n_lubricate_cycle(it/(2*half_cycle_iters)) = total_lubrication_pairs;
        toc
        
%         if isLoadTestXY
%             benchmark_data = load('test_cycle1_dt10E-3.mat');
%             bench_x = benchmark_data.particle_x;
%             bench_y = benchmark_data.particle_y;
%             
%             isEqual = isequal(particle_x, bench_x) & isequal(particle_y, bench_y)
%         end
    end
    
    
end
% Close video writer
if makeVideo
    close(vidObj);
end

%%

%
% for n = 1:N_cycles
%     startIdx = 1 + (n - 1) * 2 * half_cycle_iters;
%     endIdx = startIdx + 2 * half_cycle_iters;
%
%     msd_x_cycle(n) = sum(msd_x(startIdx:endIdx));
%     msd_y_cycle(n) = sum(msd_y(startIdx:endIdx));
% end

% plot(msd_x_cycle)
if ~continuousVideo && savedata
    if ~isLubrication
        save(sprintf('./10E-2dt/RK4/PackContactOnly/Q%d_run%d.mat',Q*1000,runNum))
    else
        save(sprintf('./10E-2dt/RK4/LubricationRes/Q%d_run%d.mat',Q*1000, runNum))
    end
end
end