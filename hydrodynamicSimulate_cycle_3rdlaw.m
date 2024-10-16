clear
close all

makeVideo = true;
continuousVideo = false;
markersize = 45;

% Define simulation constants
N_iters = 30;
half_cycle_iters = 500;
iters = 2 * half_cycle_iters * N_iters + 1;

% Set the parameters
phi = 0.8;
num_particles = 200;
% num_particles = 100;
d = 1;
epsilon = d/2;

mu = 1;        % Dynamic viscosity of the fluid
F0 = 5;        % Magnitude of the contact force (adjust as needed)
dt = 0.01;     % Time step for integration (adjust as needed)
hydro_coeff = 3 * pi * mu * d;  % Hydrodynamic drag coefficient
Vc_magnitude = 10 * d;


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

R = box_height / 2;

% Generate random positions for the particles
particle_x = rand(num_particles, 1) * box_width;
% particle_y = rand(num_particles, 1) * box_height;
particle_y = rand(num_particles, 1) * (box_height - d) + d/2;


% Plot the particles as circles
init_x = particle_x;
init_y = particle_y;



%% Particle motion
active = zeros(iters,1);
plug_areas = zeros(iters,1);
msd_x = zeros(iters,1);
msd_y = zeros(iters,1);

% Set up video writer
if makeVideo
    vidName = sprintf('Test_2DSimulation_Vc%02d_3rdLaw',Vc_magnitude);
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

for it = 1:iters
    if mod(floor((it - 1) / half_cycle_iters), 2) == 0
        Vc = Vc_magnitude;  % Positive for iterations
    else
        Vc = -Vc_magnitude;  % Negative for iterations
    end
    
    particle_x_forward = 0 * particle_x;
    
    for i = 1:num_particles
        particle_x_forward(i) = particle_x(i) + Vc * (1 - (particle_y(i) - R)^2/R^2) * dt;
    end
    
    % Check Collision before shear
    D1 = pdist([particle_x,particle_y]);
    Dsq = squareform(D1);
    ind = find(Dsq < d & Dsq ~= 0);
    % Convert linear indices to row and column indices
    [nrows, ncols] = size(Dsq);
    [row_ind_0, col_ind] = ind2sub([nrows, ncols], ind);
    
    % Check Collision after shear
    D1 = pdist([particle_x_forward,particle_y]);
    Dsq = squareform(D1);
    ind = find(Dsq < d & Dsq ~= 0);
    % Convert linear indices to row and column indices
    [nrows, ncols] = size(Dsq);
    [row_ind_1, col_ind] = ind2sub([nrows, ncols], ind);
    
    
    row_ind = unique([row_ind_0;row_ind_1]);
    active(it) = length(row_ind);
    
    
%     if active(it) == 0
%         break
%     end
    
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
    
    D1 = pdist([particle_x,particle_y]);
    Dsq = squareform(D1);
    
    total_FCX = 0;
    total_FCY = 0;
    
    F_c_x = zeros(num_particles, 1);
    F_c_y = zeros(num_particles, 1);
    % Loop over each unique pair of particles to avoid double counting
    for i = 1:num_particles-1
        for j = i+1:num_particles
            x_ij = Dsq(i, j);
            h_ij = x_ij - d;

            % Check for overlap
            if h_ij < 2 * 10^-6 * d
                % Compute distance components
                dx = particle_x(j) - particle_x(i);
                dy = particle_y(j) - particle_y(i);

                % Avoid division by zero
                if x_ij == 0
                    n_ij_x = 0;
                    n_ij_y = 0;
                else
                    % Unit vector from particle i to particle j
                    n_ij_x = dx / x_ij;
                    n_ij_y = dy / x_ij;
                end

                % Contact force magnitude
                F_contact = F0;

                % Apply forces to particle i
                F_c_x(i) = F_c_x(i) + F_contact * n_ij_x;
                F_c_y(i) = F_c_y(i) + F_contact * n_ij_y;

                % Apply equal and opposite forces to particle j
                F_c_x(j) = F_c_x(j) - F_contact * n_ij_x;
                F_c_y(j) = F_c_y(j) - F_contact * n_ij_y;
            end
        end
    end

%     total_F_c_x = sum(F_c_x);
%     total_F_c_y = sum(F_c_y);
%     fprintf('Total contact force in x-direction: %e\n', total_F_c_x);
%     fprintf('Total contact force in y-direction: %e\n', total_F_c_y);
    
    for i = 1:num_particles
        % Compute background flow velocity at particle i using parabolic profile
        u_inf_i_x = Vc * (1 - ((particle_y(i) - R) / R)^2);
        u_inf_i_y = 0; % No flow in y-direction
        
        % Calculate particle velocity from force balance
        u_i_x = u_inf_i_x - (1 / hydro_coeff) * F_c_x(i);
        u_i_y = u_inf_i_y - (1 / hydro_coeff) * F_c_y(i);
        
        % Update particle positions
        particle_x(i) = particle_x(i) + u_i_x * dt;
        particle_y(i) = particle_y(i) + u_i_y * dt;
        
        msd_x(it) = msd_x(it) + (u_i_x * dt)^2;
        msd_y(it) = msd_y(it) + (u_i_y * dt)^2;
                
        % No penetration boundary condition:
        if is_upbotwallHard
            if particle_y(i) >= box_height - d/2
                particle_y(i) = box_height - d/2;
            elseif particle_y(i) < d/2
                particle_y(i) = d/2;
            end
        end
        
        % Bounding boundary condition :
        if is_upbotwallBounce
            if particle_y(row_ind(idx)) >= box_height - d/2
                particle_y(row_ind(idx)) = box_height - d/2 - epsilon * rand;
            elseif particle_y(row_ind(idx)) < d/2
                particle_y(row_ind(idx)) = d/2 + epsilon * rand;
            end
        end
    end
    
    msd_x(it) = msd_x(it) / num_particles;
    msd_y(it) = msd_y(it) / num_particles;

    
    if mod(it, 2*half_cycle_iters) == 0 || it == 1
        [~, ~, ~, ~, p_area] = cal_grad_concentration(particle_x, particle_y, grid_size, sigma, Vc, box_height / 2, is_shearTimesConcent);
        plug_areas(it) = p_area;
    end
end
% Close video writer
if makeVideo
    close(vidObj);
end

%%
msd_x_cycle = zeros(N_iters,1);
msd_y_cycle = zeros(N_iters,1);

for n = 1:N_iters
    startIdx = 1 + (n - 1) * 2 * half_cycle_iters;
    endIdx = startIdx + 2 * half_cycle_iters;
    
    msd_x_cycle(n) = sum(msd_x(startIdx:endIdx));
    msd_y_cycle(n) = sum(msd_y(startIdx:endIdx));
end

% plot(msd_x_cycle)

% save(sprintf('./NoLubricationRes/Vc%d.mat',Vc_magnitude))
