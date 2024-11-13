clear;
close all;

% Specify the folder where the .mat files are located
folder = './10E-2dt/RK4/LubricationRes';  % Change this to your actual folder path

% Get a list of all .mat files in the folder
matFiles = dir(fullfile(folder, '*.mat'));

% Define number of groups and files per group
numGroups = 8;
filesPerGroup = 10;

% Initialize arrays to store means and standard deviations
meanPlugAreas = cell(1, numGroups);
stdPlugAreas = cell(1, numGroups);
meanMSDs = cell(1, numGroups);
stdMSDs = cell(1, numGroups);

% Loop over each group of 10 files
for groupIdx = 1:numGroups
    plugAreaData = [];
    msdData = [];

    % Loop over the 10 files in each group
    for k = 1:filesPerGroup
        fileIdx = (groupIdx - 1) * filesPerGroup + k;
        
        % Get the full file name
        fileName = fullfile(folder, matFiles(fileIdx).name);
        
        % Load the .mat file
        data = load(fileName);
        
        % Check if 'plug_areas' and 'msd_x_cycle' exist in the file
        if isfield(data, 'plug_areas') && isfield(data, 'msd_x_cycle')
            % Store 'plug_areas' and 'msd_x_cycle' in temporary arrays
            tmp = 1000:1000:length(data.plug_areas);
            index = [1 tmp];
            plugAreaData = [plugAreaData; data.plug_areas(index)'];
            msdData = [msdData; data.msd_x_cycle'];
        else
            warning(['Variables not found in ', matFiles(fileIdx).name]);
        end
    end
    
    % Calculate mean and standard deviation for each group of 10 files
    meanPlugAreas{groupIdx} = mean(plugAreaData, 1);
    stdPlugAreas{groupIdx} = std(plugAreaData, 0, 1);
    meanMSDs{groupIdx} = mean(msdData, 1);
    stdMSDs{groupIdx} = std(msdData, 0, 1);
end

%%
% Plotting the plug area averages with error bars
figure(1);
hold on;
for groupIdx = 1:numGroups
    % Define x-axis points for each group based on the index of plug_areas
    x = 1:length(meanPlugAreas{groupIdx});
    
    % Plot mean with error bars
    errorbar(x, meanPlugAreas{groupIdx}, stdPlugAreas{groupIdx}, ...
        'DisplayName', ['Group ', num2str(groupIdx)]);
end
xlabel('Index');
ylabel('Plug Area');
title('Average Plug Area with Error Bars');
legend show;
hold off;

% Plotting the MSD_y averages with error bars on a logarithmic scale
figure(2);
hold on;
for groupIdx = 7:numGroups
    % Define x-axis points for each group based on the length of msd_x_cycle
    x = 1:length(meanMSDs{groupIdx});
    
    % Plot mean with error bars
    errorbar(x, meanMSDs{groupIdx}, stdMSDs{groupIdx}, ...
        'DisplayName', ['Group ', num2str(groupIdx)]);
    
end
set(gca, 'YScale', 'log');
xlabel('Index');
ylabel('MSD_y');
title('Average MSD_y with Error Bars');
% ylim([10^-5, 100]);
legend show;
hold off;
