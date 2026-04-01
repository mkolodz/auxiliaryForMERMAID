clear
clc

%% LOAD SIMULATION DATA 
columnNames = {
    'rotationAngle'
    'axialPos'
    'rsectorID1'
    'moduleID1'
    'submoduleID1'
    'crystalID1'
    'globalPosX1'
    'globalPosY1'
    'globalPosZ1'
    'en1'
    'rsectorID2'
    'moduleID2'
    'submoduleID2'
    'crystalID2'
    'globalPosX2'
    'globalPosY2'
    'globalPosZ2'
    'en2'
    'rsectorID3'
    'moduleID3'
    'submoduleID3'
    'crystalID3'
    'globalPosX3'
    'globalPosY3'
    'globalPosZ3'
    'en3'
};

dataSim = load('/home/kolodziej/data_local/all_idealhits_120x15s_Ph_Ph.dat');
dataTableSim = array2table(dataSim, 'VariableNames', columnNames);

%% Get pairs of channels of simulation IDs in a range 0-511 (computed from rsector, module and crystal ID as below)
pairsSim = [dataTableSim.rsectorID1*128+dataTableSim.moduleID1*64+dataTableSim.crystalID1 dataTableSim.rsectorID2*128+dataTableSim.moduleID2*64+dataTableSim.crystalID2];

%% Sort the pairs, remove duplicates and get counts for each crystal pair (i.e., LOR)
pairsSortedSim = sort(pairsSim, 2);
[uniquePairsSim, ~, idSim] = unique(pairsSortedSim, 'rows');
countsSim = accumarray(idSim, 1);

%% Get counts for each pixel (i.e., crystal)
[firstIDSortedSim, ~, firstIDSim] = unique(pairsSortedSim(:,1));
[secondIDSortedSim, ~, secondIDSim] = unique(pairsSortedSim(:,2));
countsSimPerPixelTOP = accumarray(firstIDSim, 1);
countsSimPerPixelBOT = accumarray(secondIDSim, 1);
simCountsPerLOR = [uniquePairsSim countsSim];    
simCountsPerPixelTOP = [firstIDSortedSim countsSimPerPixelTOP]; 
simCountsPerPixelBOT = [secondIDSortedSim countsSimPerPixelBOT]; 

simCountsPerPixel = [simCountsPerPixelTOP; simCountsPerPixelBOT];

%% LOAD EXPERIMENTAL DATA 
dataExp = load('/home/kolodziej/data_local/coincidences_measurement_1.mat');
dataTableExp = dataExp.coincidences;

%% Get pairs of channels of raw experimental IDs in a range 0-511 (as on the PCB - unordered)
pairsExp = [dataTableExp(:,5), dataTableExp(:,10)];

%% Sort the pairs, remove duplicates and get counts for each crystal pair (i.e., LOR)
pairsSortedExp = sort(pairsExp, 2);        % Sort row-wise so (A,B) == (B,A)
[uniquePairsExp, ~, idExp] = unique(pairsSortedExp, 'rows');
countsExp = accumarray(idExp, 1);

expCountsPerLOR = [uniquePairsExp countsExp];

%% Get counts for each pixel (i.e., crystal)
[firstIDSortedExp, ~, firstIDExp] = unique(pairsSortedExp(:,1));
[secondIDSortedExp, ~, secondIDExp] = unique(pairsSortedExp(:,2));
countsExpPerPixelTOP = accumarray(firstIDExp, 1);
countsExpPerPixelBOT = accumarray(secondIDExp, 1);
expCountsPerPixelTOP = [firstIDSortedExp countsExpPerPixelTOP]; 
expCountsPerPixelBOT = [secondIDSortedExp countsExpPerPixelBOT];

expCountsPerPixel = [expCountsPerPixelTOP; expCountsPerPixelBOT];



%% Plot counts per LOR Simulation
x = simCountsPerLOR(:,1);
y = simCountsPerLOR(:,2);
z = simCountsPerLOR(:,3);

xMin = 0;          
yMin = 256;   
xIdx = x - xMin + 1;
yIdx = y - yMin + 1;
xSize = 256;                 
ySize = 256;               
zMap = NaN(ySize, xSize);

xMap = NaN(ySize, xSize);
yMap = NaN(ySize, xSize);

for i = 1:length(z)
    zMap(yIdx(i), xIdx(i)) = z(i);
    xMap(yIdx(i), xIdx(i)) = xIdx(i); % to have the same map sizes for all variables of 256*256
    yMap(yIdx(i), xIdx(i)) = yIdx(i); % idem
end
figure;
imagesc(0:255, 256:511, zMap);
axis xy;
colorbar;
axis equal tight;
ids=0:511;
n=512;
tick_step = 8;
tick_indices = 1:tick_step:n;
tick_labels = ids(tick_indices);

set(gca, 'XTick', tick_indices, 'XTickLabel', tick_labels, ...
         'YTick', tick_indices, 'YTickLabel', tick_labels);
xlabel('ID M1-M2');
ylabel('ID M3-M4');
title('2D Histogram of ID Pairs');

%% Prepare mapping for experiment (external information on the placement of the raw experimental IDs on the PCB board - hardcoded)
channel_map_cob_sipm1 = [0 2 11 4 3 13 9 7;21 1 5 27 29 17 24 39;43 45 23 49 44 35 55 59;61 56 53 62 58 57 60 63;
    26 6 22 8 14 16 20 10;28 12 33 18 30 37 15 19;25 51 31 50 34 41 46 47;52 54 38 48 42 40 36 32];
channel_map_cob_sipm2 = [0 2 5 4 3 17 9 7;21 1 23 11 13 35 24 39;43 45 27 49 44 29 55 59;56 61 62 53 57 58 63 60;
    6 26 8 22 16 14 10 20;28 12 31 18 30 41 15 19;25 51 38 33 37 40 46 47;52 54 50 48 42 34 36 32];

channel_map_bga_sipm1 = [1 0 8 10 3 15 12 5;13 14 7 11 27 6 4 9;62 58 2 61 59 38 50 57;52 56 51 55 60 54 63 53;
    16 21 18 20 23 17 22 19;25 24 34 26 28 33 29 30;31 37 32 41 40 36 42 35;44 46 39 45 47 43 48 49];
channel_map_bga_sipm2 = [1 0 7 10 3 6 12 5;13 14 2 8 15 38 4 9;62 58 11 61 59 27 50 57;56 52 55 51 54 60 53 63;
    21 16 20 18 17 23 19 22;25 24 32 26 28 36 29 30;31 37 39 34 33 43 42 35;44 46 41 45 47 40 48 49];

%% Add offsets to raw experimental IDs to have unique ID values; arrange them spatially
exp_mapping_module_0 = ([channel_map_bga_sipm1; channel_map_bga_sipm2+64]);
exp_mapping_module_1 = ([channel_map_bga_sipm1+128; channel_map_bga_sipm2+128+64]);
exp_mapping_module_2 = ([channel_map_cob_sipm1+2*128; channel_map_cob_sipm2+2*128+64]);
exp_mapping_module_3 = ([channel_map_bga_sipm1+3*128; channel_map_bga_sipm2+3*128+64]);
exp_mapping_all_modules = [exp_mapping_module_0, exp_mapping_module_1, exp_mapping_module_3, exp_mapping_module_2]; %modules 2 and 3 swapped here

%% Generate the simulation IDs with their spatial arrangement, add offsets for unique values
sim_SiPM1 = flipud(reshape(0:63, 8, 8));
sim_SiPM2 = sim_SiPM1+64;

sim_bothSiPMs_module0 = [sim_SiPM2; sim_SiPM1];
sim_bothSiPMs_module1 = [sim_SiPM2; sim_SiPM1]+128;
sim_bothSiPMs_module2 = [sim_SiPM2; sim_SiPM1]+2*128;
sim_bothSiPMs_module3 = [sim_SiPM2; sim_SiPM1]+3*128;
sim_mapping_all_modules = [sim_bothSiPMs_module0, sim_bothSiPMs_module1, sim_bothSiPMs_module2, sim_bothSiPMs_module3];

%% Create mapping table
exp_flat = exp_mapping_all_modules(:);
sim_flat = sim_mapping_all_modules(:);

mapping = sortrows([exp_flat, sim_flat], 1);

%% apply mapping to experimental data - counts per LOR
[~, idxM] = ismember(expCountsPerLOR(:,1), mapping(:,1));  % match data IDs to mapping IDs
new_columnM = mapping(idxM, 2);

[~, idxN] = ismember(expCountsPerLOR(:,2), mapping(:,1));  % match data IDs to mapping IDs
new_columnN = mapping(idxN, 2);

expCountsPerLOR = [expCountsPerLOR, new_columnM, new_columnN];

%% apply mapping to experimental data - counts per pixel
[~, idxA] = ismember(expCountsPerPixel(:,1), mapping(:,1));  % match data IDs to mapping IDs
new_columnA = mapping(idxA, 2);

expCountsPerPixel = [expCountsPerPixel, new_columnA]; %columns: raw_ID counts remapped_ID



%% plot counts per LOR Experiment
x1 = expCountsPerLOR(:,4);
y1 = expCountsPerLOR(:,5);
z1 = expCountsPerLOR(:,3); %counts

xMin1 = 0;          
yMin1 = 256;   

xIdx1 = x1 - xMin1 + 1;
yIdx1 = y1 - yMin1 + 1;
xSize1 = 256;                 
ySize1 = 256;               
zMapExp = NaN(ySize1, xSize1);
% xMapExp = NaN(ySize1, xSize1);
% yMapExp = NaN(ySize1, xSize1);

for i = 1:length(z1)
    zMapExp(yIdx1(i), xIdx1(i)) = z1(i);
    % xMapExp(yIdx1(i), xIdx1(i)) = xIdx1(i);
    % yMapExp(yIdx1(i), xIdx1(i)) = yIdx1(i);
end

figure;
imagesc(0:255, 256:511, zMapExp);
axis xy;
colorbar;
axis equal tight;
ids=0:511;
n=512;
tick_step = 8;
tick_indices = 1:tick_step:n;
tick_labels = ids(tick_indices);

set(gca, 'XTick', tick_indices, 'XTickLabel', tick_labels, ...
         'YTick', tick_indices, 'YTickLabel', tick_labels);
xlabel('ID M1-M2');
ylabel('ID M3-M4');
title('2D Histogram of ID Pairs');

%% DO NOT TRUST THE PLOTS PER PIXEL - to be checked.

%% new plots per pixel simulation
% Create SiPM mapping
sim_SiPM1 = flipud(reshape(0:63, 8, 8));
sim_SiPM2 = sim_SiPM1 + 64;
sim_bothSiPMs_module0 = [sim_SiPM2; sim_SiPM1];
sim_bothSiPMs_module1 = [sim_SiPM2; sim_SiPM1] + 128;
sim_bothSiPMs_module2 = [sim_SiPM2; sim_SiPM1] + 2*128;
sim_bothSiPMs_module3 = [sim_SiPM2; sim_SiPM1] + 3*128;

table_data = simCountsPerPixel;

% Create lookup map from ID to value
id_to_value = containers.Map(table_data(:, 1), table_data(:, 2));

% Verify each module is 16x8
modules = {sim_bothSiPMs_module1, sim_bothSiPMs_module0, ...
           sim_bothSiPMs_module2, sim_bothSiPMs_module3};
module_names = {'Module 1', 'Module 0', 'Module 2', 'Module 3'};


figure;
for idx = 1:4
    % Get module and rotate 90 degrees counterclockwise
    module = rot90(modules{idx});
    [rows, cols] = size(module);
    
    % Create color values from table lookup
    color_vals = zeros(rows, cols);
    for i = 1:rows
        for j = 1:cols
            pixel_id = module(i, j);
            color_vals(i, j) = id_to_value(pixel_id);
        end
    end
    
    subplot(2, 2, idx);
    imagesc(color_vals);
    colormap(parula);
    colorbar;

    hold on;
    
    % Add pixel ID labels
    for i = 1:rows
        for j = 1:cols
            pixel_id = module(i, j);
            text(j, i, num2str(pixel_id), 'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'middle', 'FontSize', 6, 'Color', 'white', ...
                'FontWeight', 'bold');
        end
    end
    
    title(module_names{idx}, 'FontSize', 12, 'FontWeight', 'bold');
    axis equal;
    xlim([0.5, cols+0.5]);
    ylim([0.5, rows+0.5]);
    set(gca, 'YDir', 'normal');
end

sgtitle('Simulation - counts per pixel', 'FontSize', 14, 'FontWeight', 'bold');

%% plot for exp.
table_data = expCountsPerPixel;
% id_to_value = containers.Map(table_data(:, 1), table_data(:, 2)); % for raw_exp IDs
id_to_value = containers.Map(table_data(:, 3), table_data(:, 2)); % for remapped_exp IDs

% Create figure with 4 subplots arranged as [m1 m0; m2 m3]
%figure('Position', [100, 100, 1000, 800]);
figure;

%positions = [1, 2; 3, 4]; % subplot layout: m1 m0 on top; m2 m3 on bottom
%module_list = [sim_bothSiPMs_module1, sim_bothSiPMs_module0, ...
%               sim_bothSiPMs_module2, sim_bothSiPMs_module3];
%module_names = {'Module 1', 'Module 0', 'Module 2', 'Module 3'};

for idx = 1:4
    % Get module and rotate 90 degrees counterclockwise
    module = rot90(modules{idx});
    [rows, cols] = size(module);
    
    % Create color values from table lookup
    color_vals = zeros(rows, cols);
    for i = 1:rows
        for j = 1:cols
            pixel_id = module(i, j);
            color_vals(i, j) = id_to_value(pixel_id);
        end
    end
    
    subplot(2, 2, idx);
    imagesc(color_vals);
    colormap(parula);
    colorbar;

    hold on;
    
    % Add pixel ID labels
    for i = 1:rows
        for j = 1:cols
            pixel_id = module(i, j);
            text(j, i, num2str(pixel_id), 'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'middle', 'FontSize', 6, 'Color', 'white', ...
                'FontWeight', 'bold');
        end
    end
    
    title(module_names{idx}, 'FontSize', 12, 'FontWeight', 'bold');
    axis equal;
    xlim([0.5, cols+0.5]);
    ylim([0.5, rows+0.5]);
    set(gca, 'YDir', 'normal');
end

sgtitle('Experiment - counts per pixel', 'FontSize', 14, 'FontWeight', 'bold');

% write mapping to file (optional)
writematrix(mapping, '/home/kolodziej/software/myMatlabScripts/mapping.txt', 'Delimiter', '\t');


%% Not used so far (requires other scripts to work) estimate the normalization components

m_ij = zMapExp(:);
s_ij = zMap(:);
n_crystals=512;
n_crystals_one_side=256;
i_crys = xMap(:);
j_crys = yMap(:);

% %fileID = fopen('/home/kolodziej/software/myMatlabScripts/Sij_sets.txt','w');
% %%S_ij sets - read from the above file
% %S_ij = 

[S_ij, LOR_indices] = createSetOfEquivalentLORs_v3();
[c, g, history] = estimate_normalization_components_v2(m_ij, s_ij, S_ij, LOR_indices, i_crys, j_crys, ...
     n_crystals, n_crystals_one_side, 50, 1e-3);
% % m_ij: Measured counts [LORs x 1]
% % s_ij: Simulated counts [LORs x 1]
% % i_crys, j_crys: crystal indices [LORs x 1]
% % n_crystals: total number of crystals
% % max_iter: maximum number of iterations
% % threshold: NRMS convergence threshold (e.g., 1e-3)