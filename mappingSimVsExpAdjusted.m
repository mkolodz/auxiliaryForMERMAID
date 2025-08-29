clear
clc

%LOAD SIMULATION DATA 
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
pairsSim = [dataTableSim.rsectorID1*128+dataTableSim.moduleID1*64+dataTableSim.crystalID1 dataTableSim.rsectorID2*128+dataTableSim.moduleID2*64+dataTableSim.crystalID2];

pairsSortedSim = sort(pairsSim, 2);
[uniquePairsSim, ~, idSim] = unique(pairsSortedSim, 'rows');
countsSim = accumarray(idSim, 1);

[firstIDSortedSim, ~, firstIDSim] = unique(pairsSortedSim(:,1));
[secondIDSortedSim, ~, secondIDSim] = unique(pairsSortedSim(:,2));
countsSimPerPixelTOP = accumarray(firstIDSim, 1);
countsSimPerPixelBOT = accumarray(secondIDSim, 1);

simCountsPerLOR = [uniquePairsSim countsSim];    
simCountsPerPixelTOP = [firstIDSortedSim countsSimPerPixelTOP]; 
simCountsPerPixelBOT = [secondIDSortedSim countsSimPerPixelBOT]; 
simCountsPerPixel = [simCountsPerPixelTOP; simCountsPerPixelBOT];

%LOAD EXPERIMENTAL DATA 
dataExp = load('/home/kolodziej/data_local/coincidences_measurement_1.mat');
dataTableExp = dataExp.coincidences;

pairsExp = [dataTableExp(:,5), dataTableExp(:,10)];      % Extract [chA, chB]
pairsSortedExp = sort(pairsExp, 2);        % Sort row-wise so (A,B) == (B,A)
[uniquePairsExp, ~, idExp] = unique(pairsSortedExp, 'rows');
countsExp = accumarray(idExp, 1);

expCountsPerLOR = [uniquePairsExp countsExp]; % [ch1 ch2 count]

[firstIDSortedExp, ~, firstIDExp] = unique(pairsSortedExp(:,1));
[secondIDSortedExp, ~, secondIDExp] = unique(pairsSortedExp(:,2));
countsExpPerPixelTOP = accumarray(firstIDExp, 1);
countsExpPerPixelBOT = accumarray(secondIDExp, 1);

expCountsPerPixelTOP = [firstIDSortedExp countsExpPerPixelTOP]; 
expCountsPerPixelBOT = [secondIDSortedExp countsExpPerPixelBOT]; 
expCountsPerPixel = [expCountsPerPixelTOP; expCountsPerPixelBOT];



% plot counts per LOR Simulation
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

for i = 1:length(z)
    zMap(yIdx(i), xIdx(i)) = z(i);
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

%prepare mapping for experiment
channel_map_cob_sipm1 = [0 2 11 4 3 13 9 7;21 1 5 27 29 17 24 39;43 45 23 49 44 35 55 59;61 56 53 62 58 57 60 63;
    26 6 22 8 14 16 20 10;28 12 33 18 30 37 15 19;25 51 31 50 34 41 46 47;52 54 38 48 42 40 36 32];
channel_map_cob_sipm2 = [0 2 5 4 3 17 9 7;21 1 23 11 13 35 24 39;43 45 27 49 44 29 55 59;56 61 62 53 57 58 63 60;
    6 26 8 22 16 14 10 20;28 12 31 18 30 41 15 19;25 51 38 33 37 40 46 47;52 54 50 48 42 34 36 32];

channel_map_bga_sipm1 = [1 0 8 10 3 15 12 5;13 14 7 11 27 6 4 9;62 58 2 61 59 38 50 57;52 56 51 55 60 54 63 53;
    16 21 18 20 23 17 22 19;25 24 34 26 28 33 29 30;31 37 32 41 40 36 42 35;44 46 39 45 47 43 48 49];
channel_map_bga_sipm2 = [1 0 7 10 3 6 12 5;13 14 2 8 15 38 4 9;62 58 11 61 59 27 50 57;56 52 55 51 54 60 53 63;
    21 16 20 18 17 23 19 22;25 24 32 26 28 36 29 30;31 37 39 34 33 43 42 35;44 46 41 45 47 40 48 49];

exp_mapping_module_0 = ([channel_map_bga_sipm1; channel_map_bga_sipm2+64]);
exp_mapping_module_1 = ([channel_map_bga_sipm1+128; channel_map_bga_sipm2+128+64]);
exp_mapping_module_2 = ([channel_map_cob_sipm1+2*128; channel_map_cob_sipm2+2*128+64]);
exp_mapping_module_3 = ([channel_map_bga_sipm1+3*128; channel_map_bga_sipm2+3*128+64]);
exp_mapping_all_modules = [exp_mapping_module_0, exp_mapping_module_1, exp_mapping_module_3, exp_mapping_module_2]; %modules 2 and 3 swapped here

sim_SiPM1 = flipud(reshape(0:63, 8, 8));
sim_SiPM2 = sim_SiPM1+64;

sim_bothSiPMs_module0 = [sim_SiPM2; sim_SiPM1];
sim_bothSiPMs_module1 = [sim_SiPM2; sim_SiPM1]+128;
sim_bothSiPMs_module2 = [sim_SiPM2; sim_SiPM1]+2*128;
sim_bothSiPMs_module3 = [sim_SiPM2; sim_SiPM1]+3*128;
sim_mapping_all_modules = [sim_bothSiPMs_module0, sim_bothSiPMs_module1, sim_bothSiPMs_module2, sim_bothSiPMs_module3];

exp_flat = exp_mapping_all_modules(:);
sim_flat = sim_mapping_all_modules(:);

mapping = sortrows([exp_flat, sim_flat], 1);

%apply mapping to experimental data 
[~, idxM] = ismember(expCountsPerLOR(:,1), mapping(:,1));  % match data IDs to mapping IDs
new_columnM = mapping(idxM, 2);

[~, idxN] = ismember(expCountsPerLOR(:,2), mapping(:,1));  % match data IDs to mapping IDs
new_columnN = mapping(idxN, 2);

expCountsPerLOR = [expCountsPerLOR, new_columnM, new_columnN];

% plot counts per LOR Experiment
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

for i = 1:length(z1)
    zMapExp(yIdx1(i), xIdx1(i)) = z1(i);
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

%DO NOT TRUST THE PLOTS PER PIXEL - to be checked.
% plot counts per Pixel Simulation
plotData = simCountsPerPixel(:,2);
lowerLim=2000;
upperLim=4000;

% lowerLim=0;
% upperLim=350;

n=64;
module1_1 = plotData(1:n);
module1_2 = plotData(n+1:2*n);
module2_1 = plotData(2*n+1:3*n);
module2_2 = plotData(3*n+1:4*n);
module3_1 = plotData(4*n+1:5*n);
module3_2 = plotData(5*n+1:6*n);
module4_1 = plotData(6*n+1:7*n);
module4_2 = plotData(7*n+1:8*n);
%A=1:64
module1_1_reshape=reshape(module1_1,[8,8]);
module1_2_reshape=reshape(module1_2,[8,8]);
module2_1_reshape=reshape(module2_1,[8,8]);
module2_2_reshape=reshape(module2_2,[8,8]);
module3_1_reshape=reshape(module3_1,[8,8]);
module3_2_reshape=reshape(module3_2,[8,8]);
module4_1_reshape=reshape(module4_1,[8,8]);
module4_2_reshape=reshape(module4_2,[8,8]);

module1=[module1_1_reshape;module1_2_reshape];
module2=[module2_1_reshape;module2_2_reshape];
module3=[module3_1_reshape;module3_2_reshape];
module4=[module4_1_reshape;module4_2_reshape];

figure;
subplot(2,2,1);
imagesc(module1);
axis image;
caxis([lowerLim upperLim]);
colorbar;
xlabel('Columns');
ylabel('Rows');
title('Module 1');

subplot(2,2,2);
imagesc(module2);
axis image;
caxis([lowerLim upperLim]);
colorbar;
xlabel('Columns');
ylabel('Rows');
title('Module 2');

subplot(2,2,3);
imagesc(module3);
axis image;
caxis([lowerLim upperLim]);
colorbar;
xlabel('Columns');
ylabel('Rows');
title('Module 3');

subplot(2,2,4);
imagesc(module4);
axis image;
caxis([lowerLim upperLim]);
colorbar;
xlabel('Columns');
ylabel('Rows');
title('Module 4');


% plot counts per Pixel Experiment
plotDataExp = expCountsPerPixel(:,2);
lowerLim=0;
upperLim=2500;

% lowerLim=0;
% upperLim=350;

n=64;
module1_1 = plotDataExp(1:n);
module1_2 = plotDataExp(n+1:2*n);
module2_1 = plotDataExp(2*n+1:3*n);
module2_2 = plotDataExp(3*n+1:4*n);
module3_1 = plotDataExp(4*n+1:5*n);
module3_2 = plotDataExp(5*n+1:6*n);
module4_1 = plotDataExp(6*n+1:7*n);
module4_2 = plotDataExp(7*n+1:8*n);
%A=1:64
module1_1_reshape=reshape(module1_1,[8,8]);
module1_2_reshape=reshape(module1_2,[8,8]);
module2_1_reshape=reshape(module2_1,[8,8]);
module2_2_reshape=reshape(module2_2,[8,8]);
module3_1_reshape=reshape(module3_1,[8,8]);
module3_2_reshape=reshape(module3_2,[8,8]);
module4_1_reshape=reshape(module4_1,[8,8]);
module4_2_reshape=reshape(module4_2,[8,8]);

module1=[module1_1_reshape;module1_2_reshape];
module2=[module2_1_reshape;module2_2_reshape];
module3=[module3_1_reshape;module3_2_reshape];
module4=[module4_1_reshape;module4_2_reshape];

figure;
subplot(2,2,1);
imagesc(module1);
axis image;
caxis([lowerLim upperLim]);
colorbar;
xlabel('Columns');
ylabel('Rows');
title('Module 1');

subplot(2,2,2);
imagesc(module2);
axis image;
caxis([lowerLim upperLim]);
colorbar;
xlabel('Columns');
ylabel('Rows');
title('Module 2');

subplot(2,2,3);
imagesc(module3);
axis image;
caxis([lowerLim upperLim]);
colorbar;
xlabel('Columns');
ylabel('Rows');
title('Module 3');

subplot(2,2,4);
imagesc(module4);
axis image;
caxis([lowerLim upperLim]);
colorbar;
xlabel('Columns');
ylabel('Rows');
title('Module 4');