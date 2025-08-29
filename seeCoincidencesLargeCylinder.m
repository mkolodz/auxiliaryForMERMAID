clear
clc
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
% pairsSim = [dataTableSim.moduleID1*64+dataTableSim.crystalID1 128+dataTableSim.moduleID2*64+dataTableSim.crystalID2];
% disp([dataTableSim.moduleID1, dataTableSim.crystalID1, dataTableSim.moduleID1*64+dataTableSim.crystalID1])
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

%LOAD MAPPING
%mapping_table = load("/home/kolodziej/Documents/IntroMaterialsMermaid/MAPPING_singlesgeometrySimuVsExp_rotStep0.txt");
mapping_table = readtable("/home/kolodziej/Documents/IntroMaterialsMermaid/MAPPING_singlesgeometrySimuVsExp_rotStep0.txt",'Delimiter', '\t','ReadVariableNames', true);
mapping_table.channelIDsim = mapping_table.rsectorID*128+mapping_table.moduleID*64+mapping_table.crystalID;
map_only=[mapping_table.channelIDexp mapping_table.channelIDsim];

mappingPerModule = readtable("/home/kolodziej/Documents/NormalizationPET/mapping/MappingChannelsExpVsSim.txt");
mappingPerModuleSwappedSiPMs = readtable("/home/kolodziej/Documents/NormalizationPET/mapping/MappingChannelsExpVsSim_swappedSiPMsInSimulation.txt");
mappingPerModule = table2array(mappingPerModule);
mappingPerModuleSwappedSiPMs = table2array(mappingPerModuleSwappedSiPMs);

%mappingPerModule = mappingPerModuleSwappedSiPMs;


%LOAD EXPERIMENTAL DATA 
dataExp = load('/home/kolodziej/data_local/coincidences_measurement_1.mat');
dataTableExp = dataExp.coincidences;

pairsExp = [dataTableExp(:,5), dataTableExp(:,10)];      % Extract [chA, chB]
pairsSortedExp = sort(pairsExp, 2);        % Sort row-wise so (A,B) == (B,A)
[uniquePairsExp, ~, idExp] = unique(pairsSortedExp, 'rows');
countsExp = accumarray(idExp, 1);

expCountsPerLOR = [uniquePairsExp countsExp];        % [ch1 ch2 count]

disp("dataTableSim.crystalID1 RANGE " + range(dataTableSim.crystalID1))
disp("dataTableSim.rsectorID1 RANGE " + range(dataTableSim.rsectorID1))
disp("dataTableSim.moduleID1 RANGE " + range(dataTableSim.moduleID1))


geom1 = load('/home/kolodziej/data_local/geometry.mat');


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

% Fill the zMap
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

channel_map_cob_sipm1 = [0 2 11 4 3 13 9 7;21 1 5 27 29 17 24 39;43 45 23 49 44 35 55 59;61 56 53 62 58 57 60 63;
    26 6 22 8 14 16 20 10;28 12 33 18 30 37 15 19;25 51 31 50 34 41 46 47;52 54 38 48 42 40 36 32];
channel_map_cob_sipm2 = [0 2 5 4 3 17 9 7;21 1 23 11 13 35 24 39;43 45 27 49 44 29 55 59;56 61 62 53 57 58 63 60;
    6 26 8 22 16 14 10 20;28 12 31 18 30 41 15 19;25 51 38 33 37 40 46 47;52 54 50 48 42 34 36 32];

channel_map_bga_sipm1 = [1 0 8 10 3 15 12 5;13 14 7 11 27 6 4 9;62 58 2 61 59 38 50 57;52 56 51 55 60 54 63 53;
    16 21 18 20 23 17 22 19;25 24 34 26 28 33 29 30;31 37 32 41 40 36 42 35;44 46 39 45 47 43 48 49];
channel_map_bga_sipm2 = [1 0 7 10 3 6 12 5;13 14 2 8 15 38 4 9;62 58 11 61 59 27 50 57;56 52 55 51 54 60 53 63;
    21 16 20 18 17 23 19 22;25 24 32 26 28 36 29 30;31 37 39 34 33 43 42 35;44 46 41 45 47 40 48 49];
%events = zeros(6,length(coincidence));
%events = zeros(6,1);
disp([channel_map_bga_sipm1; channel_map_bga_sipm2]);
exp_mapping_module_1 = ([channel_map_bga_sipm1; channel_map_bga_sipm2+64]);
disp(exp_mapping_module_1);
% lookup matrix for geometry position (saved from rotation_gategeometrie_newModules) for SiPM1 and SiPM2
% position of S1 and S2 coud be fliped 180° --> check in prototype setup
geom_lookup_s1 = [1:16:113;2:16:114;3:16:115;4:16:116;5:16:117;6:16:118;7:16:119;8:16:120];
geom_lookup_s2 = [9:16:121;10:16:122;11:16:123;12:16:124;13:16:125;14:16:126;15:16:127;16:16:128];

module_raw_bga=[channel_map_bga_sipm1,channel_map_bga_sipm2+64];
module_raw_cob=[channel_map_cob_sipm1,channel_map_cob_sipm2+64];
module_mapped=[geom_lookup_s1, geom_lookup_s2];
disp("\n")
sim_SiPM1 = flipud(reshape(0:63, 8, 8));
sim_SiPM2 = sim_SiPM1+64;
sim_bothSiPMs = [sim_SiPM2; sim_SiPM1];

%disp([geom_lookup_s1;geom_lookup_s2]);
disp(sim_bothSiPMs);
all_modules_raw = [module_raw_bga, module_raw_bga+128, module_raw_cob+2*128, module_raw_bga+3*128];
all_modules_mapped = [module_mapped, module_mapped+128, module_mapped+2*128, module_mapped+3*128];

raw_flat = all_modules_raw(:);
mapped_flat = all_modules_mapped(:);

% Combine into a 2-column table
%channel_mapping_table = [raw_flat, mapped_flat];
channel_mapping_sorted = sortrows([raw_flat, mapped_flat], 1);
% 
% figure;
% histogram(expCountsPerLOR(:,3));
% figure;
% histogram(simCountsPerLOR(:,3));


% plot counts per LOR Experiment ch1 - channels 0-255 (TOP side); ch2 - channels
% 256-511 (BOTTOM side)
% columns meaning: ch1raw (0-255), ch2raw (256-511), counts, ch1remapped, ch2remapped, 
% ch1remappedSwappedModules, ch2remappedSwappedModules, 
% ch1remappedSwappedModules_remappedWithinModule, ch2remappedSwappedModules_remappedWithinModule

%remappedExpCountsPerLOR = expCountsPerLOR;
[~, idx1] = ismember(expCountsPerLOR(:,1), map_only(:,1));  % match data IDs to mapping IDs
new_column1 = map_only(idx1, 2);
%remappedExpCountsPerLOR = [expCountsPerLOR, new_column];

[~, idx2] = ismember(expCountsPerLOR(:,2), map_only(:,1));  % match data IDs to mapping IDs
new_column2 = map_only(idx2, 2);
expCountsPerLOR = [expCountsPerLOR, new_column1, new_column2];

[~, idxA] = ismember(expCountsPerLOR(:,1), channel_mapping_sorted(:,1));  % match data IDs to mapping IDs
new_columnA = channel_mapping_sorted(idxA, 2);
%remappedExpCountsPerLOR = [expCountsPerLOR, new_column];

[~, idxB] = ismember(expCountsPerLOR(:,2), channel_mapping_sorted(:,1));  % match data IDs to mapping IDs
new_columnB = channel_mapping_sorted(idxB, 2);
expCountsPerLOR = [expCountsPerLOR, new_columnA, new_columnB];
% 
% [~, idxC] = ismember(expCountsPerLOR(:,6), mappingPerModule(:,1));  % match data IDs to mapping IDs
% new_columnC = mappingPerModule(idxC, 2);
% %remappedExpCountsPerLOR = [expCountsPerLOR, new_column];
% 
% [~, idxD] = ismember(expCountsPerLOR(:,7), mappingPerModule(:,1));  % match data IDs to mapping IDs
% new_columnD = mappingPerModule(idxD, 2);
% expCountsPerLOR = [expCountsPerLOR, new_columnC, new_columnD];


%expCountsPerLOR = [expCountsPerLOR, channel_mapping_sorted(1:256,2), channel_mapping_sorted(257:512,2)];

% for i = 1:height(expCountsPerLOR)
%     if (expCountsPerLOR(i,4) >= 0 && expCountsPerLOR(i,4) <= 127)
 %this below but first subtract then add back the threshold for
        %different modules.
        % idWithinModule = expCountsPerLOR(i,4);
        % [~, idx_m0] = ismember(expCountsPerLOR(:,1), mappingPerModule(:,1));  % match data IDs to mapping IDs
        % new_columnA = mappingPerModule(idx_m0, 2);
%         expCountsPerLOR(i,6) = expCountsPerLOR(i,4)+128;
%     end
%     if (expCountsPerLOR(i,4) >= 128 && expCountsPerLOR(i,4) <= 255)
%         expCountsPerLOR(i,6) = expCountsPerLOR(i,4)-128;
%     end
% end
% 
% for i = 1:height(expCountsPerLOR)
%     if (expCountsPerLOR(i,5) >= 256 && expCountsPerLOR(i,5) <= 383)
%         expCountsPerLOR(i,7) = expCountsPerLOR(i,5)+128;
%     end
%     if (expCountsPerLOR(i,5) >= 384 && expCountsPerLOR(i,5) <= 511)
%         expCountsPerLOR(i,7) = expCountsPerLOR(i,5)-128;
%     end
% end


% %remapping sim indices - rotation 180deg + swapping modules 2 and 3
% for i = 1:height(map_only)
%     if (map_only(i,2) >= 0 && map_only(i,2) <= 127)
%        map_only(i,3) = map_only(i,2)+384; %0->3
%     end
%     if (map_only(i,2) >= 128 && map_only(i,2) <= 255)
%        map_only(i,3) = map_only(i,2)+128; %1->2
%     end
%     if (map_only(i,2) >= 256 && map_only(i,2) <= 383)
%        map_only(i,3) = map_only(i,2)-256; %2->0
%     end
%     if (map_only(i,2) >= 384 && map_only(i,2) <= 511)
%        map_only(i,3) = map_only(i,2)-256; %3->1
%     end
% end

% %remapping sim indices - remapping within modules to exp convention
% for i = 1:height(map_only)
%     if (map_only(i,3) >= 0 && map_only(i,3) <= 127)
%        [~, idx_m0] = ismember(map_only(i,3)-0, mappingPerModule(:,2));  % COL2-SIM match data IDs to mapping IDs
%         map_only(i,4) = mappingPerModule(idx_m0, 1)+0; % COL1-EXP
%     end
%     if (map_only(i,3) >= 128 && map_only(i,3) <= 255)
%        [~, idx_m1] = ismember(map_only(i,3)-128, mappingPerModule(:,2));  % COL2-SIM match data IDs to mapping IDs
%         map_only(i,4) = mappingPerModule(idx_m1, 1)+128; % COL1-EXP
%     end
%     if (map_only(i,3) >= 256 && map_only(i,3) <= 383)
%        [~, idx_m2] = ismember(map_only(i,3)-128*2, mappingPerModule(:,2));  % COL2-SIM match data IDs to mapping IDs
%         map_only(i,4) = mappingPerModule(idx_m2, 1)+128*2; % COL1-EXP
%     end
%     if (map_only(i,3) >= 384 && map_only(i,3) <= 511)
%        [~, idx_m3] = ismember(map_only(i,3)-128*3, mappingPerModule(:,2));  % COL2-SIM match data IDs to mapping IDs
%         map_only(i,4) = mappingPerModule(idx_m3, 1)+128*3; % COL1-EXP
%     end
% end
% 
% [~, idxA] = ismember(expCountsPerLOR(:,1), map_only(:,1));  % match data IDs to mapping IDs
% new_columnA = map_only(idxA, 4);
% %remappedExpCountsPerLOR = [expCountsPerLOR, new_column];
% 
% [~, idxB] = ismember(expCountsPerLOR(:,2), map_only(:,1));  % match data IDs to mapping IDs
% new_columnB = map_only(idxB, 4);
% expCountsPerLOR = [expCountsPerLOR, new_columnA, new_columnB];
% 


% %below, sim?? modules are remapped to match the exp?? modules this way:1->2,
% %0->3, 2->0, 3->1 (rotation 180deg + swapping modules 2 and 3)
% 
% for i = 1:height(expCountsPerLOR)
%     if (expCountsPerLOR(i,4) >= 0 && expCountsPerLOR(i,4) <= 127)
%        expCountsPerLOR(i,6) = expCountsPerLOR(i,4)+384; %0->3
%     end
%     if (expCountsPerLOR(i,4) >= 128 && expCountsPerLOR(i,4) <= 255)
%        expCountsPerLOR(i,6) = expCountsPerLOR(i,4)+128; %1->2
%     end
%     if (expCountsPerLOR(i,5) >= 256 && expCountsPerLOR(i,5) <= 383)
%        expCountsPerLOR(i,7) = expCountsPerLOR(i,5)-256; %2->0
%     end
%     if (expCountsPerLOR(i,5) >= 384 && expCountsPerLOR(i,5) <= 511)
%        expCountsPerLOR(i,7) = expCountsPerLOR(i,5)-256; %3->1
%     end
% end

% below, exp?? modules are remapped within modules to match the sim mapping:
for i = 1:height(expCountsPerLOR)
    if (expCountsPerLOR(i,6) >= 0 && expCountsPerLOR(i,6) <= 127)
        [~, idx_m0] = ismember(expCountsPerLOR(i,6)-0, mappingPerModule(:,1));  % COL2-SIM match data IDs to mapping IDs
        expCountsPerLOR(i,9) = mappingPerModule(idx_m0, 2)+0; % COL1-EXP
    end
    if (expCountsPerLOR(i,6) >= 128 && expCountsPerLOR(i,6) <= 255)
        [~, idx_m1] = ismember(expCountsPerLOR(i,6)-128, mappingPerModule(:,1));  % COL2-SIM match data IDs to mapping IDs
        expCountsPerLOR(i,9) = mappingPerModule(idx_m1, 2)+128; % COL1-EXP
    end
    if (expCountsPerLOR(i,7) >= 256 && expCountsPerLOR(i,7) <= 383)
        [~, idx_m2] = ismember(expCountsPerLOR(i,7)-2*128, mappingPerModule(:,1));  % COL2-SIM match data IDs to mapping IDs
        expCountsPerLOR(i,8) = mappingPerModule(idx_m2, 2)+2*128; % COL1-EXP
    end
    if (expCountsPerLOR(i,7) >= 384 && expCountsPerLOR(i,7) <= 511)
        [~, idx_m3] = ismember(expCountsPerLOR(i,7)-3*128, mappingPerModule(:,1));  % COL2-SIM match data IDs to mapping IDs
        %disp(mappingPerModule(idx_m3, 1))
        expCountsPerLOR(i,8) = mappingPerModule(idx_m3, 2)+3*128; % COL1-EXP
    end
end

% sortedExpCountsPerLOR = expCountsPerLOR;
% sortedExpCountsPerLOR = sortrows(sortedExpCountsPerLOR, 4);
% sortedExpCountsPerLOR(:,4) = flipud(sortedExpCountsPerLOR(:,4));
% sortedExpCountsPerLOR = sortrows(sortedExpCountsPerLOR, 6);
% sortedExpCountsPerLOR(:,6) = flipud(sortedExpCountsPerLOR(:,6));
% swap modules 3 and 4 to match the experiment mapping


%expCountsPerLOR = [expCountsPerLOR, new_column3];

% x1 = expCountsPerLOR(:,1);
% y1 = expCountsPerLOR(:,2);
% x1 = expCountsPerLOR(:,6); %ch1remappedSwappedModules
% y1 = expCountsPerLOR(:,5); %ch2remapped
x1 = expCountsPerLOR(:,8); %6
y1 = expCountsPerLOR(:,9); %7
z1 = expCountsPerLOR(:,3); %counts

xMin1 = 0;          
yMin1 = 0;   

% xMin1 = 256;          
% yMin1 = 0; 

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

m_ij = zMapExp(:);
s_ij = zMap(:);
n_crystals=512;
% i_crys
% j_crys
% [c, g, history] = estimate_normalization_components(m_ij, s_ij, i_crys, j_crys, ...
%     n_crystals, 50, 1e-3);