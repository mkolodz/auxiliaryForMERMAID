clear
clc
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

sim_mapping_module0 = [sim_SiPM2; sim_SiPM1];
sim_mapping_module1 = [sim_SiPM2; sim_SiPM1]+128;
sim_mapping_module2 = [sim_SiPM2; sim_SiPM1]+2*128;
sim_mapping_module3 = [sim_SiPM2; sim_SiPM1]+3*128;
sim_mapping_all_modules = [sim_mapping_module0, sim_mapping_module1, sim_mapping_module2, sim_mapping_module3];

%% Create mapping table
exp_flat = exp_mapping_all_modules(:);
sim_flat = sim_mapping_all_modules(:);

mapping = sortrows([exp_flat, sim_flat], 1);
