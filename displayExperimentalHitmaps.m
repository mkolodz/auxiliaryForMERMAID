%% auxiliary for experimental mapping visualisation, no input file required

module_offset=3*128;

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
disp(exp_mapping_module_1+module_offset);
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