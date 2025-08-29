%clear
%data = load('/smbd/imageDataNucI/MERMAID/data/2025/2025-04-29/coincidences_measurement_1.mat');
%data = load('/smbd/imageDataNucI/MERMAID/data/2024/2024-08-28/F18_tube_1_mm_180s_at_14-10-50_qdc_singles.ldat');
clc
clear all
folder = '/smbd/imageDataNucI/MERMAID/data/2025/2025-04-29/';
%folder = '/smbd/imageDataNucI/MERMAID/data/2024/2024-02-13/';
filename = 'FDG_efficiency_large_cylinder_1800s_at_11-17-13_qdc.ldat';
%filename = 'Na22_1_mm_z_shift_300s_at_10-53-00_qdc_singles.ldat';
%filename = 'Na22_point_300s_at_11-29-31_qdc_singles.ldat';

[singles_raw] = readBinaryData(folder,filename);

