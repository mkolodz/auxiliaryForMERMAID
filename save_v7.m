clc
clear all

parameter_loadfolder = '/home/kolodziej/data_local/';
parameter_loadname = 'parameters_calibrationmeasurmenet_2023_07_06_offset_0.mat';
myfile=[parameter_loadfolder parameter_loadname];
olddata = load(myfile);
whos -file '/home/kolodziej/data_local/parameters_calibrationmeasurmenet_2023_07_06_offset_0.mat'
% olddata=load([parameter_loadfolder parameter_loadname]);
save('newcalibdata.mat', kev_calibration_parameterlist, '-v7');
% save('newcalibdata.mat', 'myvar', '-v7');


% umrechnung durchführen

% bei großen Daten darauf achten, dass nur das notwendigste im Workspace
% ist
%[singles_keV] = convert_DAC_to_keV_all(singles_raw,kev_calibration_parameterlist);