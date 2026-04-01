clear all;
clc;

% This script is for plotting the spectra from new (from 2026) PET modules and comparing them for different PETsys parameters. 
% It requires one function from the main MERMAID code, in the slightly edited version:   
% 1) /home/kolodziej/software/mermaid/multi_detector_scripts/energyTresh calibration gui mD/energySpectraSingles_mD_noFit.m
% it has to be run with the following path: /home/kolodziej/software/mermaid/multi_detector_scripts/energyTresh calibration gui mD/

%tunable params:
hist_high = 400;
maxbinLim = 50;
binWidth = 0.1;

folder = '/smbd/imageDataNucI/Magdalena_Kolodziej/newPETModulesTests/data/';
filename_datasetA = 'parameterScan_OV_3.0V_single.ldat';
filename_datasetB = 'parameterScan_stabilityCheckAfter_defaultParameters_single.ldat';
filename_datasetC = 'afterRecalibration_defaultParameters_single.ldat';

% filename_datasetA = 'parameterScan_OV_3.0V_single.ldat';
% filename_datasetB = 'parameterScan_lsb_t1_40_single.ldat';
% filename_datasetC = 'parameterScan_lsb_t1_60_single.ldat';

% filename_datasetA = 'parameterScan_OV_2.5V_single.ldat';
% filename_datasetB = 'parameterScan_OV_3.0V_single.ldat';
% filename_datasetC = 'parameterScan_OV_3.5V_single.ldat';

% filename_datasetA = 'parameterScan_thrE_22_single.ldat';
% filename_datasetB = 'parameterScan_thrE_23_single.ldat';
% filename_datasetC = 'parameterScan_thrE_24_single.ldat';
% filename_datasetD = 'parameterScan_thrE_25_single.ldat';
% filename_datasetE = 'parameterScan_thrE_26_single.ldat';

[singles_datasetA] = readBinaryData(folder,filename_datasetA);
[singles_datasetB] = readBinaryData(folder,filename_datasetB);
[singles_datasetC] = readBinaryData(folder,filename_datasetC);
% [singles_datasetD] = readBinaryData(folder,filename_datasetD);
% [singles_datasetE] = readBinaryData(folder,filename_datasetE);
%energy spectra evaluation for all channels:
[~,~,~,detector_module_4_A] = energySpectraSingles_mD_noFit(singles_datasetA);
[~,~,~,detector_module_4_B] = energySpectraSingles_mD_noFit(singles_datasetB);
[~,~,~,detector_module_4_C] = energySpectraSingles_mD_noFit(singles_datasetC);
% [~,~,~,detector_module_4_D] = energySpectraSingles_mD_noFit(singles_datasetD);
% [~,~,~,detector_module_4_E] = energySpectraSingles_mD_noFit(singles_datasetE);
detector_module_4_D = cell(128,2);
detector_module_4_E = cell(128,2);
module_nr = 4;

photopeak_pos_list = zeros(512,3);

channel_map_bga_sipm1 = [1 0 8 10 3 15 12 5;13 14 7 11 27 6 4 9;62 58 2 61 59 38 50 57;52 56 51 55 60 54 63 53;
                         16 21 18 20 23 17 22 19;25 24 34 26 28 33 29 30;31 37 32 41 40 36 42 35;44 46 39 45 47 43 48 49];
channel_map_bga_sipm2 = [1 0 7 10 3 6 12 5;13 14 2 8 15 38 4 9;62 58 11 61 59 27 50 57;56 52 55 51 54 60 53 63;
                         21 16 20 18 17 23 19 22;25 24 32 26 28 36 29 30;31 37 39 34 33 43 42 35;44 46 41 45 47 40 48 49];

%for first SiPM
figuretitle1 = ['Energy spectra module ', num2str(module_nr),', SiPM 1']; 
fig_1 = figure('NumberTitle','off','Name',figuretitle1);
title('Energy spectrum per channel');

%for second SiPM
figuretitle2 = ['Energy spectra module ', num2str(module_nr),', SiPM 2'];
fig_2 = figure('NumberTitle','off','Name',figuretitle2);
title('Energy spectrum per channel');


for j = 1:length(detector_module_4_A)
    % Skip empty entries in A
    if isempty(detector_module_4_A{j,1})
        continue
    end

    channelID = detector_module_4_A{j,1};

    % Find matching index in B, C, D, E by channel ID
    idx_B = find(cellfun(@(x) ~isempty(x) && x == channelID, detector_module_4_B(:,1)), 1);
    idx_C = find(cellfun(@(x) ~isempty(x) && x == channelID, detector_module_4_C(:,1)), 1);
    idx_D = find(cellfun(@(x) ~isempty(x) && x == channelID, detector_module_4_D(:,1)), 1);
    idx_E = find(cellfun(@(x) ~isempty(x) && x == channelID, detector_module_4_E(:,1)), 1);

    % Module 4 SiPM selection
    if channelID < 448
        channelID_1 = channelID - 64*6;
        set(0,'CurrentFigure', fig_1)
        for m = 1:8
            for n = 1:8
                if channel_map_bga_sipm1(n,m) == channelID_1
                    pos_on_det = (n*8) - (8-m);
                    s(pos_on_det) = subplot(8,8,pos_on_det);

                    histogram(detector_module_4_A{j,2}, 'BinWidth',binWidth, ...
                        'BinLimits',[0 maxbinLim],'DisplayStyle','stairs','EdgeColor','b');
                    hold on
                    if ~isempty(idx_B), histogram(detector_module_4_B{idx_B,2}, 'BinWidth',binWidth,'BinLimits',[0 maxbinLim],'DisplayStyle','stairs','EdgeColor','r'); end
                    if ~isempty(idx_C), histogram(detector_module_4_C{idx_C,2}, 'BinWidth',binWidth,'BinLimits',[0 maxbinLim],'DisplayStyle','stairs','EdgeColor','g'); end
                    if ~isempty(idx_D), histogram(detector_module_4_D{idx_D,2}, 'BinWidth',binWidth,'BinLimits',[0 maxbinLim],'DisplayStyle','stairs','EdgeColor','c'); end
                    if ~isempty(idx_E), histogram(detector_module_4_E{idx_E,2}, 'BinWidth',binWidth,'BinLimits',[0 maxbinLim],'DisplayStyle','stairs','EdgeColor','m'); end

                    title(s(pos_on_det), ['Channel ', num2str(channelID)]);
                    xlim([0 maxbinLim]);
                    ylim([0 hist_high]);
                end
            end
        end
    else
        channelID_2 = channelID - 64*7;
        set(0,'CurrentFigure', fig_2)
        for m = 1:8
            for n = 1:8
                if channel_map_bga_sipm2(n,m) == channelID_2
                    pos_on_det = (n*8) - (8-m);
                    s(pos_on_det) = subplot(8,8,pos_on_det);

                    histogram(detector_module_4_A{j,2}, 'BinWidth',binWidth, ...
                        'BinLimits',[0 maxbinLim],'DisplayStyle','stairs','EdgeColor','b');
                    hold on
                    if ~isempty(idx_B), histogram(detector_module_4_B{idx_B,2}, 'BinWidth',binWidth,'BinLimits',[0 maxbinLim],'DisplayStyle','stairs','EdgeColor','r'); end
                    if ~isempty(idx_C), histogram(detector_module_4_C{idx_C,2}, 'BinWidth',binWidth,'BinLimits',[0 maxbinLim],'DisplayStyle','stairs','EdgeColor','g'); end
                    if ~isempty(idx_D), histogram(detector_module_4_D{idx_D,2}, 'BinWidth',binWidth,'BinLimits',[0 maxbinLim],'DisplayStyle','stairs','EdgeColor','c'); end
                    if ~isempty(idx_E), histogram(detector_module_4_E{idx_E,2}, 'BinWidth',binWidth,'BinLimits',[0 maxbinLim],'DisplayStyle','stairs','EdgeColor','m'); end

                    title(s(pos_on_det), ['Channel ', num2str(channelID)]);
                    xlim([0 maxbinLim]);
                    ylim([0 hist_high]);
                end
            end
        end
    end
end



% %same plotting, but selected channels only:
% 
% %selected_channels = [386, 411, 504, 508];
% %selected_channels = [466, 454, 401, 402];
% selected_channels = [500, 472, 438, 410];
selected_channels = [386, 411, 504, 508];

fig_sel = figure;
sgtitle('Selected Channels');

for p = 1:length(selected_channels)
    channelID = selected_channels(p);

    % Find matching index in A, B, C, D, E by channel ID
    idx_A = find(cellfun(@(x) ~isempty(x) && x == channelID, detector_module_4_A(:,1)), 1);
    idx_B = find(cellfun(@(x) ~isempty(x) && x == channelID, detector_module_4_B(:,1)), 1);
    idx_C = find(cellfun(@(x) ~isempty(x) && x == channelID, detector_module_4_C(:,1)), 1);
    idx_D = find(cellfun(@(x) ~isempty(x) && x == channelID, detector_module_4_D(:,1)), 1);
    idx_E = find(cellfun(@(x) ~isempty(x) && x == channelID, detector_module_4_E(:,1)), 1);

    subplot(2, 2, p);
    hold on;

    if ~isempty(idx_A), histogram(detector_module_4_A{idx_A,2}, 'BinWidth',binWidth,'BinLimits',[0 maxbinLim],'DisplayStyle','stairs','EdgeColor','b'); end
    if ~isempty(idx_B), histogram(detector_module_4_B{idx_B,2}, 'BinWidth',binWidth,'BinLimits',[0 maxbinLim],'DisplayStyle','stairs','EdgeColor','r'); end
    if ~isempty(idx_C), histogram(detector_module_4_C{idx_C,2}, 'BinWidth',binWidth,'BinLimits',[0 maxbinLim],'DisplayStyle','stairs','EdgeColor','g'); end
    if ~isempty(idx_D), histogram(detector_module_4_D{idx_D,2}, 'BinWidth',binWidth,'BinLimits',[0 maxbinLim],'DisplayStyle','stairs','EdgeColor','c'); end
    if ~isempty(idx_E), histogram(detector_module_4_E{idx_E,2}, 'BinWidth',binWidth,'BinLimits',[0 maxbinLim],'DisplayStyle','stairs','EdgeColor','m'); end

    title(['Channel ', num2str(channelID)]);
    xlim([0 maxbinLim]);
    ylim([0 hist_high]);
    xlabel('Energy in DAC Units');
    ylabel('Counts');
    %legend('vth_e=22','vth_e=23','vth_e=24','vth_e=25','vth_e=26','Location','best');
    %legend('OV=2.5V', 'OV=3.0V', 'OV=3.5V','Location','best');
    %legend('vth1=vth2=20', 'vth1=vth2=22', 'vth1=vth2=24','Location','best');
    %legend('disc_lsb_vth1=51','disc_lsb_vth1=40', 'disc_lsb_vth1=60','Location','best');
    legend('\default','default repeated', 'recalibrated - wrong config', 'Location','best');
    hold off;
end