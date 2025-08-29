% clear
% % Create a random 16x8 matrix (replace this with your actual data)
% data = rand(16, 8);  % Example matrix (16 rows, 8 columns)
% 
% % Create the heatmap
% figure;
% heatmap(data);
% % Add labels for the axes (optional)
% xlabel('Columns');
% ylabel('Rows');
% title('Heatmap of 16x8 Matrix');
% drawnow;
% 
% % plot(1:10, rand(1, 10));

clear
% load("/home/kolodziej/Downloads/Measurement_Decay_Heat_Maps_Module_1_long.mat");
load("/home/kolodziej/Downloads/Measurement_Decay_Heat_Maps_Module_1.mat");
figure;
subplot(2,2,1);
imagesc(heat_map_module_1_kev);
axis equal;
caxis([0 110000]);
colorbar;
% axis([0.5 size(data,2)+0.5 0.5 size(data,1)+0.5]);
xlabel('Columns');
ylabel('Rows');
title('Module 1');

% load("/home/kolodziej/Downloads/Measurement_Decay_Heat_Maps_Module_2_long.mat");
load("/home/kolodziej/Downloads/Measurement_Decay_Heat_Maps_Module_2.mat");
subplot(2,2,2);
imagesc(heat_map_module_2_kev);
axis equal;
caxis([0 110000]);
colorbar;
% axis([0.5 size(data,2)+0.5 0.5 size(data,1)+0.5]);
xlabel('Columns');
ylabel('Rows');
title('Module 2');

% load("/home/kolodziej/Downloads/Measurement_Decay_Heat_Maps_Module_3_long.mat");
load("/home/kolodziej/Downloads/Measurement_Decay_Heat_Maps_Module_3.mat");

subplot(2,2,3);
imagesc(heat_map_module_3_kev);
axis equal;
caxis([0 110000]);
colorbar;
% axis([0.5 size(data,2)+0.5 0.5 size(data,1)+0.5]);
xlabel('Columns');
ylabel('Rows');
title('Module 3');

% load("/home/kolodziej/Downloads/Measurement_Decay_Heat_Maps_Module_4_long.mat");
load("/home/kolodziej/Downloads/Measurement_Decay_Heat_Maps_Module_4.mat");
subplot(2,2,4);
imagesc(heat_map_module_4_kev);
axis equal;
caxis([0 110000]);
colorbar;
% axis([0.5 size(data,2)+0.5 0.5 size(data,1)+0.5]);
xlabel('Columns');
ylabel('Rows');
title('Module 4');
