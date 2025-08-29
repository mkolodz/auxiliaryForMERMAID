clear
clc
data = load('/smbd/imageDataNucI/MERMAID/data/2025/2025-04-29/coincidences_measurement_1.mat');
tableData = data.coincidences;
pairs = [tableData(:,5), tableData(:,10)];      % Extract [chA, chB]
pairs_sorted = sort(pairs, 2);        % Sort so (A,B) == (B,A)

[uniquePairs, ~, idx] = unique(pairs_sorted, 'rows');
counts = accumarray(idx, 1);

result = [uniquePairs counts];        % [ch1 ch2 count]

% histogram(result(:,3));   
% xlim([-50, 50]);% Plot histogram of counts
% xlabel('# of coincidences');
% ylabel('# of channel pairs');
% title('Large cylinder, measurement 1');

% figure;
% histogram(tableData(:,5), 512); 
%histogram(tableData(:,10), 512);