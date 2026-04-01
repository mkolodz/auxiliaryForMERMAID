clear
clc

IDtype = "raw"; % options: "remapped", "raw"
folder = '/home/kolodziej/data_local/';
filename = 'alignment_info_from_homo_1_300s_lowstat.mat';
%filename = 'alignment_info_from_homo_2_1800s_highstat.mat';
%filename = 'alignment_info_high_stat_fdg_data.mat';
coincs = load([folder filename],'alignment_info');
dataset=coincs.alignment_info;

%% merge T_diff values in the main dataset
% Sort pairs and flip signs if order changed
for i = 1:length(dataset)
    originalPair = dataset{i,1};
    sortedPair = sort(dataset{i,1});
    
    % Check if pair was flipped
    if ~isequal(originalPair, sortedPair)
        dataset{i,2} = -dataset{i,2};  % Multiply list by -1
    end
    
    dataset{i,1} = sortedPair;  % Update pair to sorted version
end

% Get unique pairs and merge
[uniquePairs, ~, idx] = unique(cell2mat(dataset(:,1)), 'rows');

% Initialize result
result = cell(size(uniquePairs, 1), 3);

% Process each unique pair
for i = 1:size(uniquePairs, 1)
    % Find all rows with this pair
    matchRows = find(idx == i);
    
    % Store the pair
    result{i, 1} = uniquePairs(i, :);
    
    % Merge all lists from column 2
    mergedList = [];
    for j = matchRows'
        mergedList = [mergedList, dataset{j, 2}];
    end
    result{i, 2} = mergedList;
    
    % Set third column to 0
    result{i, 3} = mean(result{i, 2});

    % Store first number from pair in column 2
    result{i, 4} = uniquePairs(i, 1);
    
    % Store second number from pair in column 3
    result{i, 5} = uniquePairs(i, 2);

end


% 
% 
% len_dataset=length(dataset);
% 
% all_values = cell2mat(dataset(:,1));
% flat_values = all_values(:);
% unique_values = unique(flat_values);
% 
% 
% %if there are non-unique pairs, show them
% % Normalize each pair so [id1, id2] == [id2, id1]
% sorted_pairs = sort(all_values, 2);  % Each row now has smaller ID first
% 
% % Convert to strings for easy comparison
% pair_str = strcat(string(sorted_pairs(:,1)), "_", string(sorted_pairs(:,2)));
% 
% % Use tabulation to count occurrences
% [unique_keys, ~, idx_in_unique] = unique(pair_str);
% counts = histcounts(idx_in_unique, 1:max(idx_in_unique)+1);
% 
% % Find pairs that occur more than once
% duplicate_keys = unique_keys(counts > 1);
% 
% % Display duplicates and their indices
% %fprintf('Duplicate unordered pairs and their indices:\n');
% for i = 1:length(duplicate_keys)
%     key = duplicate_keys(i);
%     dup_indices = find(pair_str == key);
% 
%     % Recover the actual pair
%     pair = sorted_pairs(dup_indices(1), :);
% 
%     %fprintf('Pair [%d, %d] appears at indices: %s\n', ...
%     %    pair(1), pair(2), mat2str(dup_indices'));
% end
% 
% 
% 
% 
% sums=[];
% counters=[];
% means=[];
% for j = 1:length(unique_values)
%     sum=0;
%     counter=0;
% 
%     for i = 1:len_dataset
%         if (dataset{i,1}(1) == unique_values(j) || dataset{i,1}(2) == unique_values(j))
%             %disp(i+" "+dataset{i,1}(1)+" "+dataset{i,1}(2)+" "+dataset{i,3})
%             sum=sum+dataset{i,3};
%             counter=counter+1;
%             %disp(counter);
%         end
%     end
% 
%     if counter ~= 0
%         mean=sum/counter;
%     else
%         mean=0;
%     end
%     sums(end + 1) = sum;
%     counters(end + 1) = counter;
%     means(end + 1) = mean;
%     %fprintf('HERE %d;  %d; %d\n', sum, counter, mean);
%     %val = unique_values(j);
%     %fprintf('Element %d = %d\n', j, val);
% end
% histogram(means, 50);
% 
% 
%All channels vs all channels:
pairs = cell2mat(dataset(:,1));  % N x 2 matrix of ID pairs
values = cell2mat(dataset(:,3)); % N x 1 vector of values
ids = unique(pairs);          % All unique IDs
% 
% % Create ID to index map
% [~, ~, idxMap] = unique(ids); % idxMap(i) = index of ids(i)
% 
% % Initialize matrix
% n = length(ids);
% H = nan(n);  % Use NaN so unassigned elements appear blank in heatmap
% 
% % Fill in the matrix
% for i = 1:size(pairs,1)
%     id1 = pairs(i,1);
%     id2 = pairs(i,2);
%     val = values(i);
% 
%     idx1 = find(ids == id1);
%     idx2 = find(ids == id2);
% 
%     H(idx1, idx2) = val;
%     % Optionally: make symmetric if undirected
%     % H(idx2, idx1) = val;
% end

% % Plot
% figure;
% imagesc(H);
% xlim([255 512]);
% ylim([0 255]);
% colorbar;
% axis equal tight;
% 
% tick_step = 8;
% tick_indices = 1:tick_step:n;
% tick_labels = ids(tick_indices);
% 
% set(gca, 'XTick', tick_indices, 'XTickLabel', tick_labels, ...
%          'YTick', tick_indices, 'YTickLabel', tick_labels);
% xlabel('ID 2');
% ylabel('ID 1');
% title('2D Histogram of ID Pairs');
% % 
% % 
% % set(gca, 'XTick', 1:n, 'XTickLabel', ids, 'YTick', 1:n, 'YTickLabel', ids);
% % xlabel('ID 2');
% % ylabel('ID 1');
% % title('2D Histogram of ID Pairs');

%% read in mapping table:
mapping = readmatrix('/home/kolodziej/software/myMatlabScripts/mapping.txt');

% %% remap:
% % plotData = [ids, transpose(means)];
% plotData = [ids, cell2mat(result(:,3))]; %%this is incorrect!!!!! definition of means unclear, if needed, they have to be calculated again.
% [~, idxA] = ismember(plotData(:,1), mapping(:,1));  % match data IDs to mapping IDs
% new_columnA = mapping(idxA, 2);
% plotData = [plotData, new_columnA];
% 
% %% plotting per pixel
% % Create SiPM mapping
% sim_SiPM1 = flipud(reshape(0:63, 8, 8));
% sim_SiPM2 = sim_SiPM1 + 64;
% sim_bothSiPMs_module0 = [sim_SiPM2; sim_SiPM1];
% sim_bothSiPMs_module1 = [sim_SiPM2; sim_SiPM1] + 128;
% sim_bothSiPMs_module2 = [sim_SiPM2; sim_SiPM1] + 2*128;
% sim_bothSiPMs_module3 = [sim_SiPM2; sim_SiPM1] + 3*128;
% 
% if IDtype == "remapped"
%     id_to_value = containers.Map(plotData(:, 3), plotData(:, 2));
% elseif IDtype == "raw"
%     id_to_value = containers.Map(plotData(:, 1), plotData(:, 2));
% else 
%     disp("Error - IDtype undefined");
% end
% 
% % Verify each module is 16x8
% modules = {sim_bothSiPMs_module1, sim_bothSiPMs_module0, ...
%            sim_bothSiPMs_module2, sim_bothSiPMs_module3};
% module_names = {'Module 1', 'Module 0', 'Module 2', 'Module 3'};
% 
% figure;
% for idx = 1:4
%     % Get module and rotate 90 degrees counterclockwise
%     module = rot90(modules{idx});
%     [rows, cols] = size(module);
% 
%     % Create color values from table lookup
%     color_vals = zeros(rows, cols);
%     for i = 1:rows
%         for j = 1:cols
%             pixel_id = module(i, j);
%             color_vals(i, j) = id_to_value(pixel_id);
%         end
%     end
% 
%     subplot(2, 2, idx);
%     imagesc(color_vals);
%     colormap(parula);
%     colorbar;
% 
%     hold on;
% 
%     % Add pixel ID labels
%     for i = 1:rows
%         for j = 1:cols
%             pixel_id = module(i, j);
%             text(j, i, num2str(pixel_id), 'HorizontalAlignment', 'center', ...
%                 'VerticalAlignment', 'middle', 'FontSize', 6, 'Color', 'white', ...
%                 'FontWeight', 'bold');
%         end
%     end
% 
%     title(module_names{idx}, 'FontSize', 12, 'FontWeight', 'bold');
%     axis equal;
%     xlim([0.5, cols+0.5]);
%     ylim([0.5, rows+0.5]);
%     set(gca, 'YDir', 'normal');
% end
% 
% sgtitle('Counts per pixel', 'FontSize', 14, 'FontWeight', 'bold');

%% plotting per LOR


result_matrix = cell2mat(result(:,1));
result_matrix = [result_matrix, cell2mat(result(:,3))];

%% remap
[~, idxB] = ismember(result_matrix(:,1), mapping(:,1));  % match data IDs to mapping IDs
new_columnB = mapping(idxB, 2);
result_matrix = [result_matrix, new_columnB];

[~, idxC] = ismember(result_matrix(:,2), mapping(:,1));  % match data IDs to mapping IDs
new_columnC = mapping(idxC, 2);
result_matrix = [result_matrix, new_columnC];


if IDtype == "remapped"
    x1 = result_matrix(:,4);
    y1 = result_matrix(:,5);
    z1 = result_matrix(:,3); %counts
elseif IDtype == "raw"
    x1 = result_matrix(:,1);
    y1 = result_matrix(:,2);
    z1 = result_matrix(:,3); %counts
else 
    disp("Error - IDtype undefined");
end

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
