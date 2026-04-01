% %function [setOfEquivalentLORs] = createSetOfEquivalentLORs()
% clc; clear;
% 
% % Geometry parameters
% modules = 4;
% axCount = 8;
% transCount = 16;
% crystalsPerModule = axCount * transCount;
% N = modules * crystalsPerModule;
% 
% % Index mapping utility
% toID = @(m,a,t) (m-1)*crystalsPerModule + a*transCount + t;
% mirrorT = @(t) transCount-1-t;
% 
% % --- Compute crystal positions for plotting ---
% R_ring = 200; % mm
% moduleAnglesDeg = [-60, 60, 120, 240];
% pitchT = 5.0; pitchA = 5.0;
% halfT = (transCount-1)/2;
% halfA = (axCount-1)/2;
% pos = zeros(N,3);
% modIdx = zeros(N,1); axIdx = zeros(N,1); trIdx = zeros(N,1);
% for m = 1:modules
%     phi = deg2rad(moduleAnglesDeg(m));
%     cen = [R_ring*cos(phi), R_ring*sin(phi), 0];
%     tHat = [-sin(phi), cos(phi), 0];
%     for a = 0:axCount-1
%         z = (a - halfA) * pitchA;
%         for t = 0:transCount-1
%             xyt = cen + tHat * ((t - halfT)*pitchT);
%             id = toID(m,a,t);
%             pos(id+1,:) = [xyt(1), xyt(2), z];
%             modIdx(id+1) = m;
%             axIdx(id+1) = a;
%             trIdx(id+1) = t;
%         end
%     end
% end
% 
% % --- Build S_ij sets using full 8-LOR symmetry ---
% Sij_cell = {}; Sij_keys = {};
% 
% % --- Build ALL LORs (brute force) only between 1↔3 and 2↔4 ---
% idsM = cell(4,1);
% for m = 1:4
%     idsM{m} = find(modIdx==m) - 1;  % zero-based IDs for each module
% end
% 
% pairs = [];
% [IA, IB] = ndgrid(idsM{1}, idsM{3}); pairs = [pairs; IA(:) IB(:)];  % 1↔3
% [IA, IB] = ndgrid(idsM{2}, idsM{4}); pairs = [pairs; IA(:) IB(:)];  % 2↔4
% allLORs = pairs;  % (2*128*128) = 32768 rows
% 
% % --- Group by geometric equivalence (Kinouchi-like), invariant to swapping ---
% % Use mirrored trans index on the "right" module and take absolute diffs
% % so (aL,aR,tL,tR) and (aR,aL,tR,tL) fall in the same set.
% S_map = containers.Map('KeyType','char','ValueType','any');
% 
% for r = 1:size(allLORs,1)
%     i = allLORs(r,1); j = allLORs(r,2);      % i in {0..255}, j in {256..511} or {128..255,384..511}
%     aL = axIdx(i+1); tL = trIdx(i+1);
%     aR = axIdx(j+1); tR = trIdx(j+1);
% 
%     tR_m = mirrorT(tR);                 % align frames across opposing modules
%     dT   = tR_m - tL;
%     dA   = aR   - aL;
% 
%     % geometry-invariant key (treat swapping L/R as same geometry)
%     key  = sprintf('dT=%02d_dA=%02d', abs(dT), abs(dA));
% 
%     if ~isKey(S_map, key)
%         S_map(key) = zeros(0,2,'uint32');
%     end
%     S_map(key) = [S_map(key); uint32([i j])];
% end
% 
% % Convert map -> cell array (Sij_cell_new) and keep a parallel keys list
% Sij_keys = S_map.keys;
% Sij_cell = cell(numel(Sij_keys),1);
% for k = 1:numel(Sij_keys)
%     pairs_k = S_map(Sij_keys{k});
%     Sij_cell{k} = sortrows(double(pairs_k), [1 2]);  % keep your pair ordering style
% end
% 
% % If you still want a de-duplicated cell set in the exact way you had:
% strRepr = cellfun(@(x) mat2str(sortrows(x)), Sij_cell, 'UniformOutput', false);
% [uniqueStr, uniqueIdx] = unique(strRepr, 'stable'); %#ok<ASGLU>
% Sij_cell_unique = Sij_cell(uniqueIdx);
% 
% %%works, but without cross-axial LORs
% % for a = 0:axCount-1
% %     for tL = 0:transCount-1
% %         for tR = 0:transCount-1
% %             % Generate the 8 equivalent LORs for this combination
% %             combos = [tL, tR; mirrorT(tL), mirrorT(tR); tR, tL; mirrorT(tR), mirrorT(tL)];
% %             LORs = [];
% %             for k = 1:4
% %                 % Generate pairs for module pair 1-3
% %                 pair13 = [toID(1,a,combos(k,1)), toID(3,a,combos(k,2))];
% %                 % Generate pairs for module pair 2-4
% %                 pair24 = [toID(2,a,combos(k,1)), toID(4,a,combos(k,2))];
% %                 LORs_13(k,:) = pair13;  % store in array
% %                 LORs_24(k,:) = pair24;
% %             end
% % 
% %             % Remove duplicates within each module pair's set only
% %             LORs_13 = unique(LORs_13, 'rows');
% %             LORs_24 = unique(LORs_24, 'rows');
% % 
% %             % Concatenate both, no cross-duplication removal (because IDs different)
% %             LORs = [LORs_13; LORs_24];
% %             %disp(LORs);
% %             % Canonical key (for this axial/transaxial combo)
% %             key = sprintf('A%d_TL%d_TR%d', a, tL, tR);
% %             Sij_cell{end+1} = LORs; %#ok<AGROW>
% %             Sij_keys{end+1} = key; %#ok<AGROW>
% %         end
% %     end
% % end
% 
% for aL = 0:axCount-1
%     for aR = 0:axCount-1
%         for tL = 0:transCount-1
%             for tR = 0:transCount-1
%                 combos = [tL, tR;
%                           mirrorT(tL), mirrorT(tR);
%                           tR, tL;
%                           mirrorT(tR), mirrorT(tL)];
%                 LORs = [];
%                 for k = 1:4
%                     % Pairs for module pair 1-3
%                     pair13 = [toID(1,aL,combos(k,1)), toID(3,aR,combos(k,2))];
%                     % Pairs for module pair 2-4
%                     pair24 = [toID(2,aL,combos(k,1)), toID(4,aR,combos(k,2))];
%                     LORs_13(k,:) = pair13;
%                     LORs_24(k,:) = pair24;
%                 end
% 
%                 LORs_13 = unique(LORs_13,'rows');
%                 LORs_24 = unique(LORs_24,'rows');
%                 LORs    = [LORs_13; LORs_24];
% 
%                 key = sprintf('AL%d_AR%d_TL%d_TR%d', aL, aR, tL, tR);
%                 Sij_cell{end+1} = sortrows(LORs); %#ok<AGROW>
%                 Sij_keys{end+1} = key; %#ok<AGROW>
%             end
%         end
%     end
% end
% 
% %disp(Sij_cell{1});
% % --- Preprocess sets: sort and unique each S_ij set's crystal pairs ---
% for s = 1:numel(Sij_cell)
%      pairs = Sij_cell{s};
% %     % Sort pairs by first column ascending, then by second column ascending
%      pairs = sortrows(pairs, [1, 2]);
% %     disp("sizeBEFORE"+size(pairs));
% %     % Remove duplicate rows
% %     pairs = unique(pairs, 'rows');
% %     disp("size"+size(pairs));
% %     %disp(pairs);
% %     % Save back sorted and uniqued pairs
%      Sij_cell{s} = pairs;
% end
% 
% % Convert each cell to a string representation for comparison
% strRepr = cellfun(@(x) mat2str(sortrows(x)), Sij_cell, 'UniformOutput', false);
% 
% % Find unique string representations and their indices
% [uniqueStr, uniqueIdx] = unique(strRepr, 'stable');
% 
% % Get the unique cells based on content
% Sij_cell_unique = Sij_cell(uniqueIdx);
% 
% % Now C_unique is the cell array with duplicate 8x2 double cells removed
% fprintf('Original cells: %d\nUnique cells: %d\n', numel(Sij_cell), numel(Sij_cell_unique));
% 
% 
% % S_ij_array = cell2mat(Sij_cell);
% % S_ij_array = sortrows(S_ij_array, [1, 2]);
% % S_ij_array = unique(S_ij_array, 'rows');
% 
% % % --- Output each equivalence set as one line in TXT ---
% % fileID = fopen('/home/kolodziej/software/myMatlabScripts/Sij_sets.txt','w');
% % for s = 1:numel(Sij_cell_unique)
% %     %fprintf(fileID, '%s\t', Sij_keys{s});
% %     fprintf(fileID, '\t');
% %     pairs = Sij_cell_unique{s};
% %     for p = 1:size(pairs,1)
% %         fprintf(fileID, '%d %d\t', pairs(p,1), pairs(p,2));
% %     end
% %     fprintf(fileID, '\n');
% % end
% % fclose(fileID);
% % disp(['Wrote Sij_sets.txt with ' num2str(numel(Sij_cell)) ' equivalence sets.']);
% 
% % --- Output each equivalence set as one line in TXT ---
% fileID = fopen('/home/kolodziej/software/myMatlabScripts/Sij_sets.txt','w');
% for s = 1:numel(Sij_cell_unique)
%     fprintf(fileID, '\t');
%     pairs = Sij_cell_unique{s};
%     for p = 1:size(pairs,1)
%         fprintf(fileID, '%d %d\t', pairs(p,1), pairs(p,2));
%     end
%     fprintf(fileID, '\n');
% end
% fclose(fileID);
% disp(['Wrote Sij_sets.txt with ' num2str(numel(Sij_cell_unique)) ' equivalence sets.']);
% 
% % --- Visualization: pick a valid setIdx and plot that S_ij ---
% % setIdx = min(10, numel(Sij_cell_unique));   % pick something that exists
% setIdx = 2048;
% pairsToPlot = Sij_cell_unique{setIdx};
% 
% figure('Color','w'); hold on; grid on;
% title(sprintf('Crystals and S_{ij} set #%d', setIdx));
% plot3(pos(:,1), pos(:,2), pos(:,3), '.', 'MarkerSize', 8);
% 
% for r = 1:size(pairsToPlot,1)
%     i = pairsToPlot(r,1)+1; j = pairsToPlot(r,2)+1;
%     line([pos(i,1), pos(j,1)], [pos(i,2), pos(j,2)], [pos(i,3), pos(j,3)], ...
%          'LineWidth', 2, 'LineStyle','-');
% end
% 
% shownIDs = unique(pairsToPlot(:)) + 1;
% for idx = reshape(shownIDs,1,[])
%     text(pos(idx,1), pos(idx,2), pos(idx,3), sprintf(' %d', idx-1), ...
%          'FontSize', 7, 'Color', [0 0 0]);
% end
% 
% xlabel('x [mm]'); ylabel('y [mm]'); zlabel('z [mm]');
% axis equal; view(35,20);
% legend({'Crystals','S_{ij} set'}, 'Location','best'); hold off;
% 
% % % --- Visualization: plot setIdx-th S_ij set (change setIdx as desired) ---
% % setIdx = 2048; % Example: plot the 3rd set, which now is guaranteed 8 LORs
% % figure('Color','w'); hold on; grid on;
% % title(sprintf('Crystals and S_{ij} set #%d', setIdx));
% % 
% % 
% % % Plot all crystals
% % plot3(pos(:,1), pos(:,2), pos(:,3), '.', 'MarkerSize', 8);
% % 
% % % Draw the LORs for current S_ij set
% % pairsToPlot = Sij_cell{setIdx};
% % for r = 1:size(pairsToPlot,1)
% %     i = pairsToPlot(r,1)+1; j = pairsToPlot(r,2)+1;
% %     line([pos(i,1), pos(j,1)], [pos(i,2), pos(j,2)], [pos(i,3), pos(j,3)], ...
% %          'LineWidth', 2, 'LineStyle','-');
% % end
% % 
% % % Label only crystals involved in these LORs
% % shownIDs = unique(pairsToPlot(:)) + 1;
% % for idx = reshape(shownIDs,1,[])
% %     text(pos(idx,1), pos(idx,2), pos(idx,3), sprintf(' %d', idx-1), ...
% %          'FontSize', 7, 'Color', [0 0 0]);
% % end
% % 
% % xlabel('x [mm]'); ylabel('y [mm]'); zlabel('z [mm]');
% % axis equal; view(35,20);
% % legend({'Crystals','S_{ij} set'}, 'Location','best');
% % hold off;
% 
% % % --- Build lookup table of LOR equivalence sets ---
% % Assume cellArray is 1x576, each cell contains [Nx,2] double (variable Nx)
% 
% pairs_list = []; % Will collect: [firstElement, secondElement, cellIndex]
% 
% for cellidx = 1:numel(Sij_cell_unique)
%     pairs = Sij_cell_unique{cellidx}; % [Nx,2]
%     numPairs = size(pairs,1);
%     % Create [firstElement, secondElement, cellIndex] rows
%     if numPairs > 0
%         % Append [pair(:,1), pair(:,2), repmat(cellidx,numPairs,1)]
%         pairs_list = [pairs_list; [pairs, repmat(cellidx, numPairs, 1)]];
%     end
% end
% 
% for s = 1:numel(pairs_list)
%      tmp_pairs = pairs_list;
%      tmp_pairs = sortrows(tmp_pairs, [1, 2]);
%      pairs_list = tmp_pairs;
% end


%function [setOfEquivalentLORs] = createSetOfEquivalentLORs()
clc; clear;

% Geometry parameters
modules = 4;
axCount = 8;
transCount = 16;
crystalsPerModule = axCount * transCount;
N = modules * crystalsPerModule;

% Index mapping utility
toID = @(m,a,t) (m-1)*crystalsPerModule + a*transCount + t;
mirrorT = @(t) transCount-1-t;

% --- Compute crystal positions for plotting ---
R_ring = 200; % mm
moduleAnglesDeg = [-60, 60, 120, 240];
pitchT = 5.0; pitchA = 5.0;
halfT = (transCount-1)/2;
halfA = (axCount-1)/2;
pos = zeros(N,3);
modIdx = zeros(N,1); axIdx = zeros(N,1); trIdx = zeros(N,1);
for m = 1:modules
    phi = deg2rad(moduleAnglesDeg(m));
    cen = [R_ring*cos(phi), R_ring*sin(phi), 0];
    tHat = [-sin(phi), cos(phi), 0];
    for a = 0:axCount-1
        z = (a - halfA) * pitchA;
        for t = 0:transCount-1
            xyt = cen + tHat * ((t - halfT)*pitchT);
            id = toID(m,a,t);
            pos(id+1,:) = [xyt(1), xyt(2), z];
            modIdx(id+1) = m;
            axIdx(id+1) = a;
            trIdx(id+1) = t;
        end
    end
end

% --- Build S_ij sets using full 8-LOR symmetry ---
Sij_cell = {}; Sij_keys = {};

% --- Build S_ij sets (brute-force left->right pairs, canonicalized) ---
% Left modules: 1 and 2. Right modules: 3 and 4.
leftIDs  = find(modIdx==1 | modIdx==2) - 1;   % zero-based
rightIDs = find(modIdx==3 | modIdx==4) - 1;   % zero-based

% Map from canonical geometry key -> list of [i j] pairs (zero-based IDs)
S_map = containers.Map('KeyType','char','ValueType','any');

for ii = 1:numel(leftIDs)
    i = leftIDs(ii);
    ai = axIdx(i+1); ti = trIdx(i+1);
    for jj = 1:numel(rightIDs)
        j = rightIDs(jj);
        aj = axIdx(j+1); tj = trIdx(j+1);

        % Mirror the right-module trans index to align local frames
        tj_m = mirrorT(tj);

        % Relative displacement in the aligned frame (right minus left)
        dT = tj_m - ti;
        dA = aj   - ai;

        % Canonicalize so that (dT,dA) and (-dT,-dA) map to same key:
        % choose the representation with dA>0, or if dA==0 then dT>=0.
        if (dA < 0) || (dA == 0 && dT < 0)
            dT = -dT;
            dA = -dA;
        end

        % Build key (zero-pad for readability / stability)
        key = sprintf('dT%03d_dA%03d', dT, dA);

        if ~isKey(S_map, key)
            S_map(key) = zeros(0,2);  % initialize empty Nx2 array
        end
        % Append this left->right pair (keep zero-based IDs as in your code)
        S_map(key) = [S_map(key); double(i), double(j)];
    end
end

% Convert map -> cell arrays (sorted, unique rows inside each set)
Sij_keys = S_map.keys;
Sij_cell = cell(numel(Sij_keys),1);
for k = 1:numel(Sij_keys)
    pairs = S_map(Sij_keys{k});
    pairs = unique(pairs,'rows');        % remove any accidental duplicates
    pairs = sortrows(pairs, [1 2]);      % consistent ordering
    Sij_cell{k} = pairs;
end

% If you want a variable named Sij_cell_unique (keeps your later code unchanged)
Sij_cell_unique = Sij_cell;


% %works, but without cross-axial LORs
% for a = 0:axCount-1
%     for tL = 0:transCount-1
%         for tR = 0:transCount-1
%             % Generate the 8 equivalent LORs for this combination
%             combos = [tL, tR; mirrorT(tL), mirrorT(tR); tR, tL; mirrorT(tR), mirrorT(tL)];
%             LORs = [];
%             for k = 1:4
%                 % Generate pairs for module pair 1-3
%                 pair13 = [toID(1,a,combos(k,1)), toID(3,a,combos(k,2))];
%                 % Generate pairs for module pair 2-4
%                 pair24 = [toID(2,a,combos(k,1)), toID(4,a,combos(k,2))];
%                 LORs_13(k,:) = pair13;  % store in array
%                 LORs_24(k,:) = pair24;
%             end
% 
%             % Remove duplicates within each module pair's set only
%             LORs_13 = unique(LORs_13, 'rows');
%             LORs_24 = unique(LORs_24, 'rows');
% 
%             % Concatenate both, no cross-duplication removal (because IDs different)
%             LORs = [LORs_13; LORs_24];
%             %disp(LORs);
%             % Canonical key (for this axial/transaxial combo)
%             key = sprintf('A%d_TL%d_TR%d', a, tL, tR);
%             Sij_cell{end+1} = LORs; %#ok<AGROW>
%             Sij_keys{end+1} = key; %#ok<AGROW>
%         end
%     end
% end

% for aL = 0:axCount-1     % left crystal axial index
%     for aR = 0:axCount-1 % right crystal axial index
%         for tL = 0:transCount-1
%             for tR = 0:transCount-1
%                 % generate symmetry-equivalent combinations using axial/transaxial indices
%                 combos = [tL, tR;
%                           mirrorT(tL), mirrorT(tR);
%                           tR, tL;
%                           mirrorT(tR), mirrorT(tL)];
%                 LORs = [];
%                 for k = 1:4
%                     % Both module pairs, both axial indices
%                     pair13 = [toID(1,aL,combos(k,1)), toID(3,aR,combos(k,2))];
%                     pair24 = [toID(2,aL,combos(k,1)), toID(4,aR,combos(k,2))];
%                     LORs_13(k,:) = pair13;
%                     LORs_24(k,:) = pair24;
%                 end
%                 % Remove duplicates AND include cross-module, cross-axial, cross-transaxial pairs
%                 LORs_13 = unique(LORs_13,'rows');
%                 LORs_24 = unique(LORs_24,'rows');
%                 LORs = [LORs_13; LORs_24];
%                 key = sprintf('AL%d_AR%d_TL%d_TR%d', aL, aR, tL, tR);
%                 Sij_cell{end+1} = sortrows(LORs);
%                 Sij_keys{end+1} = key;
%             end
%         end
%     end
% end

% for aL = 0:axCount-1
%     for aR = 0:axCount-1
%         for tL = 0:transCount-1
%             for tR = 0:transCount-1
%                 combos = [tL, tR;
%                           mirrorT(tL), mirrorT(tR);
%                           tR, tL;
%                           mirrorT(tR), mirrorT(tL)];
%                 LORs = [];
%                 for k = 1:4
%                     % Pairs for module pair 1-3
%                     pair13 = [toID(1,aL,combos(k,1)), toID(3,aR,combos(k,2))];
%                     % Pairs for module pair 2-4
%                     pair24 = [toID(2,aL,combos(k,1)), toID(4,aR,combos(k,2))];
%                     LORs_13(k,:) = pair13;
%                     LORs_24(k,:) = pair24;
%                 end
% 
%                 LORs_13 = unique(LORs_13,'rows');
%                 LORs_24 = unique(LORs_24,'rows');
%                 LORs    = [LORs_13; LORs_24];
% 
%                 key = sprintf('AL%d_AR%d_TL%d_TR%d', aL, aR, tL, tR);
%                 Sij_cell{end+1} = sortrows(LORs); %#ok<AGROW>
%                 Sij_keys{end+1} = key; %#ok<AGROW>
%             end
%         end
%     end
% end

% Sij_sets = {};
% Sij_keys = {};
% 
% for mL = 1:modules
%     for mR = 1:modules
%         if mL == mR, continue; end % skip same module
%         for aL = 0:axCount-1
%             for aR = 0:axCount-1
%                 for deltaT = -(transCount-1):(transCount-1)
%                     LORs = [];
%                     % Enumerate all possible tL, tR pairs with given deltaT
%                     for tL = 0:transCount-1
%                         tR = tL + deltaT;
%                         if tR < 0 || tR >= transCount, continue; end
%                         idL = toID(mL, aL, tL);
%                         idR = toID(mR, aR, tR);
%                         % Add both directions, mirrored versions too
%                         LORs = [LORs; idL, idR; idR, idL];
%                         % Mirrored
%                         tL_m = mirrorT(tL);
%                         tR_m = mirrorT(tR);
%                         idL_m = toID(mL, aL, tL_m);
%                         idR_m = toID(mR, aR, tR_m);
%                         LORs = [LORs; idL_m, idR_m; idR_m, idL_m];
%                     end
%                     % Remove duplicates
%                     LORs = unique(LORs, 'rows');
%                     if ~isempty(LORs)
%                         key = sprintf('m%d_m%d_a%d_a%d_dT%d', mL, mR, aL, aR, deltaT);
%                         Sij_sets{end+1} = LORs;
%                         Sij_keys{end+1} = key;
%                     end
%                 end
%             end
%         end
%     end
% end


%disp(Sij_cell{1});
% --- Preprocess sets: sort and unique each S_ij set's crystal pairs ---
for s = 1:numel(Sij_cell)
     pairs = Sij_cell{s};
%     % Sort pairs by first column ascending, then by second column ascending
     pairs = sortrows(pairs, [1, 2]);
%     disp("sizeBEFORE"+size(pairs));
%     % Remove duplicate rows
%     pairs = unique(pairs, 'rows');
%     disp("size"+size(pairs));
%     %disp(pairs);
%     % Save back sorted and uniqued pairs
     Sij_cell{s} = pairs;
end

% Convert each cell to a string representation for comparison
strRepr = cellfun(@(x) mat2str(sortrows(x)), Sij_cell, 'UniformOutput', false);

% Find unique string representations and their indices
[uniqueStr, uniqueIdx] = unique(strRepr, 'stable');

% Get the unique cells based on content
Sij_cell_unique = Sij_cell(uniqueIdx);

% Now C_unique is the cell array with duplicate 8x2 double cells removed
fprintf('Original cells: %d\nUnique cells: %d\n', numel(Sij_cell), numel(Sij_cell_unique));


% S_ij_array = cell2mat(Sij_cell);
% S_ij_array = sortrows(S_ij_array, [1, 2]);
% S_ij_array = unique(S_ij_array, 'rows');

% --- Output each equivalence set as one line in TXT ---
fileID = fopen('/home/kolodziej/software/myMatlabScripts/Sij_sets.txt','w');
for s = 1:numel(Sij_cell_unique)
    %fprintf(fileID, '%s\t', Sij_keys{s});
    fprintf(fileID, '\t');
    pairs = Sij_cell_unique{s};
    for p = 1:size(pairs,1)
        fprintf(fileID, '%d %d\t', pairs(p,1), pairs(p,2));
    end
    fprintf(fileID, '\n');
end
fclose(fileID);
disp(['Wrote Sij_sets.txt with ' num2str(numel(Sij_cell)) ' equivalence sets.']);


% % --- Visualization: plot setIdx-th S_ij set (change setIdx as desired) ---
% setIdx = 3; % Example: plot the 3rd set, which now is guaranteed 8 LORs
% figure('Color','w'); hold on; grid on;
% title(sprintf('Crystals and S_{ij} set #%d', setIdx));
% 
% % Plot all crystals
% plot3(pos(:,1), pos(:,2), pos(:,3), '.', 'MarkerSize', 8);
% 
% % Draw the LORs for current S_ij set
% pairsToPlot = Sij_cell{setIdx};
% for r = 1:size(pairsToPlot,1)
%     i = pairsToPlot(r,1)+1; j = pairsToPlot(r,2)+1;
%     line([pos(i,1), pos(j,1)], [pos(i,2), pos(j,2)], [pos(i,3), pos(j,3)], ...
%          'LineWidth', 2, 'LineStyle','-');
% end
% 
% % Label only crystals involved in these LORs
% shownIDs = unique(pairsToPlot(:)) + 1;
% for idx = reshape(shownIDs,1,[])
%     text(pos(idx,1), pos(idx,2), pos(idx,3), sprintf(' %d', idx-1), ...
%          'FontSize', 7, 'Color', [0 0 0]);
% end
% 
% xlabel('x [mm]'); ylabel('y [mm]'); zlabel('z [mm]');
% axis equal; view(35,20);
% legend({'Crystals','S_{ij} set'}, 'Location','best');
% hold off;

% Example: Plot S_ij set with arbitrary index (setIdx)
setIdx = 649; % Choose the set index to plot

figure('Color','w'); hold on; grid on;
title(sprintf('Crystals and S_{ij} set #%d (%d LORs)', setIdx, size(Sij_cell_unique{setIdx},1)));

% Plot all crystals (gray dots)
plot3(pos(:,1), pos(:,2), pos(:,3), '.', 'MarkerSize', 8, 'Color', [0.7 0.7 0.7]);

% Draw all LORs for the current S_ij set
pairsToPlot = Sij_cell_unique{setIdx};
for r = 1:size(pairsToPlot, 1)
    i = pairsToPlot(r,1) + 1; % MATLAB 1-based indexing
    j = pairsToPlot(r,2) + 1;
    line([pos(i,1), pos(j,1)], [pos(i,2), pos(j,2)], [pos(i,3), pos(j,3)], ...
         'LineWidth', 2, 'Color', [0 0.5 1]);
end

% Highlight and label only crystals involved in these LORs (blue circles)
shownIDs = unique(pairsToPlot(:)) + 1;
plot3(pos(shownIDs,1), pos(shownIDs,2), pos(shownIDs,3), 'o', 'MarkerSize', 10, 'LineWidth', 2, 'Color', [0 0.3 1]);
for idx = reshape(shownIDs, 1, [])
    text(pos(idx,1), pos(idx,2), pos(idx,3), sprintf(' %d', idx-1), ...
         'FontSize', 7, 'Color', [0 0.28 0.8]);
end

xlabel('x [mm]'); ylabel('y [mm]'); zlabel('z [mm]');
axis equal; view(35,20); grid on;
legend({'All crystals','S_{ij} set LORs','S_{ij} crystals'}, 'Location','best');
hold off;


% % --- Build lookup table of LOR equivalence sets ---
% Assume cellArray is 1x576, each cell contains [Nx,2] double (variable Nx)

pairs_list = []; % Will collect: [firstElement, secondElement, cellIndex]

for cellidx = 1:numel(Sij_cell_unique)
    pairs = Sij_cell_unique{cellidx}; % [Nx,2]
    numPairs = size(pairs,1);
    % Create [firstElement, secondElement, cellIndex] rows
    if numPairs > 0
        % Append [pair(:,1), pair(:,2), repmat(cellidx,numPairs,1)]
        pairs_list = [pairs_list; [pairs, repmat(cellidx, numPairs, 1)]];
    end
end

for s = 1:numel(pairs_list)
     tmp_pairs = pairs_list;
     tmp_pairs = sortrows(tmp_pairs, [1, 2]);
     pairs_list = tmp_pairs;
end

% % --- Build lookup table of LOR equivalence sets --- %%NOT VALIDATED - to
% % be checked if working correctly!!! Still to be remapped because indices
% % here not correct.
% LOR_lookup = cell(N,N); % Preallocate 256x256 cell array
% %LOR_lookup_size = array(N,N);
% % Loop through each unique S_ij set
% for s = 1:numel(Sij_cell_unique)
%     pairs = Sij_cell_unique{s};
%     for p = 1:size(pairs,1)
%         i = pairs(p,1);
%         j = pairs(p,2);
%         % Store all other pairs in the same set (excluding the current one)
%         otherPairs = pairs(~ismember(pairs, [i, j], 'rows'), :);
%         % Append to the cell at (i+1, j+1)
%         if isempty(LOR_lookup{i+1, j+1})
%             LOR_lookup{i+1, j+1} = otherPairs;
%         else
%             LOR_lookup{i+1, j+1} = [LOR_lookup{i+1, j+1}; otherPairs];
%         end
%         LOR_lookup_size{i+1, j+1} = length(LOR_lookup{i+1, j+1});
%     end
% end

% Optional: remove duplicates in each cell
% for i = 1:N
%     for j = 1:N
%         if ~isempty(LOR_lookup{i,j})
%             LOR_lookup{i,j} = unique(LOR_lookup{i,j}, 'rows');
%         end
%     end
% end

%disp('Generated LOR equivalence lookup table.');

%setOfEquivalentLORs = LOR_lookup;