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
% % --- Build S_ij sets (brute-force left->right pairs, canonicalized) ---
% % Left modules: 1 and 2. Right modules: 3 and 4.
% leftIDs  = find(modIdx==1 | modIdx==2) - 1;   % zero-based
% rightIDs = find(modIdx==3 | modIdx==4) - 1;   % zero-based
% 
% % Map from canonical geometry key -> list of [i j] pairs (zero-based IDs)
% S_map = containers.Map('KeyType','char','ValueType','any');
% 
% for ii = 1:numel(leftIDs)
%     i = leftIDs(ii);
%     ai = axIdx(i+1); ti = trIdx(i+1);
%     for jj = 1:numel(rightIDs)
%         j = rightIDs(jj);
%         aj = axIdx(j+1); tj = trIdx(j+1);
% 
%         % Mirror the right-module trans index to align local frames
%         tj_m = mirrorT(tj);
% 
%         % Relative displacement in the aligned frame (right minus left)
%         dT = tj_m - ti;
%         dA = aj   - ai;
% 
%         % Canonicalize so that (dT,dA) and (-dT,-dA) map to same key:
%         % choose the representation with dA>0, or if dA==0 then dT>=0.
%         if (dA < 0) || (dA == 0 && dT < 0)
%             dT = -dT;
%             dA = -dA;
%         end
% 
%         % Build key (zero-pad for readability / stability)
%         key = sprintf('dT%03d_dA%03d', dT, dA);
% 
%         if ~isKey(S_map, key)
%             S_map(key) = zeros(0,2);  % initialize empty Nx2 array
%         end
%         % Append this left->right pair (keep zero-based IDs as in your code)
%         S_map(key) = [S_map(key); double(i), double(j)];
%     end
% end
% 
% % Convert map -> cell arrays (sorted, unique rows inside each set)
% Sij_keys = S_map.keys;
% Sij_cell = cell(numel(Sij_keys),1);
% for k = 1:numel(Sij_keys)
%     pairs = S_map(Sij_keys{k});
%     pairs = unique(pairs,'rows');        % remove any accidental duplicates
%     pairs = sortrows(pairs, [1 2]);      % consistent ordering
%     Sij_cell{k} = pairs;
% end
% 
% % If you want a variable named Sij_cell_unique (keeps your later code unchanged)
% Sij_cell_unique = Sij_cell;
% 
% 
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
% 
% 
% 
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
% 
% 
% % --- Output each equivalence set as one line in TXT ---
% fileID = fopen('/home/kolodziej/software/myMatlabScripts/Sij_sets.txt','w');
% for s = 1:numel(Sij_cell_unique)
%     %fprintf(fileID, '%s\t', Sij_keys{s});
%     fprintf(fileID, '\t');
%     pairs = Sij_cell_unique{s};
%     for p = 1:size(pairs,1)
%         fprintf(fileID, '%d %d\t', pairs(p,1), pairs(p,2));
%     end
%     fprintf(fileID, '\n');
% end
% fclose(fileID);
% disp(['Wrote Sij_sets.txt with ' num2str(numel(Sij_cell)) ' equivalence sets.']);
% 
% 
% 
% % Example: Plot S_ij set with arbitrary index (setIdx)
% setIdx = 649; % Choose the set index to plot
% 
% figure('Color','w'); hold on; grid on;
% title(sprintf('Crystals and S_{ij} set #%d (%d LORs)', setIdx, size(Sij_cell_unique{setIdx},1)));
% 
% % Plot all crystals (gray dots)
% plot3(pos(:,1), pos(:,2), pos(:,3), '.', 'MarkerSize', 8, 'Color', [0.7 0.7 0.7]);
% 
% % Draw all LORs for the current S_ij set
% pairsToPlot = Sij_cell_unique{setIdx};
% for r = 1:size(pairsToPlot, 1)
%     i = pairsToPlot(r,1) + 1; % MATLAB 1-based indexing
%     j = pairsToPlot(r,2) + 1;
%     line([pos(i,1), pos(j,1)], [pos(i,2), pos(j,2)], [pos(i,3), pos(j,3)], ...
%          'LineWidth', 2, 'Color', [0 0.5 1]);
% end
% 
% % Highlight and label only crystals involved in these LORs (blue circles)
% shownIDs = unique(pairsToPlot(:)) + 1;
% plot3(pos(shownIDs,1), pos(shownIDs,2), pos(shownIDs,3), 'o', 'MarkerSize', 10, 'LineWidth', 2, 'Color', [0 0.3 1]);
% for idx = reshape(shownIDs, 1, [])
%     text(pos(idx,1), pos(idx,2), pos(idx,3), sprintf(' %d', idx-1), ...
%          'FontSize', 7, 'Color', [0 0.28 0.8]);
% end
% 
% xlabel('x [mm]'); ylabel('y [mm]'); zlabel('z [mm]');
% axis equal; view(35,20); grid on;
% legend({'All crystals','S_{ij} set LORs','S_{ij} crystals'}, 'Location','best');
% hold off;
% 
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

clc; clear;

%% Geometry parameters
modules = 4;
axCount = 8;
transCount = 16;
crystalsPerModule = axCount * transCount; % 128
N = modules * crystalsPerModule;          % 512

% Index mapping utilities
toID = @(m,a,t) (m-1)*crystalsPerModule + a*transCount + t;
mirrorT = @(t) transCount-1-t;

%% Compute crystal positions (for plotting)
R_ring = 200; % mm
%moduleAnglesDeg = [-60, 60, 120, 240];
moduleAnglesDeg = [-15, +15, 165, 195];
pitchT = 5.0; pitchA = 5.0;
halfT = (transCount-1)/2;
halfA = (axCount-1)/2;
pos = zeros(N,3);
modIdx = zeros(N,1); axIdx = zeros(N,1); trIdx = zeros(N,1);

for m = 1:modules
    phi = deg2rad(moduleAnglesDeg(m));
    cen = [R_ring*cos(phi), R_ring*sin(phi), 0];
    tHat = [-sin(phi), cos(phi), 0]; % tangent direction
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

%% Build S_ij sets (brute force, canonicalized geometry)
% Left modules: 1 and 2. Right modules: 3 and 4.
% leftIDs  = find(modIdx==1 | modIdx==2) - 1;   % zero-based
% rightIDs = find(modIdx==3 | modIdx==4) - 1;   % zero-based
% Opposing module pairs only: 1↔3 and 2↔4
pairsOfModules = [1 3; 2 4];

S_map = containers.Map('KeyType','char','ValueType','any');
n = transCount; % number of transaxial crystals

for pp = 1:size(pairsOfModules,1)
    mL = pairsOfModules(pp,1);
    mR = pairsOfModules(pp,2);

    for aL = 0:axCount-1
        for aR = 0:axCount-1
            dA = abs(aR - aL); % absolute axial difference

            for tL = 0:transCount-1
                for tR = 0:transCount-1
                    % Build the 4 symmetry-related transaxial pairs
                    combos = [ tL,      tR;
                               tR,      tL;
                               (n-1-tL),(n-1-tR);
                               (n-1-tR),(n-1-tL) ];

                    for k = 1:size(combos,1)
                        tLi = combos(k,1);
                        tRi = combos(k,2);

                        % Canonicalize transaxial pair (order-independent)
                        tmin = min(tLi, tRi);
                        tmax = max(tLi, tRi);

                        % Build IDs
                        idL = toID(mL, aL, tLi);
                        idR = toID(mR, aR, tRi);

                        % Geometry key: |dA| plus canonicalized transaxial pair
                        key = sprintf('dA%02d_T%02d_%02d', dA, tmin, tmax);

                        if ~isKey(S_map,key)
                            S_map(key) = zeros(0,2);
                        end
                        S_map(key) = [S_map(key); idL, idR];
                    end
                end
            end
        end
    end
end

% Convert map -> cell arrays
Sij_keys = S_map.keys;
Sij_cell = cell(numel(Sij_keys),1);
for k = 1:numel(Sij_keys)
    pairs = unique(S_map(Sij_keys{k}),'rows');
    Sij_cell{k} = sortrows(pairs,[1 2]);
end
Sij_cell_unique = Sij_cell;

%juz spoko
% % Map from geometry key -> list of [i j] pairs
% S_map = containers.Map('KeyType','char','ValueType','any');
% 
% for pp = 1:size(pairsOfModules,1)
%     mL = pairsOfModules(pp,1);
%     mR = pairsOfModules(pp,2);
% 
%     for aL = 0:axCount-1
%         for aR = 0:axCount-1
%             dA = abs(aR - aL); % absolute axial difference
% 
%             for tL = 0:transCount-1
%                 for tR = 0:transCount-1
%                     % generate equivalent transaxial combos
%                     combos = [tL, tR;
%                               mirrorT(tL), mirrorT(tR);
%                               tR, tL;
%                               mirrorT(tR), mirrorT(tL)];
% 
%                     for k = 1:size(combos,1)
%                         tLi = combos(k,1);
%                         tRi = combos(k,2);
% 
%                         % Canonicalize transaxial pair (order-independent)
%                         tmin = min(tLi, tRi);
%                         tmax = max(tLi, tRi);
% 
%                         % Build IDs
%                         idL = toID(mL, aL, tLi);
%                         idR = toID(mR, aR, tRi);
% 
%                         % Geometry key: axial offset + unordered transaxial pair
%                         key = sprintf('dA%02d_T%02d_%02d', dA, tmin, tmax);
% 
%                         if ~isKey(S_map,key)
%                             S_map(key) = zeros(0,2);
%                         end
%                         S_map(key) = [S_map(key); idL, idR];
%                     end
%                 end
%             end
%         end
%     end
% end
% 
% % Convert to cell array
% Sij_keys = S_map.keys;
% Sij_cell = cell(numel(Sij_keys),1);
% for k = 1:numel(Sij_keys)
%     pairs = unique(S_map(Sij_keys{k}),'rows');
%     Sij_cell{k} = sortrows(pairs,[1 2]);
% end
% Sij_cell_unique = Sij_cell;

% %axials are merged into one set, but also transaxials which is bad
% for pp = 1:size(pairsOfModules,1)
%     mL = pairsOfModules(pp,1);
%     mR = pairsOfModules(pp,2);
% 
%     leftIDs  = find(modIdx==mL) - 1;   % IDs in module mL
%     rightIDs = find(modIdx==mR) - 1;   % IDs in module mR
% 
% 
%     for ii = 1:numel(leftIDs)
%         i = leftIDs(ii);
%         ai = axIdx(i+1); ti = trIdx(i+1);
%         for jj = 1:numel(rightIDs)
%             j = rightIDs(jj);
%             aj = axIdx(j+1); tj = trIdx(j+1);
% 
%             % Mirror right-module trans index to align local frames
%             tj_m = mirrorT(tj); %ok, chyba wynika z tego ze jest na kole
% 
%             % Relative displacement in aligned frame
%             dT = tj_m - ti;
%             dA = aj   - ai;
% 
%             % Canonicalize
%             if (dA < 0) || (dA == 0 && dT < 0)
%                 dT = -dT;
%                 dA = -dA;
%             end
% 
%             key = sprintf('dT%03d_dA%03d', dT, dA);
% 
%             if ~isKey(S_map,key)
%                 S_map(key) = zeros(0,2);
%             end
%             disp(key)
%             S_map(key) = [S_map(key); double(i), double(j)];
%         end
%     end
% end

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


% Convert map -> cell arrays
Sij_keys = S_map.keys;
Sij_cell = cell(numel(Sij_keys),1);
for k = 1:numel(Sij_keys)
    pairs = S_map(Sij_keys{k});
    pairs = unique(pairs,'rows');
    pairs = sortrows(pairs,[1 2]);
    Sij_cell{k} = pairs;
end

Sij_cell_unique = Sij_cell; % final set of equivalence classes

fprintf('Total unique S_ij sets: %d\n', numel(Sij_cell_unique));

%% Output each set as one line in TXT
% outFile = 'Sij_sets.txt';
% fileID = fopen(outFile,'w');
% for s = 1:numel(Sij_cell_unique)
%     pairs = Sij_cell_unique{s};
%     fprintf(fileID,'\t');
%     for p = 1:size(pairs,1)
%         fprintf(fileID,'%d %d\t',pairs(p,1),pairs(p,2));
%     end
%     fprintf(fileID,'\n');
% end
% fclose(fileID);
% disp(['Wrote ', outFile]);

outFile = 'Sij_sets.txt';
fileID = fopen(outFile,'w');
for s = 1:numel(Sij_cell_unique)
    pairs = Sij_cell_unique{s};
    key   = Sij_keys{s};
    fprintf(fileID,'%s\t', key); % print geometry key first
    for p = 1:size(pairs,1)
        fprintf(fileID,'%d %d\t', pairs(p,1), pairs(p,2));
    end
    fprintf(fileID,'\n');
end
fclose(fileID);
disp(['Wrote ', outFile]);




%% Visualization: pick a set and plot
setIdx    = 10; % choose which set to visualize
pairsToPlot = Sij_cell_unique{setIdx};
geomKey     = Sij_keys{setIdx};

% --- Extract info from key ---
tokens = regexp(geomKey,'dA(\d+)_T(\d+)_(\d+)','tokens','once');
dA_val  = str2double(tokens{1});
t1      = str2double(tokens{2});
t2      = str2double(tokens{3});
n       = transCount;

% Build the 4 symmetry-related transaxial pairs
symPairs = [ t1,      t2;
             t2,      t1;
             (n-1-t1),(n-1-t2);
             (n-1-t2),(n-1-t1) ];

symPairsStr = sprintf('(%d,%d) ', symPairs');

figure('Color','w'); hold on; grid on;
title(sprintf('S_{ij} set #%d: dA=%d, transaxial sym group: %s\n%d LORs', ...
    setIdx, dA_val, symPairsStr, size(pairsToPlot,1)));

% Plot all crystals (light gray)
plot3(pos(:,1), pos(:,2), pos(:,3), '.', 'Color', [0.8 0.8 0.8]);

% Draw LORs for current set (blue lines)
for r = 1:size(pairsToPlot,1)
    i = pairsToPlot(r,1)+1; j = pairsToPlot(r,2)+1;
    line([pos(i,1), pos(j,1)], [pos(i,2), pos(j,2)], [pos(i,3), pos(j,3)], ...
         'LineWidth', 1.5, 'Color', [0 0.5 1]);
end

% Highlight & label crystals in this set
shownIDs = unique(pairsToPlot(:)) + 1;
plot3(pos(shownIDs,1), pos(shownIDs,2), pos(shownIDs,3), ...
    'o','MarkerSize',8,'LineWidth',1.5,'Color',[0 0.3 1]);
for idx = reshape(shownIDs,1,[])
    text(pos(idx,1), pos(idx,2), pos(idx,3), sprintf(' %d', idx-1), ...
        'FontSize',7,'Color',[0 0.28 0.8]);
end

% --- Add module labels ---
for m = 1:modules
    phi = deg2rad(moduleAnglesDeg(m));
    cen = [R_ring*cos(phi), R_ring*sin(phi), 0];
    text(cen(1), cen(2), cen(3)+20, sprintf('Module %d', m), ...
        'FontSize', 10, 'FontWeight', 'bold', 'Color', [0.2 0.1 0.1], ...
        'HorizontalAlignment','center');
end

xlabel('x [mm]'); ylabel('y [mm]'); zlabel('z [mm]');
axis equal; view(35,20); grid on;
legend({'All crystals','LORs','S_{ij} crystals'},'Location','best');
hold off;