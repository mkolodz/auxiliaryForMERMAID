clear; clc;

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
for a = 0:axCount-1
    for tL = 0:transCount-1
        for tR = 0:transCount-1
            % Generate the 8 equivalent LORs for this combination
            combos = [tL, tR; mirrorT(tL), mirrorT(tR); tR, tL; mirrorT(tR), mirrorT(tL)];
            LORs = [];
            % Pair 1-3 and 2-4
            for k = 1:4
                LORs(end+1,:) = [toID(1,a,combos(k,1)), toID(3,a,combos(k,2))];
                LORs(end+1,:) = [toID(2,a,combos(k,1)), toID(4,a,combos(k,2))];
            end
            % Canonical key (for this axial/transaxial combo)
            key = sprintf('A%d_TL%d_TR%d', a, tL, tR);
            Sij_cell{end+1} = LORs; %#ok<AGROW>
            Sij_keys{end+1} = key; %#ok<AGROW>
        end
    end
end

% --- Output each equivalence set as one line in TXT ---
fileID = fopen('Sij_sets.txt','w');
for s = 1:numel(Sij_cell)
    fprintf(fileID, '%s\t', Sij_keys{s});
    pairs = Sij_cell{s};
    for p = 1:size(pairs,1)
        fprintf(fileID, '%d %d\t', pairs(p,1), pairs(p,2));
    end
    fprintf(fileID, '\n');
end
fclose(fileID);
disp(['Wrote Sij_sets.txt with ' num2str(numel(Sij_cell)) ' equivalence sets.']);

% --- Visualization: plot setIdx-th S_ij set (change setIdx as desired) ---
setIdx = 3; % Example: plot the 3rd set, which now is guaranteed 8 LORs
figure('Color','w'); hold on; grid on;
title(sprintf('Crystals and S_{ij} set #%d', setIdx));

% Plot all crystals
plot3(pos(:,1), pos(:,2), pos(:,3), '.', 'MarkerSize', 8);

% Draw the LORs for current S_ij set
pairsToPlot = Sij_cell{setIdx};
for r = 1:size(pairsToPlot,1)
    i = pairsToPlot(r,1)+1; j = pairsToPlot(r,2)+1;
    line([pos(i,1), pos(j,1)], [pos(i,2), pos(j,2)], [pos(i,3), pos(j,3)], ...
         'LineWidth', 2, 'LineStyle','-');
end

% Label only crystals involved in these LORs
shownIDs = unique(pairsToPlot(:)) + 1;
for idx = reshape(shownIDs,1,[])
    text(pos(idx,1), pos(idx,2), pos(idx,3), sprintf(' %d', idx-1), ...
         'FontSize', 7, 'Color', [0 0 0]);
end

xlabel('x [mm]'); ylabel('y [mm]'); zlabel('z [mm]');
axis equal; view(35,20);
legend({'Crystals','S_{ij} set'}, 'Location','best');
hold off;


% clear; clc;
% % --- Geometry Parameters ---
% modules = 4; axCount = 8; transCount = 16;
% crystalsPerModule = axCount * transCount;
% N = modules * crystalsPerModule;
% % Opposing pairs: [1,3], [2,4]
% pairs_mod = [1,3; 2,4]; % Each row: [left right]
% mirrorT = @(t) transCount-1-t;
% 
% % Mapping from ID to [module, axial, trans]
% toLocal = @(id) deal(floor(id / crystalsPerModule)+1, ...
%                      floor(mod(id, crystalsPerModule)/transCount), ...
%                      mod(id, transCount));
% 
% toID = @(m,a,t) (m-1)*crystalsPerModule + a*transCount + t;
% 
% % --- For each LOR between opposing modules, create its full equivalent set ---
% Sij_cell = {}; set_keys = {};
% for kpair=1:size(pairs_mod,1)
%     leftM = pairs_mod(kpair,1); rightM = pairs_mod(kpair,2);
%     for a = 0:axCount-1
%         for tL = 0:transCount-1
%             for tR = 0:transCount-1
%                 leftID = toID(leftM,a,tL);
%                 rightID = toID(rightM,a,tR);
% 
%                 % Build all four permutations for this LOR
%                 IDs = [
%                     leftID, rightID; % original
%                     %transCount-leftID, transCount-rightID;
%                     toID(leftM,a,mirrorT(tL)), toID(rightM,a,mirrorT(tR)); % mirror both
%                     rightID, leftID; % swap
%                     %transCount-rightID, transCount-leftID;
%                     toID(rightM,a,mirrorT(tR)), toID(leftM,a,mirrorT(tL)); % swap and mirror both
%                 ];
% 
%                 IDs = unique(IDs, 'rows'); % Remove duplicates
% 
%                 % Build analogous set in opposing module pair
%                 opp_leftM = mod(leftM,modules)+1; %1->2, 2->3,3->4,4->1, so just swap to the other pair
%                 opp_rightM = mod(rightM,modules)+1;
%                 opp_IDs = [
%                     toID(opp_leftM,a,tL), toID(opp_rightM,a,tR);
%                     toID(opp_leftM,a,mirrorT(tL)), toID(opp_rightM,a,mirrorT(tR));
%                     toID(opp_rightM,a,tR), toID(opp_leftM,a,tL);
%                     toID(opp_rightM,a,mirrorT(tR)), toID(opp_leftM,a,mirrorT(tL));
%                 ];
%                 opp_IDs = unique(opp_IDs, 'rows');
%                 % Full set
%                 setLORs = [IDs; opp_IDs];
% 
%                 % Use a set key based on (module pair, axial ring, transaxial positions)
%                 set_key = sprintf('L%dR%dA%dTL%dTR%d',leftM,rightM,a,tL,tR);
% 
%                 % Add if new
%                 if all(~strcmp(set_keys, set_key))
%                     Sij_cell{end+1} = setLORs; %#ok<AGROW>
%                     set_keys{end+1} = set_key; %#ok<AGROW>
%                 end
%             end
%         end
%     end
% end
% 
% % --- Write output files ---
% f1 = fopen('Sij_sets.txt','w');
% for s = 1:numel(Sij_cell)
%     pairs = Sij_cell{s};
%     fprintf(f1,'Set_%d\t',s);
%     for p = 1:size(pairs,1)
%         fprintf(f1,'%d %d\t',pairs(p,1),pairs(p,2));
%     end
%     fprintf(f1,'\n');
% end
% fclose(f1);
% 
% fprintf('Wrote Sij_sets.txt with %d unique sets\n', numel(Sij_cell));
% 
% % --- Compute crystal positions ---
% R_ring = 200; % mm
% moduleAnglesDeg = [-60, 60, 120, 240];
% pitchT = 5.0; pitchA = 5.0;
% halfT = (transCount-1)/2;
% halfA = (axCount-1)/2;
% pos = zeros(N,3);
% modIdx = zeros(N,1); axIdx = zeros(N,1); trIdx = zeros(N,1);
% 
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
% % --- Select S_ij sets to plot ---
% setIdx = 3; % Change this to plot different sets (e.g., first or chosen set)
% pairsToPlot = Sij_cell{setIdx};
% 
% figure('Color','w'); hold on; grid on;
% title(sprintf('Crystals and S_{ij} set #%d', setIdx));
% plot3(pos(:,1), pos(:,2), pos(:,3), '.', 'MarkerSize', 8);
% 
% % Draw the LORs for current S_ij set
% for r = 1:size(pairsToPlot,1)
%     i = pairsToPlot(r,1)+1; j = pairsToPlot(r,2)+1;
%     line([pos(i,1), pos(j,1)], [pos(i,2), pos(j,2)], [pos(i,3), pos(j,3)], ...
%          'LineWidth', 1.5, 'LineStyle','-');
% end
% 
% % Label only crystals involved in these LORs
% shownIDs = unique([pairsToPlot(:)]) + 1;
% for idx = reshape(shownIDs,1,[])
%     text(pos(idx,1), pos(idx,2), pos(idx,3), sprintf(' %d', idx-1), ...
%          'FontSize', 7, 'Color', [0 0 0]);
% end
% 
% xlabel('x [mm]'); ylabel('y [mm]'); zlabel('z [mm]');
% axis equal; view(35,20);
% legend({'Crystals','S_{ij} set'}, 'Location','best');
% hold off;
% % ==== S_ij (geometrically equivalent LORs) for a 4-module, 512-crystal PET ====
% % Geometry: 4 modules, each 8 (axial) x 16 (transaxial) crystals.
% % Mapping (original, row-wise):
% %   Module 1 IDs: 0..127 (rows of 16: 0..15, 16..31, ..., 112..127)
% %   Module 2 IDs: 128..255,  Module 3: 256..383,  Module 4: 384..511
% %
% % S_ij definition (Kinouchi et al., Fig. 1(b)): All LORs "geometrically equal"
% % to LOR-ij -> same relative positions within modules AND same inter-detector angle.
% % We implement this by grouping LORs with identical (Δt, Δa) after mirroring
% % the transaxial index of the opposite module (so left/right module frames align).
% %
% % Outputs:
% %   - LORs: all valid LOR pairs [i j] (i<j) between opposing modules
% %   - S: cell array; each cell is an Nx2 list of crystal ID pairs for one S_ij
% %   - Files:
% %       * Sij_sets.txt       (one line per S_ij set; pairs "i j" separated by tabs)
% %       * Sij_key_map.txt    (key index -> [dT dA] and set size)
% %
% % Visualization:
% %   - 3D plot of all crystals (with crystal ID labels for the shown sets to avoid clutter)
% %   - Lines drawn for two S_ij demo sets:
% %       * same axial ring (Δa = 0), solid lines
% %       * different rings (Δa = +1), dotted lines
% %
% % ------------------------------------------------------------------------------
% 
% clear; clc;
% 
% %% --------- Basic scanner/grid parameters ----------
% modules            = 4;
% axCount            = 8;   % axial crystals per module (rows)
% transCount         = 16;  % transaxial crystals per module (cols)
% crystalsPerModule  = axCount * transCount;           % 128
% N                  = modules * crystalsPerModule;    % 512 total
% 
% % Module angular placement (degrees):
% % Two modules on each "side", separated by 120°, and the two sides are 180° apart.
% % Opposing pairs: (1 <-> 3) and (2 <-> 4).
% moduleAnglesDeg = [-60, +60, 120, 240];
% 
% % Simple physical layout (for visualization only)
% R_ring       = 200;    % mm: radius to module center (arbitrary, for plotting)
% pitchT       = 5.0;    % mm: transaxial crystal pitch
% pitchA       = 5.0;    % mm: axial crystal pitch
% halfT        = (transCount-1)/2;   % 7.5
% halfA        = (axCount-1)/2;      % 3.5
% 
% % Opposing pairs (by module index)
% isOpposing = @(m1,m2) ((m1==1 && m2==3) || (m1==3 && m2==1) || ...
%                        (m1==2 && m2==4) || (m1==4 && m2==2));
% 
% %% --------- ID <-> (module, axial, trans) mapping ----------
% % ID in [0..511]; module m in [1..4]; a in [0..7], t in [0..15]
% toLocal = @(id) deal( ...
%     floor(id / crystalsPerModule) + 1, ...
%     floor(mod(id, crystalsPerModule) / transCount), ...
%     mod(id, transCount) );
% 
% toID = @(m,a,t) (m-1)*crystalsPerModule + a*transCount + t;
% 
% %% --------- Crystal 3D positions ----------
% % Each module has a center on the ring; local tangent = rotate radial by +90°
% pos = zeros(N,3);
% modIdx = zeros(N,1); axIdx = zeros(N,1); trIdx = zeros(N,1);
% 
% for m = 1:modules
%     phi  = deg2rad(moduleAnglesDeg(m));
%     cen  = [R_ring*cos(phi), R_ring*sin(phi), 0];
%     tHat = [-sin(phi), cos(phi), 0];   % unit tangent at module center
%     % z axis is axial
%     for a = 0:axCount-1
%         z = (a - halfA) * pitchA;
%         for t = 0:transCount-1
%             xyt = cen + tHat * ((t - halfT)*pitchT);
%             id  = toID(m,a,t);
%             pos(id+1,:) = [xyt(1), xyt(2), z];
%             modIdx(id+1) = m;
%             axIdx(id+1)  = a;
%             trIdx(id+1)  = t;
%         end
%     end
% end
% 
% %% --------- Build all valid LORs between opposing modules ----------
% pairs = [];  % [i j] with i<j
% for i = 0:N-1
%     mi = modIdx(i+1);
%     for j = i+1:N-1
%         mj = modIdx(j+1);
%         if isOpposing(mi,mj)
%             pairs(end+1,:) = [i, j]; %#ok<AGROW>
%         end
%     end
% end
% LORs = pairs;
% disp(size(LORs));
% %% --------- Geometric-equivalence key per LOR (Kinouchi S_ij) ----------
% % For opposing modules, align local frames by mirroring the "right" module's trans index.
% % Canonicalize left/right so that the "left" module is 1 or 2 (within its opposing pair).
% mirrorT = @(t) (transCount-1 - t);
% 
% keys = zeros(size(LORs,1), 2);   % columns: [dT, dA]
% for r = 1:size(LORs,1)
%     i = LORs(r,1); j = LORs(r,2);
%     [mi, ai, ti] = toLocal(i);
%     [mj, aj, tj] = toLocal(j);
% 
%     % Force order: (left,right) = (1,3) or (2,4). If we have (3,1) or (4,2), swap.
%     leftID = i; rightID = j;
%     if (mi==3 && mj==1) || (mi==4 && mj==2)
%         leftID = j; rightID = i;
%         [mi, ai, ti, mj, aj, tj] = deal(mj, aj, tj, mi, ai, ti);
%     end
% 
%     % After 180° rotation, the right module's trans index must be mirrored to align frames.
%     tj_m = mirrorT(tj);
% 
%     % Geometric key: same inter-detector angle & relative positions -> same (Δt, Δa)
%     dT = tj_m - ti;        % transaxial offset across the pair (aligned frames)
%     dA = aj   - ai;        % axial offset across the pair
%     keys(r,:) = [dT, dA];
% end
% 
% % Unique S_ij classes:
% [uniqKeys, ~, keyIdx] = unique(keys, 'rows', 'stable');
% numSets = size(uniqKeys,1);
% 
% S = cell(numSets,1);
% for k = 1:numSets
%     S{k} = LORs(keyIdx==k, :);  % Nx2 pairs for this S_ij
% end
% 
% fprintf('Total LORs (opposing modules only): %d\n', size(LORs,1));
% fprintf('Unique S_ij classes (Kinouchi geometric equivalence): %d\n', numSets);
% 
% %% --------- Save S_ij sets: one line per set ----------
% % Each line: "k  dT dA  |  i0 j0    i1 j1    i2 j2  ..."
% s1 = fopen('Sij_sets.txt','w');
% s2 = fopen('Sij_key_map.txt','w');
% fprintf(s2, '# k\t dT\t dA\t set_size\n');
% for k = 1:numSets
%     dT = uniqKeys(k,1); dA = uniqKeys(k,2);
%     P  = S{k};
%     fprintf(s2, '%d\t %d\t %d\t %d\n', k, dT, dA, size(P,1));
% 
%     fprintf(s1, '%d\t%d\t%d\t|\t', k, dT, dA);
%     % write all pairs on one line, each "i j"
%     for r = 1:size(P,1)
%         fprintf(s1, '%d %d\t', P(r,1), P(r,2));
%     end
%     fprintf(s1, '\n');
% end
% fclose(s1); fclose(s2);
% fprintf('Wrote Sij_sets.txt and Sij_key_map.txt\n');
% 
% %% --------- Pick two demo S_ij sets for plotting (like Fig. 1(b)) ----------
% % 1) same axial ring: Δa = 0, and choose Δt = 0 if it exists; else closest
% pickSame = find(uniqKeys(:,2)==0, 1, 'first'); % first with dA=0
% if isempty(pickSame), pickSame = 1; end
% 
% % 2) different rings: Δa = +1 (or any nonzero); pick first available
% pickDiff = find(uniqKeys(:,2)~=0, 1, 'first');
% if isempty(pickDiff), pickDiff = min(numSets, pickSame+1); end
% 
% S_same = S{pickSame};
% S_diff = S{pickDiff};
% 
% %% --------- 3D plot of crystals + LORs (lines) ----------
% figure('Color','w'); hold on; grid on;
% title(sprintf('Crystals and S_{ij} sets (solid: dA=%d, dotted: dA=%d)', ...
%     uniqKeys(pickSame,2), uniqKeys(pickDiff,2)));
% 
% % Plot all crystals (light markers)
% plot3(pos(:,1), pos(:,2), pos(:,3), '.', 'MarkerSize', 8);
% 
% % Draw LORs for the "same ring" set (solid)
% for r = 1:size(S_same,1)
%     i = S_same(r,1)+1; j = S_same(r,2)+1;
%     line([pos(i,1), pos(j,1)], [pos(i,2), pos(j,2)], [pos(i,3), pos(j,3)], ...
%         'LineWidth', 1.5, 'LineStyle','-');
% end
% 
% % Draw LORs for the "different rings" set (dotted)
% for r = 1:size(S_diff,1)
%     i = S_diff(r,1)+1; j = S_diff(r,2)+1;
%     line([pos(i,1), pos(j,1)], [pos(i,2), pos(j,2)], [pos(i,3), pos(j,3)], ...
%         'LineWidth', 1.5, 'LineStyle',':');
% end
% 
% % Label only crystals that participate in the shown sets (to avoid clutter)
% shownIDs = unique([S_same(:); S_diff(:)]) + 1;  % 1-based
% for idx = reshape(shownIDs,1,[])
%     text(pos(idx,1), pos(idx,2), pos(idx,3), sprintf(' %d', idx-1), ...
%         'FontSize', 7, 'Color', [0 0 0]);
% end
% 
% xlabel('x [mm]'); ylabel('y [mm]'); zlabel('z [mm]');
% axis equal; view(35,20);
% legend({'Crystals','S_{ij} (same ring)','S_{ij} (different rings)'}, 'Location','best');
% hold off;











%% --------- Also provide the flat array of unique LORs per set (as requested) ----------
% S is already a cell array of 2-column arrays of ID pairs. Example usage:
%   S{1}  -> Nx2 [id_i id_j] pairs for the first equivalence set.

% Done.


% function build_equivalent_LORs_originalMap_3D_withFilter()
% 
% % Geometry
% numModules = 4;
% axialCount = 8;    
% transCount = 16;   
% crystalsPerModule = axialCount * transCount; % 128
% totalCrystals = numModules * crystalsPerModule;
% 
% % Opposing module pairs
% facingPairs = [1 3; 2 4];
% 
% %% ---- Original crystal ID mapping (row-major) ----
% crystalID = @(m,a,t) (m-1)*crystalsPerModule + a*transCount + t;
% 
% % Precompute indices
% axIdx = zeros(totalCrystals,1);
% trIdx = zeros(totalCrystals,1);
% modIdx = zeros(totalCrystals,1);
% for m = 1:numModules
%     for a = 0:axialCount-1
%         for t = 0:transCount-1
%             gid = crystalID(m,a,t);
%             axIdx(gid+1) = a;
%             trIdx(gid+1) = t;
%             modIdx(gid+1) = m;
%         end
%     end
% end
% 
% %% ---- Build LORs between opposing modules ----
% LORs = [];
% for k = 1:size(facingPairs,1)
%     mA = facingPairs(k,1);
%     mB = facingPairs(k,2);
%     idsA = find(modIdx == mA) - 1; 
%     idsB = find(modIdx == mB) - 1;
%     [IA, IB] = ndgrid(idsA, idsB);
%     pairs = [IA(:), IB(:)];
%     pairs = pairs(pairs(:,1) < pairs(:,2), :);
%     LORs = [LORs; pairs]; %#ok<AGROW>
% end
% 
% %% ---- Group into equivalence classes (Kinouchi definition) ----
% keys = zeros(size(LORs,1), 6);
% for r = 1:size(LORs,1)
%     i = LORs(r,1)+1;
%     j = LORs(r,2)+1;
%     ai = axIdx(i); ti = trIdx(i); mi = modIdx(i);
%     aj = axIdx(j); tj = trIdx(j); mj = modIdx(j);
% 
%     % enforce consistent module ordering
%     if mi < mj
%         keys(r,:) = [mi mj ti ai tj aj];
%     else
%         keys(r,:) = [mj mi tj aj ti ai];
%     end
% end
% 
% [uniqueKeys,~,idxGroup] = unique(keys,'rows');
% S = cell(size(uniqueKeys,1),1);
% for g = 1:max(idxGroup)
%     S{g} = LORs(idxGroup==g,:);
% end
% 
% fprintf('Total unique equivalence classes (S_ij): %d\n', numel(S));
% 
% % %% ---- Group into equivalence classes (Δtrans, Δaxial) ----
% % keys = [trIdx(LORs(:,1)+1) - trIdx(LORs(:,2)+1), ...
% %         axIdx(LORs(:,1)+1) - axIdx(LORs(:,2)+1)];
% % for r = 1:size(keys,1)
% %     if keys(r,1) < 0 || (keys(r,1)==0 && keys(r,2) < 0)
% %         keys(r,:) = -keys(r,:);
% %     end
% % end
% % [uniqueKeys,~,idxGroup] = unique(keys,'rows');
% % S = cell(size(uniqueKeys,1),1);
% % for g = 1:max(idxGroup)
% %     S{g} = LORs(idxGroup==g,:);
% % end
% % fprintf('Total unique equivalence classes (S_ij): %d\n', numel(S));
% 
% %% ---- Build 3D positions of crystals ----
% R = 100;                 % detector ring radius (mm)
% dz = 5;                  % axial spacing (mm)
% pitchTx = 5;             % transaxial pitch (mm)
% modAngles = [-60, +60, +120, +240]; % deg
% 
% crystalPos = NaN(totalCrystals,3);
% for m = 1:numModules
%     ang = deg2rad(modAngles(m));
%     u_hat = [cos(ang), sin(ang), 0];    % radial
%     v_hat = [-sin(ang), cos(ang), 0];   % tangential
%     z_hat = [0 0 1];                    % axial
%     modCenter = R * u_hat;
%     for a = 0:axialCount-1
%         for t = 0:transCount-1
%             gid = crystalID(m,a,t);
%             offsetT = (t - (transCount-1)/2)*pitchTx;
%             offsetZ = (a - (axialCount-1)/2)*dz;
%             pos = modCenter + offsetT*v_hat + offsetZ*z_hat;
%             crystalPos(gid+1,:) = pos;
%         end
%     end
% end
% 
% %% ---- Plot one equivalence set in 3D ----
% exampleGroup = 1;   % pick a group index
% exampleLORs = S{exampleGroup};
% chosenRow   = 0;    % set [] = all rows, or pick 0..7
% 
% figure; hold on; grid on; axis equal;
% title(sprintf('3D Equivalence set S_{ij} #%d (Δt=%d, Δa=%d)', ...
%     exampleGroup, uniqueKeys(exampleGroup,1), uniqueKeys(exampleGroup,2)));
% xlabel('x (mm)'); ylabel('y (mm)'); zlabel('z (mm)');
% view(3);
% 
% % Draw crystals + IDs
% for gid = 1:totalCrystals
%     row = axIdx(gid);
%     if isempty(chosenRow) || row == chosenRow
%         plot3(crystalPos(gid,1),crystalPos(gid,2),crystalPos(gid,3),...
%             'bo','MarkerFaceColor','b');
%         text(crystalPos(gid,1),crystalPos(gid,2),crystalPos(gid,3)+2,...
%             num2str(gid-1),'FontSize',6,'Color','k');
%     end
% end
% 
% % Draw LORs
% for k = 1:size(exampleLORs,1)
%     i = exampleLORs(k,1)+1;
%     j = exampleLORs(k,2)+1;
%     row_i = axIdx(i);
%     row_j = axIdx(j);
%     if isempty(chosenRow) || (row_i==chosenRow && row_j==chosenRow)
%         plot3([crystalPos(i,1), crystalPos(j,1)], ...
%               [crystalPos(i,2), crystalPos(j,2)], ...
%               [crystalPos(i,3), crystalPos(j,3)], 'r-','LineWidth',1);
%     end
% end
% 
% end
% 


% function build_equivalent_LORs_originalMap()
% 
% % Geometry
% numModules = 4;
% axialCount = 8;    
% transCount = 16;   
% crystalsPerModule = axialCount * transCount; % 128
% totalCrystals = numModules * crystalsPerModule;
% 
% % Opposing module pairs
% facingPairs = [1 3; 2 4];
% 
% %% ---- Original crystal ID mapping (row-major) ----
% % Module m (1..4), axial a (0..7), trans t (0..15)
% crystalID = @(m,a,t) (m-1)*crystalsPerModule + a*transCount + t;
% 
% % Precompute indices for each crystal
% axIdx = zeros(totalCrystals,1);
% trIdx = zeros(totalCrystals,1);
% modIdx = zeros(totalCrystals,1);
% for m = 1:numModules
%     for a = 0:axialCount-1
%         for t = 0:transCount-1
%             gid = crystalID(m,a,t);
%             axIdx(gid+1) = a;
%             trIdx(gid+1) = t;
%             modIdx(gid+1) = m;
%         end
%     end
% end
% 
% %% ---- Build all LORs between opposing modules ----
% LORs = [];
% for k = 1:size(facingPairs,1)
%     mA = facingPairs(k,1);
%     mB = facingPairs(k,2);
%     idsA = find(modIdx == mA) - 1; % zero-based IDs
%     idsB = find(modIdx == mB) - 1;
%     [IA, IB] = ndgrid(idsA, idsB);
%     pairs = [IA(:), IB(:)];
%     pairs = pairs(pairs(:,1) < pairs(:,2), :);
%     LORs = [LORs; pairs]; %#ok<AGROW>
% end
% 
% %% ---- Group LORs into equivalence classes (Δtrans, Δaxial) ----
% keys = [trIdx(LORs(:,1)+1) - trIdx(LORs(:,2)+1), ...
%         axIdx(LORs(:,1)+1) - axIdx(LORs(:,2)+1)];
% for r = 1:size(keys,1)
%     if keys(r,1) < 0 || (keys(r,1)==0 && keys(r,2) < 0)
%         keys(r,:) = -keys(r,:);
%     end
% end
% [uniqueKeys,~,idxGroup] = unique(keys,'rows');
% S = cell(size(uniqueKeys,1),1);
% for g = 1:max(idxGroup)
%     S{g} = LORs(idxGroup==g,:);
% end
% 
% fprintf('Total unique equivalence classes (S_ij): %d\n', numel(S));
% 
% %% ---- Plot one equivalence set (like Fig. 1b) ----
% exampleGroup = 10; % pick a group index
% exampleLORs = S{exampleGroup};
% 
% R = 100; % ring radius
% modAngles = [-60, +60, +120, +240]; % deg
% chosenRow = 4; % axial row to display
% 
% crystalPos = NaN(totalCrystals,2);
% for m = 1:numModules
%     ang = deg2rad(modAngles(m));
%     u_hat = [cos(ang), sin(ang)];
%     v_hat = [-sin(ang), cos(ang)];
%     modCenter = R * u_hat;
%     for t = 0:transCount-1
%         gid = crystalID(m,chosenRow,t);
%         offset = (t - (transCount-1)/2);
%         pos = modCenter + offset * v_hat*5; % spacing
%         crystalPos(gid+1,:) = pos;
%     end
% end
% 
% figure; hold on; axis equal;
% title(sprintf('Equivalence set S_{ij} #%d (Δt=%d, Δa=%d)', ...
%     exampleGroup, uniqueKeys(exampleGroup,1), uniqueKeys(exampleGroup,2)));
% xlabel('x'); ylabel('y');
% 
% % Draw LORs
% for k = 1:size(exampleLORs,1)
%     i = exampleLORs(k,1)+1;
%     j = exampleLORs(k,2)+1;
%     if all(~isnan(crystalPos([i j],1)))
%         plot(crystalPos([i j],1), crystalPos([i j],2),'r-','LineWidth',1.2);
%     end
% end
% 
% % Draw crystals + labels
% for gid = 1:totalCrystals
%     if ~isnan(crystalPos(gid,1))
%         plot(crystalPos(gid,1),crystalPos(gid,2),'bo','MarkerFaceColor','b');
%         text(crystalPos(gid,1)+2,crystalPos(gid,2),num2str(gid-1),...
%             'FontSize',7,'Color','k');
%     end
% end
% 
% end



% function build_equivalent_LORs()
% 
% % Geometry
% numModules = 4;
% axialCount = 8;    % rows
% transCount = 16;   % cols
% crystalsPerModule = axialCount * transCount;
% 
% % Opposing module pairs
% facingPairs = [1 3; 2 4];
% 
% % Crystal ID mapping: (module, axial, trans) -> global ID
% crystalID = @(m,a,t) (m-1)*crystalsPerModule + a*transCount + t;
% 
% % Precompute axial & trans indices for each crystal ID
% axIdx = zeros(512,1);
% trIdx = zeros(512,1);
% modIdx = zeros(512,1);
% for m = 1:numModules
%     for a = 0:axialCount-1
%         for t = 0:transCount-1
%             gid = crystalID(m,a,t);
%             axIdx(gid+1) = a;
%             trIdx(gid+1) = t;
%             modIdx(gid+1) = m;
%         end
%     end
% end
% 
% % Collect all valid LORs (only opposing modules)
% LORs = [];
% for k = 1:size(facingPairs,1)
%     mA = facingPairs(k,1);
%     mB = facingPairs(k,2);
%     idsA = find(modIdx == mA) - 1; % zero-based IDs
%     idsB = find(modIdx == mB) - 1;
%     [IA, IB] = ndgrid(idsA, idsB);
%     pairs = [IA(:), IB(:)];
%     % enforce i<j
%     pairs = pairs(pairs(:,1) < pairs(:,2), :);
%     LORs = [LORs; pairs]; %#ok<AGROW>
% end
% 
% % --- Group into equivalent sets S_ij
% % Key = [Δtrans Δaxial]
% keys = [trIdx(LORs(:,1)+1) - trIdx(LORs(:,2)+1), ...
%         axIdx(LORs(:,1)+1) - axIdx(LORs(:,2)+1)];
% 
% % Normalize key so that direction doesn't matter
% % (Δt,Δa) and (-Δt,-Δa) are the same equivalence
% for r = 1:size(keys,1)
%     if keys(r,1) < 0 || (keys(r,1)==0 && keys(r,2) < 0)
%         keys(r,:) = -keys(r,:);
%     end
% end
% 
% % Group by unique keys
% [uniqueKeys,~,idxGroup] = unique(keys, 'rows');
% 
% % Build cell array of equivalent LOR sets
% S = cell(size(uniqueKeys,1),1);
% for g = 1:max(idxGroup)
%     S{g} = LORs(idxGroup==g,:);
% end
% 
% % Display summary
% fprintf('Total unique equivalence classes (S_ij): %d\n', numel(S));
% 
% % Example: print first 3 groups
% for g = 1:min(3,numel(S))
%     fprintf('S{%d} (Δt=%d, Δa=%d):\n', g, uniqueKeys(g,1), uniqueKeys(g,2));
%     disp(S{g}(1:max(5,end),:)); % show first few LORs
% end
% 
% % Save to MAT file
% save('Equivalent_LORs.mat','S','uniqueKeys');
% 
% end
% 
% build_equivalent_LORs()