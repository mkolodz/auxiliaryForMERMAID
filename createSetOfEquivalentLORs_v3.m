
function [setOfEqivalentLORs, setIDtoLOR_lookup] = createSetOfEquivalentLORs_v3()
clc; clear;
%% Geometry parameters
modules = 4;
axCount = 8;
transCount = 16;
crystalsPerModule = axCount * transCount; % 128
N = modules * crystalsPerModule; % 512
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
% leftIDs = find(modIdx==1 | modIdx==2) - 1; % zero-based
% rightIDs = find(modIdx==3 | modIdx==4) - 1; % zero-based
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
combos = [ tL, tR;
tR, tL;
(n-1-tL),(n-1-tR);
(n-1-tR),(n-1-tL) ];
%disp(combos)
for k = 1:size(combos,1)
tLi = combos(k,1);
tRi = combos(k,2);
% Canonicalize transaxial pair (order-independent)
tmin = min(tLi, tRi);
tmax = max(tLi, tRi);
% Build IDs
idL = toID(mL, aL, tLi);
idR = toID(mR, aR, tRi);
%disp(idL+ " "+ idR)
% Geometry key: |dA| plus canonicalized transaxial pair
% key = sprintf('dA%02d_T%02d_%02d', dA, tLi, tRi);
% Canonical geometry key: dA + transaxial symmetry group (order-independent)
symPairs = [ tLi, tRi;
tRi, tLi;
(n-1-tLi),(n-1-tRi);
(n-1-tRi),(n-1-tLi) ];
% Canonicalize: pick the lexicographically smallest representation
symPairs_sorted = sortrows(symPairs);
t1 = symPairs_sorted(1,1);
t2 = symPairs_sorted(1,2);
key = sprintf('dA%02d_combo%02d_%02d', dA, t1, t2);
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
S_ij_cell_all = Sij_cell;
%% to chyba nie robi nic
for k = 1:numel(Sij_keys)
pairs = S_map(Sij_keys{k});
pairs = unique(pairs,'rows');
pairs = sortrows(pairs,[1 2]);
Sij_cell{k} = pairs;
end
Sij_cell_unique = Sij_cell; % final set of equivalence classes
fprintf('Total unique S_ij sets: %d\n', numel(Sij_cell_unique));
%% Output each set as one line in TXT
outFile = 'C:\Users\Magdalena\Documents\LuebeckWORK\Software\auxiliaryForMERMAID\Sij_sets.txt';
fileID = fopen(outFile,'w');
for s = 1:numel(Sij_cell_unique)
pairs = Sij_cell_unique{s};
key = Sij_keys{s};
fprintf(fileID,'%s\t', Sij_keys{s}); % geometry key
%fprintf(fileID,'%s\t', key); % print geometry key first
for p = 1:size(pairs,1)
fprintf(fileID,'%d %d\t', pairs(p,1), pairs(p,2));
end
fprintf(fileID,'\n');
end
fclose(fileID);
disp(['Wrote ', outFile]);
%% Visualization: pick a set and plot
setIdx = 10; % choose which set to visualize
pairsToPlot = Sij_cell_unique{setIdx};
geomKey = Sij_keys{setIdx};
% --- Extract info from key ---
tokens = regexp(geomKey,'dA(\d+)_combo(\d+)_(\d+)','tokens','once');
dA_val = str2double(tokens{1});
t1 = str2double(tokens{2});
t2 = str2double(tokens{3});
% Build full symmetry group (for labeling in title)
n = transCount;
symPairs = [ t1, t2;
t2, t1;
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
setOfEqivalentLORs = Sij_cell_unique;
% for i=1:256
%     for j=257:512
%         LOR_pairs(i, j-256) = [i,j];
%     end
% end
A = 1:256;
B = 257:512;

% Create all combinations
[X, Y] = ndgrid(A, B);

% Put into 2-column array
LOR_pairs = [X(:), Y(:)];

setIDtoLOR_lookup = LOR_pairs;
