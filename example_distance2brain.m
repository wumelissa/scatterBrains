
load ./Subject03/Subject03_volume.mat
load ./Subject03/Subject03_locs.mat

%% distance from surface point to brain tissue

% Units are in mm

% the tissue indices correspond to:
% 1: scalp
% 2: skull
% 3: cerebrospinal fluid
% 4: grey matter
% 5: white matter

point=locs(125,:);

dist=find_distance_of_nearest_voxel(vol,point,4);

% Same result as what is stored in the distance2brain variable
dist_check=distance2brain(125);