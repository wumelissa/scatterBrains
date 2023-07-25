
load ./Subject03/Subject03_volume.mat
load ./Subject03/Subject03_locs.mat

%% distance from surface point to brain tissue

% Units are in mm

point=locs(125,:);

distance2brain=find_distance_of_nearest_voxel(vol,point,5);