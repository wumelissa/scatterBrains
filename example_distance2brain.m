
load ./Subject03/Subject03_volume.mat
load ./Subject03/Subject03_locs.mat

%%

point=locs(125,:);

distance2brain=find_distance_of_nearest_voxel(vol,point,5);