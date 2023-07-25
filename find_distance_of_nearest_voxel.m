function distance=find_distance_of_nearest_voxel(volume,point,tissue_index)

% finds distance from the point to the nearest voxel of the specified
% tissue index

% inputs:
% volume: subject 3d volume
% point: 3d coordinates of point
% tissue_index: tissue index number

% output:
% distance: distance between point and nearest voxel of tissue index 

% author: Melissa Wu (wu.melissa.m <at> gmail.com)

% this file is part of scatterBrains

[r,c,v]=ind2sub(size(volume),find(volume==tissue_index));

for I=1:length(r)
    alldistances(I)=norm(point-[r(I),c(I),v(I)]);
end

[distance,~]=min(alldistances);