function [node_closest_to_point,point_unit_vec]=getNormalVec(node,face,specified_point)

% obtains normal vector to center of volume from specified point on head surface

% input:
%   node: array containing node coordinates of mesh, dimension (nnodes,3)
%   face: array containing face coordinates of mesh, dimension (nnodes,4)
%   specified_point: array containing coordinates of specified point, dimension (1,3)

% output:
%   node_closest_to_point: coordinates of node closest to point, dimension (1,3)
%   point_unit_vec: unit normal vector from node closest to point

% author: Melissa Wu (wu.melissa.m <at> gmail.com)
% this file is part of scatterBrains
% License: GPLv3
%%

nodeindices=unique(face(:));
facenodes=node(nodeindices,:);

% finding node closest to point
all_distances_from_isosurf=[];
for idx=1:size(facenodes,1)
    all_distances_from_isosurf(idx)=norm(facenodes(idx,:)-specified_point);
end

[~,sorted_indices]=sort(all_distances_from_isosurf);
node_closest_to_point=facenodes(sorted_indices(1),:);
num_nodes_to_sample=10;

% calculating surface norms for all faces in mesh
snorm=surfacenorm(node,face,'Normalize',1);

% finding all faces for sampled nodes
idx=1;
all_face_arr=[];
for single_node=1:num_nodes_to_sample
    I=sorted_indices(single_node);
    for face_column=1:3
        [row,~]=find(ismember(node,facenodes(I,:)));
        faces_with_source=find(face(:,face_column)==row(1));
        all_face_arr=cat(1,all_face_arr,faces_with_source);
    end
end

% finding normal vectors for all unique faces
unique_faces=unique(all_face_arr);
for f_index=1:length(unique_faces)
    all_point_snorm(idx,:)=snorm(unique_faces(f_index),:);
    idx=idx+1;
end

% making sure normal vectors point inward
midpoint=mean(node);
comparison_vec=midpoint-specified_point;

for I=1:size(all_point_snorm,1)
    dotprods(I)=dot(comparison_vec,all_point_snorm(I,:));
end

point_unit_vec=mean(all_point_snorm(dotprods>0,:));
point_unit_vec=point_unit_vec/sqrt(sum(point_unit_vec.^2));
