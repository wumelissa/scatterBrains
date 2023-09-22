%% add paths

% make sure to install iso2mesh version 1.9.0 (century egg) and mmc 
% version 1.9 (v2020, Moon Cake - beta) in the same folder.

% author: Melissa Wu (wu.melissa.m <at> gmail.com)

% All units are in mm

% EDIT to scatterBrains path directory
scatterBrains_path='.';
cd(scatterBrains_path)

% EDIT to directory where iso2mesh and mmc are installed
iso2mesh_path=[scatterBrains_path filesep 'iso2mesh'];
mmc_path=[scatterBrains_path filesep 'mmc'];

addpath(genpath(iso2mesh_path))
addpath(genpath(mmc_path))

%% load mesh and locations

subject_num='03';

load([scatterBrains_path filesep 'Subject' subject_num filesep 'Subject' subject_num '_mesh.mat'])
load([scatterBrains_path filesep 'Subject' subject_num filesep 'Subject' subject_num '_locs.mat'])

%% visualize tissue layers

% Zoom in and out of the mesh to see the tissue structure inside!

% the tissue indices correspond to:
% 1: scalp
% 2: skull
% 3: cerebrospinal fluid
% 4: grey matter
% 5: white matter

colors=[235 204 171;
    192 192 192;
    153 204 255;
    255 153 204;
    229 76 76];

figure
hold on
for I=1:5
    trimesh(surfaceElem{I}(:,1:3),surfaceNodes{I}(:,1),surfaceNodes{I}(:,2),surfaceNodes{I}(:,3),...
        'FaceAlpha',0.5,'EdgeColor',colors(I,:)/256,'FaceColor','none')
end
axis equal
view(-37.5,30)
zoom(1.5)
title('Zooming in and out lets you see the tissue layers')

%% select a point

[user_point,user_pt_idx]=get_surface_loc_multitissue(surfaceElem,surfaceNodes,locs);
