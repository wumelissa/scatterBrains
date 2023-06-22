%% add paths

% make sure to install iso2mesh version 1.9.0 (century egg) and mmc 
% version 1.9 (v2020, Moon Cake - beta) in the same folder.

iso2mesh_path=['.' filesep 'iso2mesh'];
mmc_path=['.' filesep 'mmc'];

addpath(genpath(iso2mesh_path))
addpath(genpath(mmc_path))

%% load mesh

subject_num='03';

load(['.' filesep 'Subject' subject_num filesep 'Subject' subject_num '_mesh.mat'])

%% visualize tissue layers

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
