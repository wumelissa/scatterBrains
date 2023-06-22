% make sure to install iso2mesh version 1.9.0 (century egg) and mmc 
% version 1.9 (v2020, Moon Cake - beta) in the same folder.

iso2mesh_path=['.' filesep 'iso2mesh'];
mmc_path=['.' filesep 'mmc'];

addpath(genpath(iso2mesh_path))
addpath(genpath(mmc_path))

%% convert volume to mesh, save mesh

subject_num='03';
load(['.' filesep 'Subject' subject_num filesep 'Subject' subject_num '_volume.mat'])

newvol=uint8(vol);
[node,elem,face]=v2m(newvol,[],5,200,'cgalmesh');

savemmcmesh(['subject' subject_num],node(:,1:3),elem);

%% generate surface meshes for visualization

for tiss_type=1:5;
    [surfaceNodes{tiss_type},surfaceElem{tiss_type}]=v2s(newvol,tiss_type,1.5);
end

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

%%

save(['.' filesep 'Subject' subject_num filesep 'Subject' subject_num 'mesh.mat'],'node','elem','face','surfaceElem','surfaceNodes')