function varargout=get_surface_loc(face,node,locs)

% finds surface location nearest to user-specificied point
% inputs:
% face: nface x 3 array of face indices
% node: nnodes x 3 array of node coordinates
% locs: 805 x 3 array of predefined locations

% outputs:
% point (first output): 3d coordinates of surface location
% pt_idx (second output): index of surface location in locs variable

% author: Melissa Wu, wu.melissa.m <at> gmail.com

% this file is part of scatterBrains

arguments
    face (:,3) double
    node (:,3) double
    locs (:,3) double
end
%%

% plot mesh

fig=figure('units','normalized','outerposition',[0 0 1 1]);
trimesh(face(:,1:3),node(:,1),node(:,2),...
    node(:,3),'facecolor','none','FaceAlpha',0.5,'EdgeColor',[192 192 192]/256)
hold on; axis equal
plot3(locs(:,1),locs(:,2),locs(:,3),'k.','MarkerSize',20)
title({'Click on a surface location and then press return',...
    'The program will find the predefined location closest to your selection'})

% get cursor info
datacursormode on
dcm_obj = datacursormode(fig);

pause 

cursor_info = getCursorInfo(dcm_obj);

user_point=cursor_info.Position;

% find closest surface location
for I=1:size(locs,1)
    pt_distances(I)=norm(user_point-locs(I,:));
end

[~,pt_idx]=min(abs(pt_distances));
point=locs(pt_idx,:);

varargout{1}=point;
varargout{2}=pt_idx;

fprintf(['Location is at %.2f %.2f %.2f\n'],point)