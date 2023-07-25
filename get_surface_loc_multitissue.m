function varargout=get_surface_loc_multitissue(surfaceElem,surfaceNodes,locs)

% finds surface location nearest to user-specified point
% inputs:
% surfaceElem: cell array of surface elements for each tissue type
% surfaceNodes: cell array of surface nodes for each tissue type
% locs: 805 x 3 array of predefined locations

% outputs:
% point (first output): 3d coordinates of surface location
% pt_idx (second output): index of surface location in locs variable

% author: Melissa Wu, wu.melissa.m <at> gmail.com

% this file is part of scatterBrains
%%

% plot mesh

colors=[235 204 171;
    192 192 192;
    153 204 255;
    255 153 204;
    229 76 76];

fig=figure;
hold on
for I=1:length(surfaceElem)
    trimesh(surfaceElem{I}(:,1:3),surfaceNodes{I}(:,1),surfaceNodes{I}(:,2),surfaceNodes{I}(:,3),...
        'FaceAlpha',0.5,'EdgeColor',colors(I,:)/256,'FaceColor','none')
end
axis equal
plot3(locs(:,1),locs(:,2),locs(:,3),'k.','MarkerSize',20)
view(-37.5,30)
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