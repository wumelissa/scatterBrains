function [source,det_arr,source_unit_vec,det_unit_vec]=wrapProbe(node,face,point,rz,lz,det_distances,varargin)

% returns source and detector coordinates as well as normal vectors for
% source and detector array

% inputs:
% node: nnodes x 3 array of node coordinates
% face: nface x 3 array of face indices
% point: point of interest on surface of head
% rz: right ear point
% lz: left ear point
% det_distances: desired distances in mm of detectors
% varargin: set flag to 1 to show plot

% outputs:
% source: source coordinates (1,3)
% det_arr: detector coordinates (ndet,3)
% source_unit_vec: unit normal vector from source (1,3)
% det_unit_vec: unit normal vector from each detector (1,3)

% part of this code was taken from iso2mesh demo scripts:
% demo_qmeshcut_ex1.m, demo_surf2vol_ex1.m
% of which the author is Qianqian Fang <q.fang at neu.edu>

% author: Melissa Wu, wu.melissa.m <at> gmail.com

% this file is part of scatterBrains

% if nargin<7, showfig=0;end
if nargin<6, det_distances=5:5:40;end

% if isempty(varargin)
%     showfig=0;
% else
%     showfig=varargin{1};
% end

%% find outline of intersection of plane + volume

% taken from iso2mesh demo code - start
[cut_pos,~,cut_edges]=qmeshcut(face,node,node(:,1),[rz;lz;point]);
[cut_pos,cut_edges]=removedupnodes(cut_pos,cut_edges);
cut_loop=extractloops(cut_edges);

cut_loop(isnan(cut_loop))=[];
% taken from iso2mesh demo code - end  

outer_line=cut_pos(cut_loop,:);

% remove any points that are on the ears
% right_ear_indices=outer_line(:,1)<(rz(1)-1) & outer_line(:,2) > (rz(2)-10) & outer_line(:,3) < rz(3);
% left_ear_indices=outer_line(:,1)>(lz(1)-1) & outer_line(:,2) > (lz(2)-10) & outer_line(:,3) < lz(3);
% all_indices=boolean(right_ear_indices+left_ear_indices);
% outer_line(all_indices,:)=[];
% outer_line=outer_line(outer_line(:,1)>(rz(1)-1) & outer_line(:,2) < (rz(2)-10),:);
% outer_line=outer_line(outer_line(:,1)<(lz(1)+1) & outer_line(:,2) > (lz(2)-10),:);

outer_line=unique(outer_line,'rows','stable');

%% find closest outline point to POI

distance_to_start=norm(point-outer_line(1,:));

for I=1:size(outer_line,1)
    distance_to_points(I)=norm(point-outer_line(I,:));
end

[~,closestLinePt]=min(distance_to_points);

%% adjust "starting point" of the outline so it's far from POI

shift_tolerance=max(det_distances);
shift_length=floor(size(outer_line,1)/2);

% shift if too close to fiducial
if distance_to_start<shift_tolerance || closestLinePt <shift_tolerance || (size(outer_line,1) - closestLinePt)< shift_tolerance
    outer_line_shift=outer_line(shift_length:end,:);
    outer_line_shift=cat(1,outer_line_shift,outer_line(1:(shift_length-1),:));
else
    outer_line_shift=outer_line;
end

% check for jumps from plane intersections
% distances_outerline=sqrt(sum(abs(diff(outer_line_shift).^2),2));

% jump_tol=50;
% if any(distances_outerline>jump_tol)
%     deletion_indices=find(distances_outerline>jump_tol);
%     delarr=[];
%     while 
%     for delidx=1:2:length(deletion_indices)
%         delarr=cat(2,delarr,deletion_indices(delidx):deletion_indices(delidx+1));
%     end
%     outer_line_shift(delarr,:)=[];
% end

for pt_idx=1:size(outer_line_shift,1)
    all_distances_outer_line_shift(pt_idx)=norm(point-outer_line_shift(pt_idx,:));
end

[~,index]=min(all_distances_outer_line_shift);

%% get short segment of outline

desired_interval=ceil(max(det_distances)/2)+5;
distance_between_pts=sqrt(sum(abs(diff(outer_line_shift)).^2,2));

left_ptsToAdd=1;
left_distance=distance_between_pts(index);
while left_distance<desired_interval
    left_ptsToAdd=left_ptsToAdd+1;
    left_distance=sum(distance_between_pts(index:(index+left_ptsToAdd)));
end

right_ptsToAdd=1;
right_distance=distance_between_pts(index-1);
while right_distance<desired_interval
    right_ptsToAdd=right_ptsToAdd+1;
    right_distance=sum(distance_between_pts((index-right_ptsToAdd):(index-1)));
end

% approx_distance_between_pts=median(sqrt(sum(abs(diff(outer_line_shift)).^2,2)));
% interv=floor(shift_tolerance/approx_distance_between_pts);

[x,xa,~]=unique(outer_line_shift((index-right_ptsToAdd):(index+left_ptsToAdd),1));
y=outer_line_shift((index-right_ptsToAdd):(index+left_ptsToAdd),2);
y=y(xa);
z=outer_line_shift((index-right_ptsToAdd):(index+left_ptsToAdd),3);
z=z(xa);

%% take section of the outer line surrounding the fiducial, and interpolate then smooth it

orig_points=[x,y,z];
allranges=[range(x) range(y) range(z)];
[~,max_idx]=max(allranges);
orig_points=sortrows(orig_points,max_idx);

CS = cat(1,0,cumsum(sqrt(sum(diff(orig_points,[],1).^2,2))));

[~,ia,~]=unique(CS);

% interp_points_jagged = interp1(CS(ia), orig_points(ia,:), unique([CS(ia)' linspace(0,CS(end),1000)]),'spline');
interp_points_jagged = interp1(CS(ia), orig_points(ia,:), unique([CS(ia)' linspace(0,CS(end),1000)]),'pchip');
[interp_points,~]=smoothdata(interp_points_jagged,1,'sgolay');

%% find source and detector locations

sd_sep=det_distances(end)/2;
 
for pt=1:size(interp_points,1)
    distance_from_fiducial(pt)=norm(point-interp_points(pt,:));
end

% ensure subject left-right directionality in probe placement
% [~,source_start_idx]=find(abs(distance_from_fiducial-sd_sep)<0.15,1);
[~,source_start_idx]=min(abs(distance_from_fiducial-sd_sep));
source=interp_points(source_start_idx,:);

% find detector locations
if source_start_idx<length(interp_points)/2
    temp_interp_points=interp_points(source_start_idx:end,:);
else
    temp_interp_points=interp_points(1:source_start_idx,:);
    temp_interp_points=flipud(temp_interp_points);
end

distance_between_dets=sqrt(sum(abs(diff(temp_interp_points).^2),2));

sum_dist=cumsum(distance_between_dets);
for idx=1:length(det_distances)
    [~,index_values(idx)]=min(abs(sum_dist-det_distances(idx)));
end

det_arr=temp_interp_points(index_values,:);
 
%% combining into single array and finding normals for source and each detector

% find source normal
[source_unit_vec,~]=getNormalVec(node,face,source);

% find detector normals
for det_idx=1:length(det_distances)
    [det_unit_vec(det_idx,:),~]=getNormalVec(node,face,squeeze(det_arr(det_idx,:)));
end

%% plot

% if showfig
%     plotSDmesh(node,face,point,source,source_unit_vec,det_arr,det_unit_vec);
% end

