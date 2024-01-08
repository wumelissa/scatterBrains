function [total_distance,interp_points]=find_spherical_distance(node,face,point,point2,refpoint)

% find distance along head surface between two points

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

%% find outline of intersection of plane + volume

% taken from iso2mesh demo code - start
[cut_pos,~,cut_edges]=qmeshcut(face,node,node(:,1),[point;point2;refpoint]);
[cut_pos,cut_edges]=removedupnodes(cut_pos,cut_edges);
cut_loop=extractloops(cut_edges);

cut_loop(isnan(cut_loop))=[];
% taken from iso2mesh demo code - end  

outer_line=cut_pos(cut_loop,:);

outer_line=unique(outer_line,'rows','stable');

%% find closest outline point to POI

distance_to_start=norm(point-outer_line(1,:));

for I=1:size(outer_line,1)
    distance_to_points(I)=norm(point-outer_line(I,:));
end

[~,closestLinePt]=min(distance_to_points);

%% adjust "starting point" of the outline so it's far from POI

euclidean=norm(point-point2);
shift_tolerance=max(euclidean);
shift_length=floor(size(outer_line,1)/2);

% shift if too close to fiducial
if distance_to_start<shift_tolerance || closestLinePt <shift_tolerance || (size(outer_line,1) - closestLinePt)< shift_tolerance
    outer_line_shift=outer_line(shift_length:end,:);
    outer_line_shift=cat(1,outer_line_shift,outer_line(1:(shift_length-1),:));
else
    outer_line_shift=outer_line;
end

for pt_idx=1:size(outer_line_shift,1)
    all_distances_outer_line_shift(pt_idx)=norm(point-outer_line_shift(pt_idx,:));
    all_distances_outer_line_shift_pt2(pt_idx)=norm(point2-outer_line_shift(pt_idx,:));
end

[~,index1]=min(all_distances_outer_line_shift);
[~,index2]=min(all_distances_outer_line_shift_pt2);

%% get short segment of outline

[x,xa,~]=unique(outer_line_shift(index1:index2,1));
if isempty(x)
    [x,xa,~]=unique(outer_line_shift(index2:index1,1));
    y=outer_line_shift(index2:index1,2);
    y=y(xa);
    z=outer_line_shift(index2:index1,3);
    z=z(xa);
else
    y=outer_line_shift(index1:index2,2);
    y=y(xa);
    z=outer_line_shift(index1:index2,3);
    z=z(xa);
end

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

CS2 = cat(1,0,cumsum(sqrt(sum(diff(interp_points,[],1).^2,2))));
total_distance=CS2(end);
