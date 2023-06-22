
%% add paths

% make sure to install iso2mesh version 1.9.0 (century egg) and mmc 
% version 1.9 (v2020, Moon Cake - beta) in the same folder.

% author: Melissa Wu (wu.melissa.m <at> gmail.com)

iso2mesh_path=['.' filesep 'iso2mesh'];
mmc_path=['.' filesep 'mmc'];

addpath(genpath(iso2mesh_path))
addpath(genpath(mmc_path))

%% read subject mesh

% we will work with subject 3 in this example. here we load the mesh

subject_num='03';

node=readmmcnode(['.' filesep 'Subject' subject_num filesep 'node_subject' subject_num '.dat']);
elem=readmmcelem(['.' filesep 'Subject' subject_num filesep 'elem_subject' subject_num '.dat']);
face=readmmcface(['.' filesep 'Subject' subject_num filesep 'face_subject' subject_num '.dat']);
load(['.' filesep 'Subject' subject_num filesep 'Subject' subject_num '_locs.mat'])

%% visualize all of the points

% we will use point 125 in this example. we can visualize it in red

figure
trimesh(face(:,1:3),node(:,1),node(:,2),...
    node(:,3),'facecolor','none','FaceAlpha',0.5,'EdgeColor',[192 192 192]/256)
hold on; axis equal
plot3(locs(:,1),locs(:,2),locs(:,3),'k.','MarkerSize',20)
plot3(locs(125,1),locs(125,2),locs(125,3),'r.','MarkerSize',20)

%% wrap the probe

% wrap the probe with detector distances of choice

point=locs(125,:);
det_distances=[5 10 25 50];

[source,det_arr,source_unit_vec,~]=wrapProbe(node,face,point,rz,lz,det_distances);

% original point is black, detectors are blue, source is red, source launch vector is red

figure
trimesh(face,node(:,1),node(:,2),node(:,3),...
    'facecolor','none','FaceAlpha',0.5,'EdgeColor',[192 192 192]/256)
hold on
axis equal
plot3(point(1),point(2),point(3),'ko','LineWidth',1,'MarkerFaceColor','k')
plot3(source(1),source(2),source(3),'ro','LineWidth',1,'MarkerFaceColor','r')
plot3(det_arr(:,1),det_arr(:,2),det_arr(:,3),'bo','LineWidth',1,'MarkerFaceColor','b')
quiver3(source(1),source(2),source(3),source_unit_vec(1)*30,...
    source_unit_vec(2)*30,source_unit_vec(3)*30,'LineWidth',2,'Color','r')

%% ensure source is in mesh

% mmc requires the mesh element number in which the source is enclosed.
% sometimes the source doesn't fall perfectly in the mesh element such that
% mmc is able to detect it. In these cases, we can move the source
% incrementally by 0.01 mm steps until mmc can detect it.

eid=tsearchn(node(:,1:3),elem(:,1:4),source);
increment=0.01;
stretch=0;
while isnan(eid)
    stretch=stretch+increment;
    eid=tsearchn(node(:,1:3),elem(:,1:4),source+stretch*source_unit_vec);
end

source=source+stretch*source_unit_vec;

%% write to json file

% now we are ready to run mmc. this cell will generate an input json file.
% Please see mmc documentation for all of the input parameter descriptions.

% define data structure
data.Domain.MeshID=['subject' subject_num];
data.Domain.InitElem=eid;
data.Session.ID=['subject' subject_num];
data.Session.Photons=1e7; % photons launched, user edit
data.Session.Seed=567621849; % simulation seed, user edit
data.Forward.T0=0; % time start, user edit
data.Forward.T1=9e-9; % time end, user edit
data.Forward.Dt=9e-9; % time step, user edit
data.Optode.Source.Pos=source;
data.Optode.Source.Dir=source_unit_vec;
for I=1:size(det_arr,1)
    data.Optode.Detector(I).Pos=squeeze(det_arr(I,:));
    data.Optode.Detector(I).R=1; % detector radii, user edit
end

% write to json
str=jsonencode(data,'PrettyPrint',true);
fileID=fopen(['.' filesep 'Subject' subject_num filesep 'example.json'],'w');
fwrite(fileID,str)
fclose(fileID);

%% write optical properties file

% each column corresponds to a tissue layer, i.e. mua(1) has the absorption
% coefficient for tissue 1, mua(2) has the optical property for tissue
% 2, etc. All of the below quantities are in mm-1.
%
% mua: absorption coefficient
% mus: scattering coefficient
% g: anisotropy
% n: refractive index

mua=[0.0164 0.0115 0.0017 0.017 0.017];
mus=[0.7475 0.8182 0.0101 1.1717 1.1717];
g=[0.01 0.01 0.01 0.01 0.01];
n=[1.37 1.37 1.37 1.37 1.37];

write_dat_file(mua,mus,g,n,['subject' subject_num])

%% run mmc

% please run mmc with your method of choice. 
% If you are working with DCS, be sure to turn the "savedet"
% and "momentum" flags. My typical command line mmc run looks something
% like this (make sure the paths are correct):

% <mmc_executable> -f example.json -D P --momentum 1 -b 1 -d 1 -x 1
    
%% calculate g1, g2

% calculate g1 and g2 from the mmc output file.

[g2,g1,tau]=calculate_g2_g1(['.' filesep 'Subject' subject_num filesep 'subject' subject_num '.mch']);

figure
subplot(121)
semilogx(tau,g1)
grid on; grid minor
xlabel('\tau'); ylabel('g_1')
title('g_1')
subplot(122)
semilogx(tau,g2)
xlabel('\tau'); ylabel('g_2')
title('g_2')
grid on; grid minor
for I=1:length(det_distances); leg_arr{I}=[num2str(det_distances(I)) ' mm'];end
legend(leg_arr)

%% fit the data

% fit the g2 using the semi-infinite diffusion correlation equation.
% the input to semi_infinite_g2 is [beta, BFi*1e9]; the BFi is multiplied
% by a factor of 1e9 for a more stable fit. Thus, the output is divided by
% 1e9.

fit_options.mu_a=0.015; % mm-1
fit_options.mu_s=1; % mm-1
fit_options.lambda=785e-6; % mm-1
fit_options.alpha=1; % probability of scattering from moving scatterer
fit_options.n=1.37; % tissue index of refraction

% calculate "true" distances based on Euclidean distance metric
for detidx=1:size(det_arr,1)
    rhos(detidx)=norm(source-det_arr(detidx,:));
end

for detidx=1:size(g2,2)
    fit_options.rho=rhos(detidx);
    [fitted_values, ~] = lsqcurvefit(@(x,taus)semi_infinite_g2(x,tau,fit_options),[.5 1e3],tau,g2(:,detidx));
    beta(detidx)=fitted_values(1);
    BFi(detidx)=fitted_values(2)/1e9;
    figure(100)
    computed_g2=semi_infinite_g2(fitted_values,tau,fit_options);
    figure(100);
    semilogx(tau,computed_g2)
    hold on
    semilogx(tau,g2(:,detidx),'--')
    grid on; grid minor
    ylabel('g_2')
    xlabel('\tau')
    title([num2str(rhos(detidx)) ' mm'])
    hold off
    pause(0.5)
end