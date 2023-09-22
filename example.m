
%% add paths

% make sure to install iso2mesh version 1.9.0 (century egg) and mmc 
% version 1.9 (v2020, Moon Cake - beta) in the same folder, OR edit below
% variables (iso2mesh_path and mmc_path) to point to software downloads

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

%% read subject mesh

% we will work with subject 3 in this example. here we load the mesh

subject_num='03';

node=readmmcnode([scatterBrains_path filesep 'Subject' subject_num filesep 'node_subject' subject_num '.dat']);
elem=readmmcelem([scatterBrains_path filesep 'Subject' subject_num filesep 'elem_subject' subject_num '.dat']);
face=readmmcface([scatterBrains_path filesep 'Subject' subject_num filesep 'face_subject' subject_num '.dat']);
load([scatterBrains_path filesep 'Subject' subject_num filesep 'Subject' subject_num '_locs.mat'])

%% visualize all of the points

% the below code allows you to visualize all of the points and select one interactively. 

[user_point,user_pt_idx]=get_surface_loc(face,node,locs);

%% use example point

% for this example we'll use point 125

pt_idx=125;
point=locs(pt_idx,:);
fprintf('Extracerebral thickness at this point is %.2f mm\n',distance2brain(pt_idx))

figure(149)
trimesh(face(:,1:3),node(:,1),node(:,2),...
    node(:,3),'facecolor','none','FaceAlpha',0.5,'EdgeColor',[192 192 192]/256)
hold on; axis equal
plot3(locs(:,1),locs(:,2),locs(:,3),'k.','MarkerSize',20)
plot3(point(1),point(2),point(3),'r.','MarkerSize',20)

%% wrap the probe

% wrap the probe with detector distances of choice. The original point will
% be halfway in between the source and the maximum detector distance

det_distances=[5 10 25 30];

[source,det_arr,source_unit_vec,~]=wrapProbe(node,face,point,rz,lz,det_distances);

% original point is black, detectors are blue, source is red, source launch vector is red

figure(150)
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

% In addition to detector locations, the radius of each detector can also
% be specified. Detectors closer to the source will detect more photons, so
% sometimes it helps to increase the detector radius for those that are
% further away from the source, to balance the number of photons across
% detectors. However, a detector size too large will impact its "true"
% distance from the source. Here we set each detector radius to 1 mm but
% the user can edit as needed; approximately 5000 photons are needed for a
% reasonably accurate g1 computation, if doing DCS

% edit detector radiii
detector_radii=ones(size(det_arr,1),1);

% define data structure
data.Domain.MeshID=['subject' subject_num];
data.Domain.InitElem=eid;
data.Session.ID=['subject' subject_num];
data.Session.Photons=1e8; % photons launched, user edit
data.Session.Seed=567621849; % simulation seed, user edit
data.Forward.T0=0; % time start, user edit
data.Forward.T1=9e-9; % time end, user edit
data.Forward.Dt=9e-9; % time step, user edit
data.Optode.Source.Pos=source;
data.Optode.Source.Dir=source_unit_vec;
for I=1:size(det_arr,1)
    data.Optode.Detector(I).Pos=squeeze(det_arr(I,:));
    data.Optode.Detector(I).R=detector_radii(I); 
end

% write to json
str=jsonencode(data,'PrettyPrint',true);
fileID=fopen([scatterBrains_path filesep 'Subject' subject_num filesep 'example.json'],'w');
fwrite(fileID,str);
fclose(fileID);

%% write optical properties file

% each column corresponds to a tissue layer, i.e. mua(1) has the absorption
% coefficient for tissue 1, mua(2) has the optical property for tissue
% 2, etc. All of the below quantities are in mm-1. The example properties
% listed here are simplified to demonstrate accurate BFi recovery with the
% semi-infinite model below.
%
% mua: absorption coefficient
% mus: scattering coefficient (not reduced)
% g: anisotropy
% n: refractive index

% the tissue indices correspond to:
% 1: scalp
% 2: skull
% 3: cerebrospinal fluid
% 4: grey matter
% 5: white matter

mua=ones(1,5)*0.015;
mus=ones(1,5);
g=ones(1,5)*0.01;
n=ones(1,5)*1.37;

% note to user: below are more realistic optical properties from Wu et al, Biomed Opt Exp 2022
% however these will not guarantee accurate BFi recovery when fitting with
% the semi-infinite model
% mua=[0.0164 0.0115 0.0017 0.017 0.017];
% mus=[0.7475 0.8182 0.0101 1.1717 1.1717];
% g=[0.01 0.01 0.01 0.01 0.01];
% n=[1.37 1.37 1.37 1.37 1.37];

write_prop_file(mua,mus,g,n,[scatterBrains_path filesep 'Subject' subject_num],['subject' subject_num])

%% run mmc

% please run mmc with your method of choice. 
% If you are working with DCS, be sure to turn the "savedet"
% and "momentum" flags. My typical command line mmc run looks something
% like this (make sure the paths are correct):

% <mmc_executable> -f example.json -D P --momentum 1 -b 1 -d 1
    
%% calculate g1, g2

% calculate g1 and g2 from the mmc output file. Only the first argument
% (history filename) is necessary; the other arguments will default to the
% values shown below if not inputted

% calculation options
Db=6e-6; % diffusion coefficient, mm^2/s
tau=logspace(-8,0,200); % delays, seconds
lambda=850e-6; % wavelength, mm
max_photons=1e5; % photon detection cap to save computation
beta=0.5; % coherence factor in Siegert relation; one for each detector or one value for all detectors

% calculate g2, g1
[g2,g1,tau]=calculate_g2_g1([scatterBrains_path filesep 'Subject' subject_num filesep 'subject' subject_num '.mch'],Db,tau,lambda,max_photons,beta);

% plot
figure(151)
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
% the input to semi_infinite_g2 is [beta, BFi]. Units for BFi are mm^2/s

% semi-infinite model parameters
fit_options.mu_a=0.015; % absorption, mm-1
fit_options.mu_s=1; % reduced scattering, mm-1
fit_options.lambda=850e-6; % mm-1
fit_options.alpha=1; % probability of scattering from moving scatterer
fit_options.n=1.37; % tissue index of refraction

% fitting parameters
x0=[0.5 6e-6]; % initial guesses
lb=[0 6e-8]; % lower bounds
ub=[0.6 6e-4]; % upper bounds
opts = optimset('Display','off'); % fitting options

% calculate "true" distances based on Euclidean distance metric
for detidx=1:size(det_arr,1)
    rhos(detidx)=norm(source-det_arr(detidx,:));
end

% loop through detectors and fit
for detidx=1:size(g2,2)
    fit_options.rho=rhos(detidx);
    [fitted_values, ~] = lsqcurvefit(@(x,taus)semi_infinite_g2(x,tau,fit_options),x0,tau,g2(:,detidx),lb,ub,opts);
    beta(detidx)=fitted_values(1);
    BFi(detidx)=fitted_values(2);
    figure(100)
    computed_g2=semi_infinite_g2(fitted_values,tau,fit_options);
    figure(152);
    semilogx(tau,computed_g2)
    hold on
    semilogx(tau,g2(:,detidx),'--')
    grid on; grid minor
    ylabel('g_2')
    xlabel('\tau')
    title([num2str(rhos(detidx),'%.1f') ' mm'])
    hold off
    pause(0.5)
end