function [g2,g1,tau]=calculate_g2_g1(history_file,options)

% this function calculates g2 and g1 

% inputs:
% history_file: history filename (.mch file)
% Db: brownian motion coefficient, mm^2/s
% tau: array of delays, seconds
% lambda: wavelength, in mm
% max_photons: max number of photons to be stored per detector
% beta: coherence factor, one for each detector or one value for all
% detectors

% outputs:
% g2: array of (ntau x ndetectors)
% g1: array of (ntau x ndetectors)
% tau: array of delays

% author: Melissa Wu (wu.melissa.m <at> gmail.com)
% contributing author: Stefan Carp (stefan.carp <at> mgh.harvard.edu)

% this file is part of scatterBrains

% check arguments and define variables
arguments
    history_file char
    options.Db (1,:) double = 6e-6 % mm^2/s
    options.tau (1,:) double = logspace(-8,0,200)
    options.lambda (1,1) double = 850e-6 % mm
    options.max_photons (1,1) double = 1e5
    options.beta (1,1) double = 0.5
end

Db=options.Db;
tau=options.tau;
lambda=options.lambda;
max_photons=options.max_photons;
beta=options.beta;

%% load photon history and optical properties

his_temp=getPhotonHistory(history_file);

tmp=strfind(history_file,filesep);
if isempty(tmp)
    history_file=[pwd filesep history_file]; 
end
tmp=strfind(history_file,filesep);
lastslash=tmp(end);
sim_label=history_file(lastslash+1:end-4);

[mua,~,~,n]=load_mc_prop([history_file(1:lastslash) filesep 'prop_' sim_label '.dat']);

%% getting detector indices

num_dets=length(unique(his_temp(:,1)));

for det_idx=1:num_dets
    photon_indices{det_idx}=find(his_temp(:,1)==(det_idx),max_photons);
end

num_tissue_layers=(size(his_temp,2)-1)/2;

%% calculate g2 + g1

if length(beta)==1
    beta=beta*ones(num_dets,1);
end
k0=2*pi*n/lambda;

fprintf('Calculating g1...\n')
[g2,g1]=mc_g2_Db_beta(his_temp,photon_indices,Db,beta,tau,num_dets,num_tissue_layers,k0,mua);
fprintf('Done.\n')