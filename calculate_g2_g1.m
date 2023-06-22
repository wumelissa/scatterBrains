function [g2,g1,tau]=calculate_g2_g1(fhist,Db,tau,lambda,max_photons,varargin)

% this function calculates g2 and g1 

% inputs:
% fhist: history filename (.mch file)
% Db: brownian motion coefficient, m^2/s
% lambda: wavelength, in m
% max_photons: max number of photons to be stored per detector

% outputs:
% g2: array of (ntau x ndetectors)
% g1: array of (ntau x ndetectors)
% tau: ntau

% author: Melissa Wu (wu.melissa.m <at> gmail.com)
% contributing author: Stefan Carp (stefan.carp <at> mgh.harvard.edu)

% this file is part of scatterBrains
% License: GPLv3

if nargin<5, max_photons=1e5; end
if nargin<4, lambda=785e-6; end
if nargin<3, tau=logspace(-8,0,200); end
if nargin<2, Db=6e-6; end

%% load photon history and optical properties

% here
his_temp=getPhotonHistory(fhist);

temp=strfind(fhist,filesep);
if isempty(temp), fhist=[pwd filesep fhist]; end
temp=strfind(fhist,filesep);
lastslash=temp(end);
sim_label=fhist(lastslash+1:end-4);
[mua,~,~,n]=load_mc_prop([fhist(1:lastslash) filesep 'prop_' sim_label '.dat']);

%% reordering photon history data

photon_counts=histcounts(his_temp(:,1),'BinMethod','integer');
photon_counts(photon_counts>max_photons)=max_photons;
photon_indices(1,:)=[1 photon_counts(1)];

for det_idx=2:length(photon_counts)
    photon_indices(det_idx,:)=[1+sum(photon_counts(1:(det_idx-1))) sum(photon_counts(1:det_idx))];
end

his_data=zeros(sum(photon_counts),size(his_temp,2));
for det_idx=1:length(photon_counts),
    allphot_idx=find(his_temp(:,1)==(det_idx),max_photons);
    curr_idxs=photon_indices(det_idx,1):photon_indices(det_idx,2);
    his_data(curr_idxs,:)=his_temp(allphot_idx,:);
end

num_tissue_layers=(size(his_data,2)-1)/2;

%% summing across tissue layers

num_dets=length(photon_counts);

for detector=1:num_dets
    idx_start=photon_indices(detector,1);
    idx_end=photon_indices(detector,2);
    mtransfer=his_data(idx_start:idx_end,(num_tissue_layers+2):end);
    path_length=his_data(idx_start:idx_end,2:(num_tissue_layers+1));
    exp1=repmat(exp(-mua*path_length'),[length(tau) 1])';
    g1_norm=sum(exp(-mua*path_length'));
    mts{detector}=mtransfer;
    exps{detector}=exp1;
    g1_norms{detector}=g1_norm;
    fprintf('Detector %d: %d photons\n',detector, idx_end-idx_start)
end

%% calculate g2 + g1

beta=0.5*ones(num_dets,1);
k0=2*pi*n/lambda;
k0=k0(1);
expg2=0;

fprintf('Calculating g1...\n')
[g2,g1]=mc_g2_Db_beta(Db,beta,tau,expg2,num_dets,mts,exps,g1_norms,...
    num_tissue_layers,k0);
fprintf('Done.\n')