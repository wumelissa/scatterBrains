function [g2,g1]=mc_g2_Db_beta(Db,beta,tau,num_dets,...
    num_layers,k0,mua,his_data,photon_indices,varargin)

% inputs:
% Db: brownian motion coefficient mm^2/s
% beta: coherence factor
% tau: array of delays in seconds
% num_dets: number of detectors
% num_layers: number of tissue layers
% k0: wavenumber mm-1
% mua: absorption coefficient mm-1, one for each tissue layer
% his_data: photon history array
% photon_indices: indices of the rows corresponding to each detector
% varargin: set flag as 1 to plot the curve

% outputs:
% g2: ntau x ndetectors
% g1: ntau x ndetectors

% author: Melissa Wu (wu.melissa.m <at> gmail.com)
% contributing author: Stefan Carp (stefan.carp <at> mgh.harvard.edu)

% this file is part of scatterBrains

%% calculating g1, g2

if ~isempty(varargin), expg2=varargin{1}; else, expg2=[]; end

if length(Db)==1, Db=ones(1,num_layers)*Db; end

for I=1:num_dets
    mtransfer=his_data(photon_indices{I},(num_layers+2):end);
    path_length=his_data(photon_indices{I},2:(num_layers+1));
    parfor J=1:length(tau)
        rmsdisp=6*Db.*tau(J);
        g1(J,I)=sum(exp(-(k0.^2.*rmsdisp/3)*mtransfer'-mua*path_length'));
    end
    g1_norm=double(sum(exp(-mua*path_length')));
    g1(:,I)=g1(:,I)/g1_norm;
    g2(:,I)=1+beta(I)*g1(:,I).^2;
    fprintf(['Detector %d out of %d finished: %d photons detected\n'],I,num_dets,size(mtransfer,1))
end
 
 %% plotting
 
 if ~isempty(expg2)
     figure(149);
     hold off;
     semilogx(repmat(tau,[1 num_dets]),expg2(:));
     hold on
     semilogx(repmat(tau,[1 num_dets]),g2(:),'r');
     title({sprintf('%0.2f %0.2f %0.2f\n',beta), sprintf('%1.2e %1.2e %1.2e %1.2e\n',Db([1 end]))})
     ylim([0.8 1.7]); grid on; grid minor
     drawnow
 end
 