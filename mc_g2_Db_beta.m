function [g2,g1]=mc_g2_Db_beta(Db,beta,tau,expg2,num_dets,mts,exps,g1_norms,...
    num_layers,k0,varargin)

% calculates g2; see lines
% 54-65 of calculate_g2_g1.m function for example

% inputs:
% Db: brownian motion coefficient
% beta: coherence factor
% tau: array of delays
% expg2: experimental g2 (for plotting purposes; can be set to empty)
% num_dets: number of detectors
% mts: cell array of momentum transfers per detector
% exps: cell array of exp(-mua*pl) per detector
% g1_norms: cell array of normalization factors per detector
% num_layers: number of tissue layers
% k0: wavenumber
% varargin: set flag as 1 to plot the curve

% outputs:
% g2: ntau x ndetectors
% g1: ntau x ndetectors

% author: Melissa Wu (wu.melissa.m <at> gmail.com)
% contributing author: Stefan Carp (stefan.carp <at> mgh.harvard.edu)

% this file is part of scatterBrains
% License: GPLv3

%% adding up g2

if ~isempty(varargin), showfig=varargin{1}; else, showfig=0; end

if length(Db)==1, Db=ones(num_layers,1)*Db; end

for I=1:num_dets
    temp_big=exp(-(k0.^2.*2*mts{I}*(Db*tau))).*exps{I};
    g1(:,I)=sum(temp_big)'/g1_norms{I};
    g2(:,I)=1+beta(I)*(g1(:,I).^2);
end
  
g2=gather(g2);
 
 %% plotting
 
 if showfig
     figure(149);
     hold off;
     semilogx(repmat(tau,[1 num_dets]),expg2(:));
     hold on
     semilogx(repmat(tau,[1 num_dets]),g2(:),'r');
     title({sprintf('%0.2f %0.2f %0.2f\n',beta), sprintf('%1.2e %1.2e %1.2e %1.2e\n',Db([1 end]))})
     ylim([0.8 1.7]); grid on; grid minor
     drawnow
 end
 