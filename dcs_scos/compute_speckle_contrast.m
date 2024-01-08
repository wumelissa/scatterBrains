function kf2=compute_speckle_contrast(g1,tau,expT)
%% Computes the fundamental speckle contrast

% author: Mitchell Robinson (robinson.mitchell.b <at> gmail.com)

% Inputs
% g1: Matrix of g1 curves of size [nTau,nSDS]
% tau: Array of autocorrelation lags of size [1,nTau]
% expT: Exposure time used to make the contrast measurement

% Outputs
% kf2: Fundamental squared speckle contrast of size [1,nSDS]

arguments
    g1 double
    tau double
    expT (1,1) double
end

if ~(size(g1,1)==size(tau,2))
    error('Length of the correlation function array should match the length of the tau array for speckle contrast calculations');
end
% Determine whether the exposure time is a member of the tau array
if sum(tau==expT)==0
    % Add the exposure time into the tau array
    tauNew=sort([tau,expT]);
    g1New=zeros([length(tauNew),size(g1,2)]);
    for g1_iter=1:size(g1,2)
        g1New(:,g1_iter)=interp1(tau,g1(:,g1_iter),tauNew);
    end
    
    % Replace the tau and g1 arrays
    tau=tauNew;
    g1=g1New;
end

% Find the index for the exposure time
[~,exptInd]=min(abs(tau-expT));

% Separate the parts of the tau and g1 arrays within the exposure time
selIndex=1:exptInd;
selTau=tau(selIndex);
selg1=g1(selIndex,:);

% Generate the integration ramp function
integrationRamp=1-selTau/expT;

% Compute the square of the speckle contrast
kf2=2./expT*trapz(selTau',selg1.^2.*integrationRamp',1);
end