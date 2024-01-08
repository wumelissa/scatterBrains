function g2Noise=generate_dcs_g2_noise(g2,beta,CR,tau,avgT)
%% Computes the noise of the intensity autocorrelation function, g2, based
% on the noise model expressed in "Koppel, D. (1974). Statistical accuracy in FCS. Physical Review A"

% author: Mitchell Robinson (robinson.mitchell.b <at> gmail.com)

% Inputs
% g2: Matrix of g2 curves of size [nTau,nSDS]
% beta: Coherence factor of the measurement, if it is not of size [1,nSDS],
% the same coherence factor value will be used for all measurements
% CR: Array of count rate values in units of photons per second of size
% [1,nSDS]
% tau: Array of autocorrelation lags
% avgT: Averaging time of the measurement

% Outputs
% g2Noise: Matrix of the values of the standard deviation of the g2 curves,   
% is of size [nTau,nSDS]

%% Argument management 
arguments
    g2 double
    beta double
    CR double
    tau (1,:) double
    avgT (1,1) double
end
if ~(length(beta)==size(g2,2))
	beta=beta(1)*ones([1,size(g2,2)]);
end
if ~(size(g2,1)==size(tau,2))
    error('Length of the correlation function array should match the length of the tau array for noise calculation.');
end
if ~(size(g2,2)==numel(CR))
    error('Each correlation function requires a photon count rate for noise calculation.')
end
%% Setting up computations common to all g2s that have been input


% Computing the width of the tau bins
binWidth=diff(tau);
binWidth=[binWidth(1),binWidth];

% Initializing the noise output
g2Noise=zeros(size(g2));
 
%% Loop over each correlation function supplied
for g2Iter=1:size(g2,2)
    
    % First fit the first 50% of the g2 as a single exponential decay
    g2Norm=(g2(:,g2Iter)-min(g2(:,g2Iter)))./range(g2(:,g2Iter));
    [~,ind50]=min(abs(g2Norm-0.5));
    
    % Compute the decay value (slope of the log of the g2 decay)
    gamma=-sum(tau(1:ind50)'.*log(g2Norm(1:ind50)))./sum(tau(1:ind50).^2);

    % Terms for the noise
    CR_array=CR(g2Iter).*binWidth; % average number of photons in each bin
    
    term1=beta(g2Iter)^2.*((1+exp(-2*gamma*binWidth)).*(1+exp(-2*gamma*tau))+2*tau./binWidth.*(1-exp(-2*gamma*binWidth)).*(exp(-2*gamma*tau)))./(1-exp(-2*gamma*binWidth));

    term2=2./CR_array.*beta(g2Iter).*(1+exp(-2*gamma*tau));

    term3=1./CR_array.^2.*(1+beta(g2Iter)*exp(-gamma*tau));

    g2Noise(:,g2Iter)=sqrt(binWidth./avgT).*sqrt(term1+term2+term3);
end