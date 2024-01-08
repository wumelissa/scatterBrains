function dcs_bfi_noise=compute_dcs_bfi_noise(tau,g2,g2Noise,sds,beta,options)
%% Estimates the BFi noise for a set of correlation functions and the noise properties determined by the count rate.
% author: Mitchell Robinson (robinson.mitchell.b <at> gmail.com)

% Inputs
% tau: An array of autocorrelation lags in seconds
% g2: A matrix of values containing the autocorrelation functions
% g2Noise: A matrix of values containing the standard deviations for values
% in the corresponding g2 matrix
% mua: Mua used for fitting in the semi-infinite model in units of
% mm^-1
% msp: Musp used for fitting in the semi-infinite model in units of
% mm^-1
% wv: Wavelength used for fitting in units of mm
% sds: source-detector separation used for fitting in units of mm
% n: Index of refraction used for fitting
% beta: Coherence parameter of the measurement
% fitting_range: Range of the proportion of the g2 decay that is fit
% lb: lower bounds for the fitting for parameters [bfi] in units of
% [mm^2/s]*1e6
% ub: upper bounds for the fitting for parameters [bfi] in units of
% [mm^2/s]*1e6
% fit_options: Options used for the fitting function
% nObs: number of noisy observations for each measurement 
% snr_limit: Limit of the SNR that is acceptable for fitting. This
% eliminates the BFi "spiking" that will occur at low SNR.
% noise_range: Region of the g2 curve that is used for the SNR estimation



% Outputs
% dcs_bfi_noise: Estimated standard deviation of BFi for each input

%% Handling input arguments
arguments
tau double
g2 double
g2Noise double
sds double
beta double
options.mua (1,1) double = 0.017
options.msp (1,1) double = 0.80
options.wv (1,1) double = 850e-6
options.n (1,1) double =1.4
options.bfi_lb (1,1) double = 0
options.bfi_ub (1,1) double = 20
options.g2_lb (1,1) double = 0.01
options.g2_ub (1,1) double = .999
options.snr_limit (1,1) double = 20
options.nObs (1,1) double = 500
options.noise_range (1,2) double = [0.99,0.9];
options.lsqCurveFitOptions = optimoptions('lsqcurvefit','display','none','OptimalityTolerance',1e-12,...
        'FunctionTolerance',1e-12,'StepTolerance',1e-12);
end

if ~(size(g2,1)==length(tau))
    error('Tau array and g2 array are not the same length.');
end
if ~(size(g2Noise,1)==length(tau))
    error('Tau array and g2 noise array are not the same length.');
end
if ~(size(g2,2)==size(g2Noise,2))
    error('g2 array and g2 noise array do not have the same number of source-detector separations.');
end
if ~(size(g2,2)==length(sds))
    error('Number of correlation functions is not equal to the number of given source-detector separations.');
end
%% Internal function settings used for estimating noise
mua=options.mua;
msp=options.msp;
wv=options.wv;
n=options.n;
bfi_lb=options.bfi_lb;
bfi_ub=options.bfi_ub;
g2_lb=options.g2_lb;
g2_ub=options.g2_ub;
nObs=options.nObs;
noise_range=options.noise_range;
snr_limit=options.snr_limit;
fit_options=options.lsqCurveFitOptions;


%% Performing the BFi noise estimation
dcs_bfi_noise=zeros(size(sds));
for sdsIter=1:length(dcs_bfi_noise)
    % Select the correlation function, noise, and source-detector separation
    sel_g2=g2(:,sdsIter);
    sel_noise=g2Noise(:,sdsIter);
    sel_sds=sds(sdsIter);
    
    % Determine the fitting indices
    [~,fitI(1)]=min(abs(sel_g2-1-g2_ub));
    [~,fitI(2)]=min(abs(sel_g2-1-g2_lb));
    fitInds=fitI(1):fitI(2);
    
    % Determine the noise estimation indices
    [~,noiseI(1)]=min(abs(sel_g2-1-noise_range(1)));
    [~,noiseI(2)]=min(abs(sel_g2-1-noise_range(2)));
    noiseInds=noiseI(1):noiseI(2);
    
    % Compute the snr scaling required to reach the desired noise
    % properties of the g2 curve. This step would be equivalent to
    % increasing the averaging time by a factor of snr_scale^2
    snr_scale=snr_limit./mean(beta*(sel_g2(noiseInds)-1)./sel_noise(noiseInds));
    
    % Scale the selected noise by the computed SNR scale
    sel_noise=sel_noise./snr_scale;

    % Generate many noisy observations of the correlation function and fit
    % them
    noisy_g2s=1+beta*(sel_g2(fitInds)-1)+(sel_noise(fitInds).*randn([length(fitInds),nObs]));
    noisy_bfi=zeros([1,nObs]);
    
    for obIter=1:nObs
        noisy_bfi(obIter)=g2_fit_fixBeta(tau(fitInds),noisy_g2s(:,obIter)',...
            sel_sds,beta,'mua',mua,'msp',msp,'wv',wv,'n',n,...
            'bfi_lb',bfi_lb,'bfi_ub',bfi_ub,'lsqCurveFitOptions',fit_options,...
            'g2_lb',-inf,'g2_ub',inf);
    end
    dcs_bfi_noise(sdsIter)=std(noisy_bfi)*snr_scale;
end
end