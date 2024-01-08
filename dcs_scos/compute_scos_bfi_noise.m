function scos_bfi_noise=compute_scos_bfi_noise(expT,kf2,kf2_noise,sds,beta,options)
%% Estimates the BFi noise for a set of correlation functions and the noise properties determined by the count rate.
% author: Mitchell Robinson (robinson.mitchell.b <at> gmail.com)

% Inputs
% expT: The exposure time of the simulated measurement
% kf2: An array of values containing the fundamental, squared speckle
% contrast
% kf2_noise: An array of values containing the standard deviation of the
% corresponding values in the kf2 array
% fit_mua: Mua used for fitting in the semi-infinite model in units of
% mm^-1
% fit_msp: Musp used for fitting in the semi-infinite model in units of
% mm^-1
% fit_wv: Wavelength used for fitting in units of mm
% fit_sds: source-detector separation used for fitting in units of mm
% fit_n: Index of refraction used for fitting
% beta: Coherence parameter of the measurement
% fitting_range: Range of the proportion of the g2 decay that is fit
% lb: lower bounds for the fitting for parameters [bfi] in units of
% [mm^2/s]*1e6
% ub: upper bounds for the fitting for parameters [bfi] in units of
% [mm^2/s]*1e6
% lsqCurveFitOptions: Options used for the fitting function
% snr_limit: Scaling factor for SNR to ensure a reasonable estimate. 
% Will not allow for the accounting of BFi spikes caused by fitting 
% errors in favor of a more stable estimate.

% Outputs
% scos_bfi_noise: Estimated standard deviation of BFi for each input

%% Handling input arguments
arguments
expT (1,1) double
kf2 double
kf2_noise double
sds double
beta double
options.mua (1,1) double = 0.017
options.msp (1,1) double = 0.80
options.wv (1,1) double = 850e-6
options.n (1,1) double =1.4
options.bfi_lb (1,1) double = 0
options.bfi_ub (1,1) double = 20
options.snr_limit (1,1) double = 20
options.lsqCurveFitOptions = optimoptions('lsqcurvefit','display','none','OptimalityTolerance',1e-12,...
        'FunctionTolerance',1e-12,'StepTolerance',1e-12);
end

if or(~(length(kf2)==length(kf2_noise)),~(length(kf2)==length(sds)))
    error('The lengths of the kf2 array, the kf2 noise array, and the source-detector separation array should be equal.')
end

%% Simplifying the options inputs
mua=options.mua;
msp=options.msp;
wv=options.wv;
n=options.n;
bfi_lb=options.bfi_lb;
bfi_ub=options.bfi_ub;
snr_limit=options.snr_limit;
fit_options=options.lsqCurveFitOptions;

%% Performing the BFi noise estimation
scos_bfi_noise=zeros(size(sds));
for sdsIter=1:length(sds)
    kf2_sel=kf2(sdsIter);
    noise_sel=kf2_noise(sdsIter);
    sel_sds=sds(sdsIter);
    
    % Scale the kf2 based on the beta of the measurement
    kf2_sel=kf2_sel*beta;
    
    % Compute the snr scale
    snr_scale=max([snr_limit/(kf2_sel/noise_sel),1]);
    
    noise_sel=noise_sel/snr_scale;
    
    % Plot the distribution of scos values
    kf2_range=(-snr_limit/2:.1:snr_limit/2)*noise_sel+kf2_sel;
    kf2_prob=exp(-(kf2_range-kf2_sel).^2./(2*noise_sel^2));
    
    % Fit the 
    bfi_range=zeros(size(kf2_range));
    
    for i=1:length(bfi_range)
         bfi_range(i)=kf2_fit_fixBeta(expT,kf2_range(i),sel_sds,beta,'mua',mua,'msp',msp,'wv',wv,'n',n,...
             'bfi_lb',bfi_lb,'bfi_ub',bfi_ub,'lsqCurveFitOptions',fit_options);
    end
    
    bfi_range=bfi_range((length(bfi_range)+1)/2)+(bfi_range-bfi_range((length(bfi_range)+1)/2))*snr_scale;
    bfi_m=sum(bfi_range.*kf2_prob)./sum(kf2_prob);
    bfi_std=sqrt(sum((bfi_range-bfi_m).^2.*kf2_prob)./sum(kf2_prob));
    scos_bfi_noise(sdsIter)=bfi_std;
end

end