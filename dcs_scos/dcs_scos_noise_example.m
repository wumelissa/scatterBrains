%% Full pipeline to take photon trajectory data generated from Monte Carlo simulation 
%  to estimates of DCS and SCOS measurement characteristics.
% author: Mitchell Robinson (robinson.mitchell.b <at> gmail.com)

%% 0. Run MMC simulation to generate the diffuse reflectance profile and the history of detected photons
% An example of the mmc command after running the example script in the
% base scatterbrains implementation when in the folder of the .json file would be:
% MMC_DIRECTORY\mmc.exe -f example.json -n 1e8 -D P -b 1 -d 1 --saveref 1 --momentum 1
%% 1. Setting up the simulation parameters

% Add paths to main scatter brains, mmc, and isotomesh
mmcPath=['..',filesep,'mmc',filesep];
scatterbrainsPath='..';
iso2meshPath=['..',filesep,'iso2mesh'];

addpath(genpath(mmcPath));
addpath(genpath(scatterbrainsPath));
addpath(genpath(iso2meshPath));

% Data directory
data_directory=['..',filesep,'Subject03',filesep];

% History file
mch_history_file=[data_directory,'subject03.mch'];

% g2 Settings
tauRange=[1e-9,1e-1]; % s
tauN=501;

% Count rate reference
sds_reference=25; % mm
cr_reference=10e3; % counts per second per mode
ansi_reference=38; % mW
wv_reference=850; % nm

% Laser settings
beta0=1; % Coherence factor of the laser before sending to the tissue (value between 0 and 1)
nPolarizations=2; % Number of polarizations collected (value between 1 and 2)
laserMaxPower=300; % mW
nSources=1; % Gives the option to simulate a multi-source arrangement

% DCS detector settings
dcsQE=0.55; 
dcsNFiberModes=1; % Number of spatial fiber modes

% SCOS detector settings
scosQE=0.20;
scosFiberModes=1.5e7; % Determined by the v-number of the fiber used as well as the number of fibers in a bundle
scosPixelNumber=2.5e6;
scosReadNoise=2; % e-
scosSPRatio=1;
scosMaxFrameRate=150; % Hz
scosExpT=1e-3; % s

% Measurement conditions
fs_bfi=10; % Hz
tissueBFi=[1e-6,1e-8,1e-6,6e-6,6e-6]; % mm^2/s, given for the 5 tissue layers (scalp, skull, CSF, gray matter, white matter)
brainPerturb=0.25; % Fractional change in brain BFi (gray and white matter)
sensitivityPerturbation=[1,1,1,1+brainPerturb,1+brainPerturb];

% Fitting settings for semi-infinite fitting model
fit_mua=0.017; % mm^-1 
fit_msp=0.80; % mm^-1
fit_n=1.4;
fit_range=[.9999,.05]; % Values of normalized g2 to include in the fit

fprintf('Section 1 complete.\n');
%% 2. Generate the variables based on the settings given above

if tauRange(2)<scosExpT
    tauRange(2)=scosExpT;
end
tau=[0,logspace(log10(tauRange(1)),log10(tauRange(2)),tauN)];

% Calculation of DCS beta value
dcsNModes=nPolarizations*dcsNFiberModes;
dcsBeta=beta0/dcsNModes;

% Calculation of SCOS beta value
if scosSPRatio<sqrt(scosPixelNumber/scosFiberModes)
    scosSPRatio=sqrt(scosPixelNumber/scosFiberModes);
end
scosNModes=(1+scosSPRatio.^2)./(scosSPRatio.^2)*nPolarizations;
scosBeta=beta0/scosNModes;

if laserMaxPower/ansi_reference<=nSources
    nSources=laserMaxPower/ansi_reference;
elseif nSources<1
    nSources=1;
end

fprintf('Section 2 complete.\n');
%% 3. Estimate the count rate at each of the simulated source-detector separations

% Compute the intensity distribution based on the diffuse reflectance from
% the simulation
[intDist,intScale]=compute_intensity_scale(data_directory,sds_reference,cr_reference);

% Scale the DCS and SCOS count rate based upon the number of polarization
% states, and the s/p ratio/number of dcs modes
dcsCR=intScale*dcsQE*dcsNModes;
scosCR=intScale*scosQE*nPolarizations/scosSPRatio.^2;

fprintf('Section 3 complete: Measurement count rates estimated.\n');
%% 4. Generate the correlation functions from the photon history files

% Baseline blood flow
[g2_baseline,g1_baseline,tau]=calculate_g2_g1(mch_history_file,...
 'Db',tissueBFi,'tau',tau,'lambda',wv_reference*1e-6,'max_photons',1e9,'beta',1);

% Perturbed blood flow
[g2_perturb,g1_perturb,tau]=calculate_g2_g1(mch_history_file,...
 'Db',tissueBFi.*sensitivityPerturbation,'tau',tau,'lambda',wv_reference*1e-6,'max_photons',1e9,'beta',1);

fprintf('Section 4 complete: DCS correlation functions calculated.\n');
%% 5. Generate the speckle contrast for each measurement

% Baseline squared fundamental speckle contrast
kf2_baseline=compute_speckle_contrast(g1_baseline,tau,scosExpT);

% Perturbed squared fundamental speckle contrast
kf2_perturb=compute_speckle_contrast(g1_perturb,tau,scosExpT);

fprintf('Section 5 complete: SCOS kf2 calculated.\n');
%% 6. Compute the noise of the g2 curve

% Compute the noise for the baseline and the perturbed case for the CW
% laser
g2_baseline_noise_cw=generate_dcs_g2_noise(g2_baseline,dcsBeta,dcsCR*nSources,tau,1/fs_bfi);
g2_perturb_noise_cw=generate_dcs_g2_noise(g2_perturb,dcsBeta,dcsCR*nSources,tau,1/fs_bfi);

% Compute the noise for the baseline and the perturbed case for the pulsed
% laser case
g2_baseline_noise_pwm=generate_dcs_g2_noise(g2_baseline,dcsBeta,dcsCR*laserMaxPower/ansi_reference,tau,min([1,(ansi_reference*nSources)/laserMaxPower])/fs_bfi);
g2_perturb_noise_pwm=generate_dcs_g2_noise(g2_perturb,dcsBeta,dcsCR*laserMaxPower/ansi_reference,tau,min([1,(ansi_reference*nSources)/laserMaxPower])/fs_bfi);

% The percent difference in the noise between the two cases is quite small,
% and so their average will be used for noise estimation
g2_noise_cw=g2_baseline_noise_cw/2+g2_perturb_noise_cw/2;
g2_noise_pwm=g2_baseline_noise_pwm/2+g2_perturb_noise_pwm/2;

fprintf('Section 6 complete: Noise of the DCS correlation functions calculated.\n');
%% 7. Compute the noise of the squared speckle contrast 

% Estimate the frame rate and the power input of each strategy for SCOS
% Each array has 3 entries, the first is CW illumination without a pwm
% strategy, the second is a pwm strategy where the power input is limited,
% the third is a pwm strategy  where the frame rate is limited.
powerIn=[ansi_reference*nSources,min([laserMaxPower,ansi_reference*nSources/min([scosMaxFrameRate,1./scosExpT])/scosExpT]),laserMaxPower];
frameRate=[min([scosMaxFrameRate,1./scosExpT]),min([scosMaxFrameRate,1./scosExpT]),min([scosMaxFrameRate,ansi_reference*nSources/scosExpT/laserMaxPower])];

scosCRScale=powerIn./ansi_reference;
% Compute the noise for the baseline and the perturbed case for the CW
% illumination case
kf2_baseline_noise_cw=generate_scos_kf2_noise(kf2_baseline,scosBeta,scosCR*scosCRScale(1),scosExpT,scosSPRatio,scosReadNoise,scosPixelNumber);
kf2_perturb_noise_cw=generate_scos_kf2_noise(kf2_perturb,scosBeta,scosCR*scosCRScale(1),scosExpT,scosSPRatio,scosReadNoise,scosPixelNumber);

% Compute the noise for the baseline and the perturbed case for the pwm
% strategies
% First is the power input modulating strategy
kf2_baseline_noise_pwm1=generate_scos_kf2_noise(kf2_baseline,scosBeta,scosCR*scosCRScale(2),scosExpT,scosSPRatio,scosReadNoise,scosPixelNumber);
kf2_perturb_noise_pwm1=generate_scos_kf2_noise(kf2_perturb,scosBeta,scosCR*scosCRScale(2),scosExpT,scosSPRatio,scosReadNoise,scosPixelNumber);

% Second is the frame rate modulating strategy
kf2_baseline_noise_pwm2=generate_scos_kf2_noise(kf2_baseline,scosBeta,scosCR*scosCRScale(3),scosExpT,scosSPRatio,scosReadNoise,scosPixelNumber);
kf2_perturb_noise_pwm2=generate_scos_kf2_noise(kf2_perturb,scosBeta,scosCR*scosCRScale(3),scosExpT,scosSPRatio,scosReadNoise,scosPixelNumber);

% The percent difference in the noise between the two cases is quite small,
% and so their average will be used for noise estimation
kf2_noise_cw=kf2_baseline_noise_cw/2+kf2_perturb_noise_cw/2;
kf2_noise_pwm1=kf2_baseline_noise_pwm1/2+kf2_perturb_noise_pwm1/2;
kf2_noise_pwm2=kf2_baseline_noise_pwm2/2+kf2_perturb_noise_pwm2/2;

fprintf('Section 7 complete: SCOS kf2 noise estimates calculated.\n');
%% 8. Fit the clean autocorrelation and speckle contrast for the sensitivity estimate
%  All fitting is performed with a fixed beta, though the function for
%  fitting beta is included

% Initialize the bfi fitting arrays
dcs_bfi_baseline=zeros(size(intScale));
dcs_bfi_perturb=zeros(size(intScale));

scos_bfi_baseline=zeros(size(intScale));
scos_bfi_perturb=zeros(size(intScale));

% Begin fitting
for sdsIter=1:length(intScale)
   % Fit the DCS baseline results   
   dcs_bfi_baseline(sdsIter)=g2_fit_fixBeta(tau,g2_baseline(:,sdsIter)',...
       intDist(sdsIter),1,'mua',fit_mua,'msp',fit_msp,'n',fit_n,...
       'wv',wv_reference*1e-6,'g2_lb',fit_range(2),'g2_ub',fit_range(1));
   
   % Fit the DCS perturbed results   
   dcs_bfi_perturb(sdsIter)=g2_fit_fixBeta(tau,g2_perturb(:,sdsIter)',...
       intDist(sdsIter),1,'mua',fit_mua,'msp',fit_msp,'n',fit_n,...
       'wv',wv_reference*1e-6,'g2_lb',fit_range(2),'g2_ub',fit_range(1));
   
   % Fit the SCOS baseline data
   scos_bfi_baseline(sdsIter)=kf2_fit_fixBeta(scosExpT,kf2_baseline(sdsIter),...
       intDist(sdsIter),1,'mua',fit_mua,'msp',fit_msp,'n',fit_n,...
       'wv',wv_reference*1e-6);
   
   % Fit the SCOS perturbation data
   scos_bfi_perturb(sdsIter)=kf2_fit_fixBeta(scosExpT,kf2_perturb(sdsIter),...
       intDist(sdsIter),1,'mua',fit_mua,'msp',fit_msp,'n',fit_n,...
       'wv',wv_reference*1e-6);
end

% Estimate the sensitivity of DCS and SCOS to the cerebral blood flow
% change
dcs_sensitivity=(dcs_bfi_perturb-dcs_bfi_baseline)./dcs_bfi_baseline/brainPerturb;
scos_sensitivity=(scos_bfi_perturb-scos_bfi_baseline)./scos_bfi_baseline/brainPerturb;

fprintf('Section 8 complete: Clean SCOS and DCS data fit, and estimates of cerebral sensitivity calculated.\n');
%% 9. Fit noisy observations of the autocorrelation function and kf2 to estimate the 
%     noise properties of the signal.

% Estimate the DCS bfi noise for the cw power delivery strategy
dcs_bfi_std_cw=compute_dcs_bfi_noise(tau,g2_baseline,g2_noise_cw,intDist,dcsBeta,...
    'mua',fit_mua,'msp',fit_msp,'n',fit_n,'wv',wv_reference*1e-6,'g2_lb',fit_range(2),'g2_ub',fit_range(1));
% Estimate the bfi noise for the pwm power delivery strategy                                        
dcs_bfi_std_pwm=compute_dcs_bfi_noise(tau,g2_baseline,g2_noise_pwm,intDist,dcsBeta,...
    'mua',fit_mua,'msp',fit_msp,'n',fit_n,'wv',wv_reference*1e-6,'g2_lb',fit_range(2),'g2_ub',fit_range(1));

% Estimate the SCOS bfi noise for the cw power delivery strategy
scos_bfi_std_cw=compute_scos_bfi_noise(scosExpT,kf2_baseline,kf2_noise_cw,intDist,scosBeta,...
    'mua',fit_mua,'msp',fit_msp,'n',fit_n,'wv',wv_reference*1e-6);
% Estimate the SCOS bfi noise for the first pwm power delivery strategy
scos_bfi_std_pwm1=compute_scos_bfi_noise(scosExpT,kf2_baseline,kf2_noise_pwm1,intDist,scosBeta,...
    'mua',fit_mua,'msp',fit_msp,'n',fit_n,'wv',wv_reference*1e-6);
% Estimate the SCOS bfi noise for the second pwm power delivery strategy
scos_bfi_std_pwm2=compute_scos_bfi_noise(scosExpT,kf2_baseline,kf2_noise_pwm2,intDist,scosBeta,...
    'mua',fit_mua,'msp',fit_msp,'n',fit_n,'wv',wv_reference*1e-6);

% Calculation of the coefficient of variation for the measurements, and
% accounting for frame averaging (SCOS)
dcs_cov_cw=mean([dcs_bfi_std_cw./dcs_bfi_baseline;dcs_bfi_std_cw./dcs_bfi_perturb],1);
dcs_cov_pwm=mean([dcs_bfi_std_pwm./dcs_bfi_baseline;dcs_bfi_std_pwm./dcs_bfi_perturb],1);

scos_cov_cw=mean([scos_bfi_std_cw./scos_bfi_baseline;scos_bfi_std_cw./scos_bfi_perturb],1)./sqrt(frameRate(1)/fs_bfi);
scos_cov_pwm1=mean([scos_bfi_std_pwm1./scos_bfi_baseline;scos_bfi_std_pwm1./scos_bfi_perturb],1)./sqrt(frameRate(2)/fs_bfi);
scos_cov_pwm2=mean([scos_bfi_std_pwm2./scos_bfi_baseline;scos_bfi_std_pwm2./scos_bfi_perturb],1)./sqrt(frameRate(3)/fs_bfi);

fprintf('Section 9 complete: Noisy SCOS and DCS data fit, and estimates of coefficient of variation calculated.\n');
%% 10. Compute the contrast-to-noise ratio

% Calculation of dcs cnr
dcs_cnr_cw=dcs_sensitivity./dcs_cov_cw;
dcs_cnr_pwm=dcs_sensitivity./dcs_cov_pwm;

% Calcultion of scos cnr 
scos_cnr_cw=scos_sensitivity./scos_cov_cw;
scos_cnr_pwm1=scos_sensitivity./scos_cov_pwm1;
scos_cnr_pwm2=scos_sensitivity./scos_cov_pwm2;

fprintf('Section 10 complete: SCOS and DCS contrast-to-noise ratio calculated.\n');
%% 11. Plot the results
plotResults(dcs_sensitivity,scos_sensitivity,dcs_cov_cw,dcs_cov_pwm,scos_cov_cw,scos_cov_pwm1,scos_cov_pwm2,...
dcs_cnr_cw,dcs_cnr_pwm,scos_cnr_cw,scos_cnr_pwm1,scos_cnr_pwm2,intDist,fs_bfi)

fprintf('Section 11 complete: Results plotted.\n');
%% 12. Clean up path dependencies
rmpath(genpath(mmcPath));
rmpath(genpath(scatterbrainsPath));
rmpath(genpath(iso2meshPath));