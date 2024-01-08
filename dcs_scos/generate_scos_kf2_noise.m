function kf2Noise=generate_scos_kf2_noise(kf2,beta,CR,expT,spRatio,sigma_r,nPixels)
%% Computes the standard deviation of the estimated kf2 after corrections to
%  variance are performed. Based on the model presented in 
% "Zilpelwar, S., et al.(2022). A model of dynamic speckle evolution for evaluating laser speckle contrast measurements of tissue dynamics. Biomedical Optics Express, 13(12), 6533–6549."

% author: Mitchell Robinson (robinson.mitchell.b <at> gmail.com)

% Inputs
% kf2: The fundamental speckle contrast not yet scaled by a beta value 
% beta: Coherence factor of the measurment, based on coherence of the
% laser, polarization state, and s/p ratio.
% CR: Count rate of the measurement for each of the k2f's 
% expT: Exposure time used
% spRatio: Value quantifying the ratio of Speckle diameter / pixel side
% length
% sigma_r: Temporal read noise
% nPixels: Number of pixels in the camera sensor

% Outputs
% k2fNoise: standard deviation of the squared speckle contrast

%% Argument handling
arguments
        kf2 double
        beta double
        CR double
        expT (1,1) double
        spRatio (1,1) double
        sigma_r (1,1) double
        nPixels (1,1) double
end
if ~(length(kf2)==length(CR))
    error('Each squared speckle contrast requires a photon count rate for noise calculation.')
end
%% Noise calculation
% Constants
c_Ks=1.9; % Constant relating the measured shot noise variance to the standard deviation of the shot noise variance
c_Kr=1.47; % Constant relating the measured read noise variance to the standard deviation of the read noise variance
c_sp=0.292; % Dimensional factor that relates the s/p ratio to the spatial extent of the speckle on the pixels

% Adjustment of the k2f to account for the beta of the SCOS measurement
kf2=kf2*beta;

% Computing the scaling factor for converting the number of speckles
% measured to the number of independent speckle observations
f_sp=1./(1+c_sp*spRatio^2);

% Calculation of the number of independent speckle observations
NIO=nPixels*f_sp;

% Calculation of the shot and read noise variance contributions based on
% the camera (read noise) and the singal intensity (shot noise)
ks2=1./(CR*expT);
kr2=(sigma_r./(CR*expT)).^2;

sigma_kf2=kf2.*sqrt(2.*(1+2.*kf2)./NIO);
sigma_ks2=c_Ks.*ks2.*sqrt(1./nPixels);
sigma_kr2=c_Kr.*kr2.*sqrt(1./nPixels);

kf2Noise=sqrt(sigma_kf2.^2+sigma_ks2.^2+sigma_kr2.^2);
end