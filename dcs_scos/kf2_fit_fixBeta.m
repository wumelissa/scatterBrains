function [bfi]=kf2_fit_fixBeta(expT,kf2,sds,beta,options)
%% Fits for the BFi of a kf2 value with a supplied beta value.
% author: Mitchell Robinson (robinson.mitchell.b <at> gmail.com)

% Inputs
% expT: The exposure time used to generate the simulated contrast value
% kf2: The fundamental squared speckle contrast of the measurement
% mua: Mua used for fitting in the semi-infinite model in units of
% mm^-1
% msp: Musp used for fitting in the semi-infinite model in units of
% mm^-1
% wv: Wavelength used for fitting in units of mm
% sds: source-detector separation used for fitting in units of mm
% n: Index of refraction used for fitting
% bfi_lb: lower bounds for the fitting for parameters [bfi] in units of
% [mm^2/s]*1e6
% bfi_ub: upper bounds for the fitting for parameters [bfi] in units of
% [mm^2/s]*1e6
% lsqCurveFitOptions: Options used for the fitting function

% Outputs
% bfi: Fit value for the flow measured from the kf2 value

%% Handling input arguments
arguments
expT double
kf2 double
sds double
beta double
options.mua (1,1) double = 0.017
options.msp (1,1) double = 0.80
options.wv (1,1) double = 850e-6
options.n (1,1) double =1.4
options.bfi_lb (1,1) double = 0
options.bfi_ub (1,1) double = 20
options.lsqCurveFitOptions = optimoptions('lsqcurvefit','display','none','OptimalityTolerance',1e-12,...
        'FunctionTolerance',1e-12,'StepTolerance',1e-12);
end
x0=options.bfi_ub/2+options.bfi_lb/2;


%% Simplifying the options inputs
mua=options.mua;
msp=options.msp;
wv=options.wv;
n=options.n;
bfi_lb=options.bfi_lb;
bfi_ub=options.bfi_ub;
fit_options=options.lsqCurveFitOptions;


%% Data are fit using the generated model
kf2Model=make_kf2_model([mua,msp,n,wv,sds]);
fitOut=lsqcurvefit(@(x,xdata)kf2Model([x,beta],xdata),x0,expT,kf2,bfi_lb,bfi_ub,fit_options);

bfi=fitOut(1)*1e-6; % mm^2/s
end