function [bfi]=g2_fit_fixBeta(tau,g2,sds,beta,options)
%% Fits for the BFi of a g2 curve with a supplied beta value.
% author: Mitchell Robinson (robinson.mitchell.b <at> gmail.com)

% Inputs
% tau: An array of autocorrelation lags in seconds
% g2: An array of values containing the autocorrelation function
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
% g2_lb: Lower cutoff value for the fitting region of the correlation
% function
% g2_ub: Upper cutoff value for the fitting region of the correlation
% function
% lsqCurveFitOptions: Options used for the fitting function

% Outputs
% bfi: Fit value for the flow measured from the g2 curve

%% Handling input arguments
arguments
tau double
g2 double
sds double
beta double
options.mua (1,1) double = 0.017
options.msp (1,1) double = 0.80
options.wv (1,1) double = 850e-6
options.n (1,1) double =1.4
options.bfi_lb (1,1) double = 0
options.bfi_ub (1,1) double = 20
options.g2_lb (1,1) double = 0.05
options.g2_ub (1,1) double = .999
options.lsqCurveFitOptions = optimoptions('lsqcurvefit','display','none','OptimalityTolerance',1e-12,...
        'FunctionTolerance',1e-12,'StepTolerance',1e-12);
end
x0=options.bfi_ub/2+options.bfi_lb/2;

if ~(length(g2)==length(tau))
    error('Tau array and g2 array are not the same length.');
end
%% Simplifying the options inputs
mua=options.mua;
msp=options.msp;
wv=options.wv;
n=options.n;
bfi_lb=options.bfi_lb;
bfi_ub=options.bfi_ub;
g2_lb=options.g2_lb;
g2_ub=options.g2_ub;
fit_options=options.lsqCurveFitOptions;

%% Find the fitting index for the correlation function based on the g2_lb and g2_ub values
if or(g2_ub==inf,g2_lb==-inf)
    fitInds=1:length(g2);
else
    [~,fitInd(1)]=min(abs(g2-1-g2_ub));
    [~,fitInd(2)]=min(abs(g2-1-g2_lb));
    fitInds=fitInd(1):fitInd(2);
end
tau=tau(fitInds);
g2=g2(fitInds);

%% Data are fit using the generated model
g2Model=make_g2_model([mua,msp,n,wv,sds]);
fitOut=lsqcurvefit(@(x,xdata)g2Model([x,beta],xdata),x0,tau,g2,bfi_lb,bfi_ub,fit_options);
bfi=fitOut(1)*1e-6; % mm^2/s
end