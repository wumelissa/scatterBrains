
function y = semi_infinite_g2(x,taus,fit_options)

% returns autocorrelation using semi-infinite diffuse correlation equation

% input:
%   x: array with [beta, Db] value
%   taus: array of taus
%   fit_options: structure with fields
%       mu_a: absorption, mm-1
%       mu_s: scattering, mm-1
%       alpha: probability of scattering from a moving scatterer
%       lambda: wavelength, mm-1
%       rho: source-detector distance, mm
%       n: index of refraction

% author: Melissa Wu (wu.melissa.m <at> gmail.com)
% contributing author: Stefan Carp (stefan.carp <at> mgh.harvard.edu)

% this file is part of scatterBrains
% License: GPLv3

%% variable assignment

mu_a = fit_options.mu_a;
mu_s = fit_options.mu_s;
alpha = fit_options.alpha;
lambda = fit_options.lambda;
rho = fit_options.rho;
n0=fit_options.n;
beta = x(1);
D = x(2)/1e9;

%%

%equations

R = -1.440./n0^2+0.710/n0+0.668+0.0636*n0; % Effective reflection coefficient
ze = 2/3*(1+R)/(1-R); 
z_b = ze/mu_s;
z_0 = (mu_a + mu_s)^-1;
r_1 =(rho^2 + z_0^2)^0.5;
r_2 = (rho^2 + (z_0 + 2*z_b)^2)^0.5;
k_0 = 2*pi*fit_options.n/lambda;
delta_r2 = 6*D*taus;
K = (3*mu_a*mu_s + alpha*mu_s^2*k_0^2*delta_r2).^(0.5);
K_tau0 = (3*mu_a*mu_s)^0.5;
G1 = (3*mu_s/(4*pi))*(exp(-K*r_1)/r_1 - exp(-K*r_2)/r_2);
G1_tau0=(3*mu_s/(4*pi))*(exp(-K_tau0*r_1)/r_1 - exp(-K_tau0*r_2)/r_2);
small_g1 = abs(G1./G1_tau0);
y = (1+beta*(small_g1.^2))'; %final g2 equation

end

