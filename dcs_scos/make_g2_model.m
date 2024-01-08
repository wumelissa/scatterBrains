function g2Model=make_g2_model(optical_properties)
%% Builds the g2 model based on the semi-infinite solution to the correlation 
%  diffusion equation for fitting for BFi.

% author: Mitchell Robinson (robinson.mitchell.b <at> gmail.com)

% Optical properties are an array where
% Element 1: mua
% Element 2: msp
% Element 3: n
% Element 4: wv
% Element 5: sds

% Constants for the model fit are calculated based on the inputs
k0=2*pi/optical_properties(4)*optical_properties(3);
Reff=-1.44*optical_properties(3)^-2+0.71*optical_properties(3)^-1+0.668+0.0636*optical_properties(3);
ltr=1/optical_properties(2);
zb=2*ltr*(1+Reff)/(3*(1-Reff));
r1=sqrt(ltr^2+optical_properties(5)^2);
rb=sqrt((2*zb+ltr)^2+optical_properties(5)^2);
scaleVal=exp(-sqrt(3*optical_properties(1)*optical_properties(2))*r1)/r1-...
    exp(-sqrt(3*optical_properties(1)*optical_properties(2))*rb)/rb;

% Calculation of more coefficients to reduce the length of the model
% function
cfs=[3*optical_properties(1)*optical_properties(2),...
    6*k0^2*optical_properties(2)^2*1e-6,...
    r1,...
    rb,...
    scaleVal];

g2Model=@(x,xdata)1+x(2)*(exp(-sqrt(cfs(1)+cfs(2)*x(1)*xdata)*cfs(3))/(cfs(3)*cfs(5))-...
    exp(-sqrt(cfs(1)+cfs(2)*x(1)*xdata)*cfs(4))/(cfs(4)*cfs(5))).^2;
end