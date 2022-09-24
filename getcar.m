function [AM,SR,car] = getcar(CaCO3,dZ,age_top,age_bottom)
 
% The array AM contains the age model for a segment long length(CaCO3)
% SR contains the sedimentation rate for depths in between the CaCO3
% measurements. Therefore length(SR) = length(CaCO3)-1.
% car contains a single value and represents the pelagic sedimentation
% rate. The sedimentation rate (SR) at a given depth in the core is
% calculated as SR = CAR/(CaCO3/100)
% Furthermore psr is obtained from the input data for which the CaCO3 data,
% the sample spacing (dZ)and the total amount of time 
% (age_bottom - age_top) represented by the segment are needed. The most
% crucial formula in this function is: car =sum(midpointCaCO3values)*(dZ/T)
% which determines the car value = a constant, for the segment considered.
% The function getcar requires input of equal spaced data, dZ = constant.
 
% =============================================================
% Dr. Frank Peeters, Faculty of Earth- and Life Sciences, Vrije
% Universiteit Amsterdam, De Boelelaan 1085, 1081 HV Amsterdam,
% The Netherlands. email: f.j.c.peeters@vu.nl
% =============================================================
 
if nargin < 4
    error('Not enough input parameters')
end
 
CaCO3 = CaCO3/100; % convert to fraction
ni = length(CaCO3);
T = age_bottom - age_top;
midpointCaCO3values = zeros(ni-1,1);
AM = zeros(length(CaCO3),1);
 
for i = 1:ni-1
    midpointCaCO3values(i) = (CaCO3(i)+CaCO3(i+1))/2;
end
 
car =sum(midpointCaCO3values)*(dZ/T);
SR = car./midpointCaCO3values;
dt = dZ./SR;
AM(1) = age_top;
AM(2:end) = cumsum(dt)+age_top;
end
