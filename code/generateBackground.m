function bg_data = generateBackground(bg_density,membrane_width)
% This function provides n_bg in terms of number of photons per um2.
%
% Inputs:
%   bg_density: density of background photons per um2
%   membrane_width: width of the membrane in nm
%
% Outputs:
%   bg_data: 2-column matrix with the x and y coordinates of the background
%   photons
% 
%
% Created by Weihong Yeo, Northwestern University, 2022-07-29.
% Last modified by Weihong Yeo, Northwestern University, 2022-10-02.
% 
% ######################################################################### 
% Changelog
% ######################################################################### 
% 2022-10-02: corrected the n_bg formula
% 

% default membrane width is 1000 um
if nargin < 2
    membrane_width = 1000;
end

% generate background data
n_bg = round(bg_density * (membrane_width / 1000)^2);
bg_data = rand(n_bg,2)*membrane_width;

end
