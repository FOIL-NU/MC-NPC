function bg_data = generateBackground(n_bg,membrane_width)
% provide n_bg in terms of number of photons per um2.
% This function generates the background based on a given background density.
%
%
% Created by Weihong Yeo, Northwestern University, 2022-08-19.
%

if nargin < 1
    ng_density = 40;
end

if nargin < 2
    membrane_width = 1000;
end

n_bg = round(bg_density * (membrane_width / 1000));
bg_data = rand(n_bg,2)*membrane_width;

end