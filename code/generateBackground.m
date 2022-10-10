function bg_data = generateBackground(bg_density,membrane_width)
% This function provides n_bg in terms of number of photons per um2.
%
%
% Created by Weihong Yeo, Northwestern University, 2022-07-29.
% Last modified by Weihong Yeo, Northwestern University, 2022-10-02.
% 
% # Changelog
% ## 2022-10-02
% - corrected n_bg formula
% 

if nargin < 2
    membrane_width = 1000;
end

n_bg = round(bg_density * (membrane_width / 1000)^2);
bg_data = rand(n_bg,2)*membrane_width;

end
