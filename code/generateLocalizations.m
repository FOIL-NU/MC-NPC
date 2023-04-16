function localizations = generateLocalizations(n_blinks,loc_dist)
% This function generates localizations based on the expected number of
% blinks given by n_blinks.
% 
% Inputs:
%   n_blinks: The expected number of blinks. If n_blinks is a scalar, the
%       function will use the value as the expected number of blinks. If
%       n_blinks is the nup_param structure, the function will sample from
%       the distribution of the number of blinks. Otherwise, the function
%       will use the default value of 12.
%   loc_dist: A distribution object or a scalar value. If a distribution
%       object is provided, the function will sample from the distribution
%       to generate localization uncertainty. If a scalar value is provided,
%       the function will use the value as a constant value for
%       localization uncertainty. If loc_dist is not provided, the function
%       will use the default lognormal distribution.
%
% Outputs:
%   localizations: A 3-by-n_blinks matrix containing the localizations of
%       the fluorophores attached on the simulated NUPs. The rows contains
%       the x, y, and z coordinates of the localizations, respectively.
%
%
% Created by Weihong Yeo, Northwestern University, 2022-08-19.
% Last modified by Weihong Yeo, Northwestern University, 2023-04-15.
%
% ######################################################################### 
% Changelog
% ######################################################################### 
% 2023-04-15: Added the option to use a scalar value or leave blank for the
% input of the number of blinks.
%
% 2023-04-09: Added the option to use a distribution object as the input for
% localization uncertainty.
%
% 2022-08-19: Created the function.
%

% Check inputs
if nargin < 1 || isempty(n_blinks)
    n_blinks = 12;
elseif isstruct(n_blinks)
    % Calculate the expected number of blinks from the nbinrnd distribution
    r = round((n_blinks.avg_sml * exp(-0.5)) / (1 - exp(-0.5)));
    n_blinks = round(nbinrnd(r, exp(-0.5)));
end

if nargin < 2 || isempty(loc_dist)
    % If loc_dist is empty or not provided, use the default lognormal distribution
    loc_uncertainty = random('lognormal', 2.29608, 0.401825, 1, n_blinks);
elseif isscalar(loc_dist)
    % If loc_dist is a scalar, use it as a constant value for loc_uncertainty
    loc_uncertainty = repmat(loc_dist, 1, n_blinks);
else
    % Sample from the input distribution object
    loc_uncertainty = random(loc_dist, 1, n_blinks);
end

localizations = [randn(2, n_blinks) .* loc_uncertainty; zeros(1, n_blinks)];

end
