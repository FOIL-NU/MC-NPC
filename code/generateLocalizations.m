function localizations = generateLocalizations(nup_param,loc_uncertainty)
% This function generates localizations based on the expected number of
% blinks given by nup_params.
% 
% We assume a lognormally distributed localization uncertainty of 2.29608,0.401825 nm by default.
%
%
% Created by Weihong Yeo, Northwestern University, 2022-08-19.
%

%% check inputs
if nargin < 1
    n_blinks = 12;
else
    r = round((nup_param.avg_sml * exp(-0.5)) / (1 - exp(-0.5)));
    n_blinks = round(nbinrnd(r,exp(-0.5)));
end
if nargin < 2
    loc_uncertainty = repmat(random('lognormal',2.29608,0.401825,1,n_blinks),2,1);
end

localizations = [randn(2,n_blinks).*loc_uncertainty; zeros(1,n_blinks)];

end