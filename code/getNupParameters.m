function nup_parameters = getNupParameters(nup_species)
% This function gives the parameters of the NUP protein in the NPC used in
% simulations. 
%
% It provides location information, radial distances for each labeling
% site, probability of labeling onto each site, and average single molecule
% localization.
%
% Input:
%   nup_species: string, name of the NUP protein.
%
% Output:
%   nup_parameters: struct, containing the parameters of the NUP protein.
%       nup_parameters.location_name: string, name of the location of the
%           NUP protein.
%       nup_parameters.minor_angle: scalar, minor angle of the site
%           relative to the inner torus curvature of the NPC.
%       nup_parameters.radial_dist: 1D array, radial distances of each
%           labeling site.
%       nup_parameters.prob_label: 1D array, probability of labeling onto
%           each site.
%       nup_parameters.avg_sml: scalar, average single molecule
%           localization.
%
% 
% Created by Weihong Yeo, Northwestern University, 2022-07-29.
% Last modified by Weihong Yeo, Northwestern University, 2023-04-15.
% 
% ######################################################################### 
% Changelog
% ######################################################################### 
% 2023-04-15: added the minor angle parameter to the csv file.
%
% 2023-04-09: changed code to read a csv file instead of hardcoding the
% parameters.
%
% 2022-08-19: included more NUPs in the database, and added default values
% for prob_label and avg_sml.
% 
% 2022-07-29: created.
%

% Load data from CSV file
data = readtable('nup_parameters.csv');

% Find the row corresponding to the input nup_species
nup_idx = find(strcmpi(data.Nup_species, nup_species), 1);
assert(~isempty(nup_idx), 'Invalid Nup species: %s', nup_species);

% Extract parameters from row
nup_parameters.location_name = data.Location_name{nup_idx};
nup_parameters.minor_angle = data.Minor_angle(nup_idx);
nup_parameters.radial_dist = data.Radial_distance(nup_idx);
nup_parameters.prob_label = data.Probability_label(nup_idx);
nup_parameters.avg_sml = data.Average_SML(nup_idx);

end