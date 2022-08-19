function nup_parameters = getNupParameters(nup_species)
% This function gives the parameters of the NUP protein in the NPC used in
% simulations. 
%
% It provides location information, radial distances for each labeling
% site, probability of labeling onto each site, and average single molecule
% localization.
%
% Here, we assume that each nuclear pore complex can only be in one location,
% either outer, inner, transmembrane, or nuclear.
%
%
% Created by Weihong Yeo, Northwestern University, 2022-07-29.
% Last modified by Weihong Yeo, Northwestern University, 2022-08-19.
% 
% # Changelog
% ## 2022-08-19
% - included more NUPs in the database
% - added default values for prob_label and avg_sml.

%% Store a list of known NUPs at the associated binding sites
nup_index = ["fgrepeats","nup93","nup133","nup210"];
location_names = ["outer","inner","transmembrane","nuclear"];

%% Check input is within the list
assert(any(strcmpi(nup_index,nup_species)==1));

%% Simulation parameters
locations_index = [ % outer, inner, transmembrane, nuclear
      2, ... % fgrepeats
      2, ... % nup93
      1, ... % nup133
      3  ... % nup210
    ];

radial_distances = [ % nm, radial distance relative to centroid
      10.0, ... % fgrepeats
      40.0, ... % nup93
      53.5, ... % nup133
      80.0  ... % nup210
    ];

nup_parameters.location_index = locations_index(strcmpi(nup_index,nup_species));
nup_parameters.location_name = location_names(nup_parameters.location_index);
nup_parameters.radial_dist = radial_distances(strcmpi(nup_index,nup_species));
nup_parameters.prob_label = 0.6;
nup_parameters.avg_sml = 10;

end