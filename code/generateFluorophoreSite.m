function fluor_pos = generateFluorophoreSite(N,npc_params,igg_params)
% This function generates the fluorophore site relative to a label site.
% 
% N gives the number of igg to generate at the same time.
% npc_params will orient the antibodies in the correct way for simulation.
% igg_params are to be given as a structure array. read the code to see what changes.
% 
% fluor_pos gives a 3 by N matrix of fluorophore label locations.
% 
%
% Created by Weihong Yeo, Northwestern University, 2022-08-19.
% Last modified by Weihong Yeo, Northwestern University, 2022-08-19.
% 
% # Changelog
% ## 2022-08-19
% - added rotations for labeling on different locations

%% Set default parameters
% The default IgG sites are 2 half spheres with the northern hemisphere oriented 
% along the z direction.

igg1_el_lim = deg2rad(90);   % degrees of freedom in elevation for the primary antibody
igg1_az_lim = deg2rad(360);  % degrees of freedom in azimuth for the primary antibody
igg1_len = 15;

igg2_el_lim = deg2rad(90);   % degrees of freedom in elevation for the secondary antibody
igg2_az_lim = deg2rad(360);  % degrees of freedom in azimuth for the secondary antibody
igg2_len_lim = 15;

if (nargin == 0) || isempty(N)
    N = 1;
end

if (nargin <= 1) || isempty(npc_params.location_index)
    npc_params.location_index = "outer_membrane";
end
    
if nargin == 3
    if isfield(igg_params,'igg1')
        if isfield(igg_params.igg1,'el')
            igg1_el_lim = igg_params.igg1.el;
        elseif isfield(igg_params.igg1,'el_lim')
            igg1_el_lim = igg_params.igg1.el_lim;
        end
    
        if isfield(igg_params.igg1,'az')
            igg1_az_lim = igg_params.igg1.az;
        elseif isfield(igg_params.igg1,'az_lim')
            igg1_az_lim = igg_params.igg1.az_lim;
        end
    
        if isfield(igg_params.igg1,'len')
            igg1_len = igg_params.igg1.len;
        end
    end
    
    if isfield(igg_params,'igg2')
        if isfield(igg_params.igg2,'el')
            igg2_el_lim = igg_params.igg2.el;
        elseif isfield(igg_params.igg2,'el_lim')
            igg2_el_lim = igg_params.igg2.el_lim;
        end
    
        if isfield(igg_params.igg2,'len')
            igg2_len_lim = igg_params.igg2.len;
        elseif isfield(igg_params.igg2,'len_lim')
            igg2_len_lim = igg_params.igg2.len_lim;
        end
    end
end

%% 
% Generate the random fluorophore location

AZ = rand(N,2).*[igg1_az_lim,igg2_az_lim];
EL = asin(rand(N,2).*sin([igg1_el_lim,igg2_el_lim]));
R  = [repmat(igg1_len,N,1),rand(N,1).*igg2_len_lim];

[X,Y,Z] = sph2cart(AZ,EL,R);

X = sum(X,2);
Y = sum(Y,2);
Z = sum(Z,2);

fluor_pos = [X,Y,Z]';

%% orientate the primary and secondary antibodies to based on location of label
% we do not need to do anything if location is 'outer' or 'nuclear'.

if any(npc_params.location_index == [2,3])
    fluor_pos = roty(90) * fluor_pos;
end

end