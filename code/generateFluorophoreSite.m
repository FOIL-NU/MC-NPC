function [fluor_pos,segment_igg] = generateFluorophoreSite(N,npc_params)
% This function generates the fluorophore site relative to a label site.
% 
% N gives the number of igg to generate at the same time.
% npc_params will orient the antibodies in the correct way for simulation.
% igg_params are to be given as a structure array. read the code to see what changes.
% 
% fluor_pos gives a 3 by N matrix of fluorophore label locations.
% 
% Created by Weihong Yeo, Northwestern University, 2022-08-19.
% Last modified by Weihong Yeo, Northwestern University, 2022-10-02.
% 
% # Changelog
% ## 2022-08-19
% - added rotations for labeling on different locations
%
% ## 2022-10-02
% - changed codebase to a more realistic model found in
%   https://pubs.acs.org/doi/10.1021/acsnano.1c03677
%
% ## 2022-10-05
% - updated code to using random labeling on secondary antibodies, since we
%   do not label onto Glu in our antibodies.

%% Set default parameters

if (nargin == 0) || isempty(N)
    N = 1;
end

if (nargin <= 1) || isempty(npc_params.location_index)
    npc_params.location_index = 1; % outer_membrane
end

npc_params.location_index = 1;
min_dist_tol = 3;

P = [0.31, 0.31, 0.38];

% primary antibody values
L_igg1 = {[6.527, 6.669, 6.819, 7.157]; [2.520, 7.239]};
phi_igg1 = {[-22.5, 22.5]; [15.0, 127.6]};
theta_igg1 = {[-22.5, 22.5]; [0, 360]};

% secondary antibody values
% values if AF647 is on Fab1
L_Fab1 = [L_igg1; {[0.0, 0.7]; 2.3; 0}];
phi_Fab1 = [phi_igg1; {[67.5, 112.5]; 90; 0}];
theta_Fab1 = [theta_igg1; {[-22.5, 22.5]; [0, 360]; 0}];
     
 % values if AF647 is on Fab2
L_Fab2 = [L_igg1; {[6.527, 6.669, 6.819, 7.157]; [0.0, 0.7]; 2.3}];
phi_Fab2 = [phi_igg1; {[67.5, 112.5]; 'phi_4'; 90}];
theta_Fab2 = [theta_igg1; {[-22.5, 22.5]; [0, 360]; [0, 360]}];

 % values if AF647 is on Fc
L_Fc = [L_igg1; {[6.527, 6.669, 6.819, 7.157]; [0.0, 8.5]; 2.1}];
phi_Fc = [phi_igg1; {[67.5, 112.5]; [15.0, 127.6]; 90}];
theta_Fc = [theta_igg1; {[-22.5, 22.5]; [0, 360]; [0, 360]}];

L_igg = {L_Fab1, L_Fab2, L_Fc};
phi_igg = {phi_Fab1, phi_Fab2, phi_Fc};
theta_igg = {theta_Fab1, theta_Fab2, theta_Fc};

segment_igg = randsample(3,N,true,P);

xyz_full = zeros(N,3,6);

for i_igg = 1:3
    N_segment = sum(segment_igg == i_igg);  % find the number of fluorophores that lie on each segment
    if N_segment > 0
        % get the relevant angle ranges.
        L_segment = L_igg{i_igg};
        phi_segment = phi_igg{i_igg};
        theta_segment = theta_igg{i_igg};
        
        xyz = zeros(N_segment,3,6);
        
        min_dist = zeros(N_segment,1);
        steric_npc = false(N_segment,1);

        while any((min_dist < min_dist_tol) | steric_npc)
            update_index = (min_dist < min_dist_tol) | steric_npc;
            
            N_segment_ = sum(update_index);
            
            xyz_ = zeros(N_segment_,3,6);
            
            for i_segment = length(L_segment):-1:1
                % get random angles
                L_ = get_rv(L_segment{i_segment},N_segment_);
                phi_ = get_rv(phi_segment{i_segment},N_segment_);
                theta_ = get_rv(theta_segment{i_segment},N_segment_);

                xyz_temp = get_xyz(L_,phi_,theta_);

                if i_segment < length(L_segment) % skip first rotate joint because the value was generated relative to that
                    xyz_(:,:,i_segment+1:end) = rotate_joint(phi_,theta_,xyz_(:,:,i_segment+1:end)) + xyz_temp;
                else
                    xyz_(:,:,i_segment+1) = xyz_temp;
                end
            end
            
            xyz(update_index,:,:) = xyz_;

            % compute pairwise distances
            min_dist = min_dist_tol * ones(N_segment,1);
            steric_npc = false(N_segment,1);
            
            min_dist_ = zeros(N_segment_,1);
            steric_npc_ = false(N_segment_,1);
            for ii = 1:N_segment_
                if i_igg == 1
                    distances = triu(squareform(pdist(squeeze(xyz_(ii,:,1:end-2))')),2);
                else
                    distances = triu(squareform(pdist(squeeze(xyz_(ii,:,1:end-1))')),2);
                end
                min_dist_(ii) = min(distances(distances ~= 0)); % get the minimum non-zero distance
                steric_npc_(ii) = any(xyz_(ii,2,:) < 0); % check if any y < 0
            end
            min_dist(update_index) = min_dist_;
            steric_npc(update_index) = steric_npc_;
        end
        
        assert(any(any(xyz(:,2,:) < 0)) == 0);
        
        xyz_full(segment_igg == i_igg,:,:) = xyz;
    end
end

fluor_pos = squeeze(xyz_full(:,:,end))';

%% orientate the primary and secondary antibodies to based on location of label
% we do not need to do anything if location is 'outer' or 'nuclear'.

if ~any(npc_params.location_index == [2,3])
    fluor_pos = rotx(90) * fluor_pos;
end

end

%% helper functions

function xyz = rotate_joint(phi,theta,xyz)
assert(size(xyz,2) == 3);

for i1 = 1:size(xyz,1)
    for i2 = 1:size(xyz,3)
        xyz(i1,:,i2) = (roty(theta(i1)) * rotx(phi(i1)) * xyz(i1,:,i2)')';
    end
end
end

function xyz = get_xyz(L,phi,theta)

x = sind(phi) .* sind(theta);
y = cosd(phi);
z = sind(phi) .* cosd(theta);

xyz = L.*[x,y,z];

end

function rv = get_rv(val,n)
if nargin < 2
    n = 1;
end

if strcmpi(val, 'phi_4')
    rv = get_xi(n);
elseif length(val) == 1
    rv = repmat(val,n,1);
elseif length(val) == 2
    rv = rand(n,1) .* (val(2) - val(1)) + val(1);
else
    rv = datasample(val,n)';
end

end

% pdf of this function based on Bongini et. al., Proc. Natl. Acad. Sci. U.S.A., 2004, Appendix 2.
function xi = get_xi(n)

sigma_theta = 91.6;
theta0 = 180;
phi_max = 120;  % corrected this to 120 because values beyond 120 always give complex numbers.
phi_min = 15;

xi = ones(n,1).* 1i;

while any(imag(xi) ~= 0)
    n = sum(imag(xi) ~= 0);
    
    phi1 = rand(n,1) * (phi_max - phi_min) + phi_min;
    phi2 = rand(n,1) * (phi_max - phi_min) + phi_min;
    theta = theta0 - sigma_theta .* tand((1 - 2.*rand(n,1)).*atand(theta0./sigma_theta));
    
    cos_phi1 = cosd(phi1);
    cos_phi2 = cosd(phi2);
    cos2_phi1 = cos_phi1 .* cos_phi1;
    cos2_phi2 = cos_phi2 .* cos_phi2;
    
    xi_temp = acosd(sqrt(1+cos_phi1-2.*cos2_phi1).*sqrt(1+cos_phi2-2.*cos2_phi2).*cosd(theta)+cos_phi1.*cos_phi2);
    
    phi1(imag(xi) ~= 0) = phi1;
    phi2(imag(xi) ~= 0) = phi2;
    theta(imag(xi) ~= 0) = theta;
    xi(imag(xi) ~= 0) = xi_temp;
end

end
