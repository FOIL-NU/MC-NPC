function [fluor_pos, local_quaternions] = generateFluorophoreSites(site_xyz,site_quat)
% This function generates the fluorophore site relative to a label site
% using a forward-propagation method of the antibodies.
%
% Inputs:
%   site_xyz: a Nx3 matrix, where each row is the x, y, and z global
%       coordinates of the label site.
%   site_quat: a Nx1 cell array, where each row is the quaternion of the
%       label site.
%
% Outputs:
%   fluor_pos: a 3xN matrix, where each column is the position of the label
%       site relative to the NPC label site. The rows are the x, y, and z
%       global coordinates, respectively.
%   segment_igg: a 1xN vector, where each element is the segment of the IgG
%       that the fluorophore is attached to.
%
%
% Created by Weihong Yeo, Northwestern University, 2023-04-10.
% Last modified by Weihong Yeo, Northwestern University, 2023-04-10.
% 
% ######################################################################### 
% Changelog
% #########################################################################
% 2023-04-10: Created file.
%

% Set default parameters.
if (nargin < 1) || isempty(site_xyz)
    xyz = zeros(1,3);
    N = 1;
else
    N = size(site_xyz,1);
    xyz = nan(N,3,6);
    xyz(:,:,1) = site_xyz;
end

local_quaternions = cell(N,6);
iterations = 0;

if (nargin < 2) || isempty(site_quat)
    for ii = 1:N
        local_quaternions{ii,1} = quaternion(1,0,0,0);
    end
else
    assert(size(site_quat,1) == N, ...
        sprintf('The number of NPC label sites (%d) does not match the number of IgG (%d).', ...
        size(site_quat,1),N));
    local_quaternions(:,1) = site_quat;
end

min_dist_tol = 3; % Minimum distance between the label site and the NPC label site.

% primary antibody values
L_igg1 = {[6.527, 6.669, 6.819, 7.157]; [2.520, 7.239]};
phi_igg1 = {[-22.5, 22.5]; [15.0, 127.6]};
theta_igg1 = {[-22.5, 22.5]; [0, 360]};

% secondary antibody values
% values if AF647 is on Fab1
L_Fab1 = {[0.0, 7.0]; 2.3; 0};
phi_Fab1 = {[67.5, 112.5]; 90; 0};
theta_Fab1 = {[0, 360]; [0, 360]; 0};
     
% values if AF647 is on Fab2
L_Fab2 = {[6.527, 6.669, 6.819, 7.157]; [0.0, 7.0]; 2.3};
phi_Fab2 = {[67.5, 112.5]; 'phi4'; 90};
theta_Fab2 = {[0, 360]; [0, 360]; [0, 360]};

% values if AF647 is on Fc
L_Fc = {[6.527, 6.669, 6.819, 7.157]; [0.0, 8.5]; 2.1};
phi_Fc = {[67.5, 112.5]; [15.0, 127.6]; 90};
theta_Fc = {[0, 360]; [0, 360]; [0, 360]};

P = [0.31, 0.31, 0.38]; % probability of each segment of the IgG being labeled

L_igg2 = {L_Fab1, L_Fab2, L_Fc};
phi_igg2 = {phi_Fab1, phi_Fab2, phi_Fc};
theta_igg2 = {theta_Fab1, theta_Fab2, theta_Fc};
% we randomly sample the branch of the IgG that is labeled.
branch_igg2 = randsample(3,N,true,P);

for i_segment = 1:5 % loop through each segment of the IgG
    if i_segment <= 2
        % if the current segment is the primary IgG.
        L_ = get_rv(L_igg1{i_segment},N);
        phi_ = get_rv(phi_igg1{i_segment},N);
        theta_ = get_rv(theta_igg1{i_segment},N);
        quaternions_ = local_quaternions(:,i_segment);

        % get the next point and orientation of the next point
        [next_point, next_quaternion] = generateNextPoint(L_,phi_,theta_, ...
            quaternions_);
        
        xyz(:,:,i_segment+1) = xyz(:,:,i_segment) + next_point';
        local_quaternions(:,i_segment+1) = next_quaternion;

    else
        % if the current segment is the secondary IgG.
        next_quaternions = local_quaternions(:,i_segment+1);
        next_xyz = xyz(:,:,i_segment+1);

        min_dist = zeros(N,1);
        steric_npc = false(N,1);

        while any((min_dist < min_dist_tol) | steric_npc) && iterations < 10

            for i_branch = 1:3 % loop for each potential branch of the IgG
                % get the number of segments for the current branch
                update_sel = (branch_igg2 == i_branch) & ((min_dist < min_dist_tol) | steric_npc);
                N_sel = sum(update_sel);

                if N > 0
                    % get the current segment of the seonndary IgG
                    L_ = get_rv(L_igg2{i_branch}{i_segment-2},N_sel);
                    phi_ = get_rv(phi_igg2{i_branch}{i_segment-2},N_sel);
                    theta_ = get_rv(theta_igg2{i_branch}{i_segment-2},N_sel);
                    quaternions_ = local_quaternions(update_sel,i_segment);

                    % get the next point and orientation of the next point
                    [next_point, next_quaternion] = generateNextPoint(L_,phi_,theta_, ...
                        quaternions_);

                    next_xyz(update_sel,:) = xyz(update_sel,:,i_segment) + next_point';
                    next_quaternions(update_sel) = next_quaternion;
                end
            end
            
            % check the segment distances
            distances = check_distances(xyz,i_segment,next_xyz);
            min_dist = squeeze(min(distances,[],[1,3]))';
            iterations = iterations + 1;

            % check for any non-feasible locations
            steric_npc = check_steric_npc(xyz,local_quaternions,next_xyz);
        end

        if iterations >= 10
            warning('Could not find a feasible solution for the secondary antibody');
            update_sel = ~((min_dist < min_dist_tol) | steric_npc);
            xyz(update_sel,:,i_segment+1) = next_xyz(update_sel,:);
            local_quaternions(update_sel,i_segment+1) = next_quaternions(update_sel);
            
        else
            xyz(:,:,i_segment+1) = next_xyz;
            local_quaternions(:,i_segment+1) = next_quaternions;
        end

        % reset the iterations
        iterations = 0;
    end
end
fluor_pos = xyz;

end

%% helper functions
function distances = check_distances(xyz,curr_segment,next_xyz)
% This function checks the distances of the last line segments against the older line segments.
N_chains = size(xyz,1);
distances = zeros(N_chains,N_chains,curr_segment-1);

for i_chain = 1:N_chains
    for j_chain = 1:N_chains
        for i_segment = 1:curr_segment-1
            if i_chain == j_chain
                distances(i_chain,j_chain,i_segment) = Inf;
            else
                distances(i_chain,j_chain,i_segment) = shortest_distance_segment(...
                    xyz(j_chain,:,curr_segment),next_xyz(j_chain,:),...
                    xyz(i_chain,:,i_segment),next_xyz(i_chain,:));
            end
        end
    end
end

end

function steric = check_steric_npc(xyz,local_quaternions,next_xyz)
% this function checks if the antibody is stericly clashing with the NPC.
% we mark those points that are stericly clashing with the NPC as true.

N_chains = size(xyz,1);
steric = false(N_chains,1);

for i_chain = 1:N_chains
    % Get the quaternion and coordinates of the binding site
    q = local_quaternions{i_chain,1};
    p0 = xyz(i_chain,:,1);

    % point to check
    r = next_xyz(i_chain,:);

    % Convert the quaternion to a directional vector of the global coordinates
    v = quatrotate(compact(q), [0,0,1]);

    % Compute the signed distance of the point v from the plane
    signed_distance = dot(v,r-p0) / norm(v);

    % Check if the point is within the NPC
    if signed_distance <= 1.5
        steric(i_chain) = true;
    end
end
end

function rv = get_rv(val,n)
% This function returns a random value from a given range of values.
if nargin < 2
    n = 1;
end

if strcmpi(val, 'phi4') % special case for phi4
    rv = get_phi4(n);
elseif length(val) == 1 % return that value
    rv = repmat(val,n,1);
elseif length(val) == 2 % returns from the range of values
    rv = rand(n,1) .* (val(2) - val(1)) + val(1);
else % return a random value from the list
    rv = datasample(val,n)';
end

end

% pdf of this function based on Bongini et. al., 
% Proc. Natl. Acad. Sci. U.S.A., 2004, Appendix 2.
function phi4 = get_phi4(n)

sigma_theta = 91.6;
theta0 = 180;
phi_min = 15;
phi_max = 120; % corrected to 120, as phi_max > 120 is always complex.

phi4 = ones(n,1).* 1i;

while any(imag(phi4) ~= 0)
    n = sum(imag(phi4) ~= 0);
    
    phi1 = rand(n,1) * (phi_max - phi_min) + phi_min;
    phi2 = rand(n,1) * (phi_max - phi_min) + phi_min;
    theta = theta0 - sigma_theta .* tand((1 - 2.*rand(n,1)).*atand(theta0./sigma_theta));
    
    cos_phi1 = cosd(phi1);
    cos_phi2 = cosd(phi2);
    cos2_phi1 = cos_phi1 .* cos_phi1;
    cos2_phi2 = cos_phi2 .* cos_phi2;
    
    phi4_temp = acosd(sqrt(1+cos_phi1-2.*cos2_phi1).*sqrt(1+cos_phi2-2.*cos2_phi2).*cosd(theta)+cos_phi1.*cos_phi2);
    
    phi4(imag(phi4) ~= 0) = phi4_temp;
end

end

function shortest_distance = shortest_distance_segment(A, B, C, D)
% This function returns the shortest distance between two line segments.
% A, B, C, and D are the endpoints of the two line segments.
R = B - A;
S = D - C;
W = A - C;

t_numerator = dot(W, S) * dot(S, S) - dot(W, R) * dot(S, R);
t_denominator = dot(R, R) * dot(S, S) - dot(R, S) * dot(R, S);

if t_denominator ~= 0
    t = t_numerator / t_denominator;
    t = min(max(0, t), 1);
else
    t = 0;
end

s_numerator = t * dot(R, S) - dot(W, S);
s_denominator = dot(S, S);

if s_denominator ~= 0
    s = s_numerator / s_denominator;
    s = min(max(0, s), 1);
else
    s = 0;
end

closest_point_on_segment1 = A + t * R;
closest_point_on_segment2 = C + s * S;
shortest_distance = norm(closest_point_on_segment1 - closest_point_on_segment2);

end
