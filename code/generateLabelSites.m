function fluor_sites = generateLabelSites(npc_params,n_npc,n_sites)
% This function generates the label sites for a given list of NPC
% centroids.
% 
%
% Created by Weihong Yeo, Northwestern University, 2022-07-29.
% Last modified by Weihong Yeo, Northwestern University, 2023-04-15.
% 
%
% ######################################################################### 
% Changelog
% ######################################################################### 
% 2023-04-15:
% - changed code to use minor angle and radius instead of using location
%
% 2022-10-02:
% - changed code to accomodate updated generateFluorophoreSite
% - added input parameter to pass n_sites
%
% 2022-08-19: added transmembrane and other NUPs
% 

% Generate fluorophore label sites on each NPC
% If the target site is on the coat NPC

% parse location by location.
fluor_sites = nan(3,8,n_npc);

if nargin < 3
    n_sites = binornd(8,npc_params.prob_label,1,n_npc);
else
    assert(length(n_sites) == n_npc);
end

% loop for each NPC
for i_npc = 1:n_npc
    n_site = n_sites(i_npc);

    if n_site > 0
        % get the number of sites for each NPC
        idx_sites = randperm(8,n_site);
        npc_ang = (idx_sites + rand)*45;

        % generate coodinates for each label site
        label_sites_ = npc_params.radial_dist .* ...
            [sind(npc_ang); cosd(npc_ang); zeros(size(npc_ang))]';

        % determine the quaternion for each site
        label_quaternions_ = cell(n_site,1);
        q_t = quaternion(cosd(npc_params.minor_angle/2), ...
                        -sind(npc_params.minor_angle/2),0,0);
        for i1 = 1:n_site
            q_z = quaternion(cosd(npc_ang(i1)/2), ...
                            0,0,sind(npc_ang(i1)/2));
            label_quaternions_{i1} = q_t * q_z;
        end
        
        % generate fluorophore sites
        fluor_sites_ = generateFluorophoreSites(label_sites_,label_quaternions_);
        
        % swap the order of axis 1 and 2, and save it to fluor_sites
        fluor_sites_ = permute(fluor_sites_,[2 1 3]);
        fluor_sites(:,1:n_site,i_npc) = fluor_sites_(:,1:n_site,end);
    end
end

end