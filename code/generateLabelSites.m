function [fluor_sites,label_sites] = generateLabelSites(npc_params,n_npc,n_sites)
% This function generates the label sites for a given list of NPC
% centroids.
% 
%
% Created by Weihong Yeo, Northwestern University, 2022-07-29.
% Last modified by Weihong Yeo, Northwestern University, 2022-10-02.
% 
% # Changelog
% ## 2022-08-19
% - added transmembrane and other NUPs
% 
% ## 2022-10-02
% - changed code to accomodate updated generateFluorophoreSite
% - added input parameter to pass n_sites

%%
% Generate fluorophore label sites on each NPC
% If the target site is on the coat NPC

%% parse location by location.
fluor_sites = nan(3,8,n_npc);
label_sites = nan(3,8,n_npc);

if nargin < 3
    n_sites = binornd(8,npc_params.prob_label,1,n_npc);
else
    assert(length(n_sites) == n_npc);
end

for i_npc = 1:n_npc
    n_site = n_sites(i_npc);
    idx_sites = randperm(8,n_site);
    npc_ang = (idx_sites + rand)*45;
    label_sites_ = npc_params.radial_dist .* [cosd(npc_ang); sind(npc_ang); zeros(size(npc_ang))];

    if any(npc_params.location_index == 2)  % for inner ring
        fluor_sites_ = nan(3,n_site);
        for i_site = 1:n_site
            fluor_sites_(:,i_site) = rotz(npc_ang(i_site)+90) * generateFluorophoreSite(1,npc_params);
        end
        
    elseif any(npc_params.location_index == 3)  % for transmembrane
        fluor_sites_ = nan(3,n_site);
        for i_site = 1:n_site
            fluor_sites_(:,i_site) = rotz(npc_ang(i_site)-90) * generateFluorophoreSite(1,npc_params);
        end
    else
        fluor_sites_ = generateFluorophoreSite(n_site,npc_params);
    end
    
    fluor_sites(:,1:n_site,i_npc) = fluor_sites_;
    label_sites(:,1:n_site,i_npc) = label_sites_;
end

end