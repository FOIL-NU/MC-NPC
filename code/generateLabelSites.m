function [fluor_sites,label_sites] = generateLabelSites(npc_params,n_npc)
% This function generates the label sites for a given list of NPC
% centroids.
% 
%
% Created by Weihong Yeo, Northwestern University, 2022-07-29.
% Last modified by Weihong Yeo, Northwestern University, 2022-08-19.
% 
% # Changelog
% ## 2022-08-19
% - added transmembrane and other NUPs

%%
% Generate fluorophore label sites on each NPC
% If the target site is on the coat NPC

%% parse location by location.
fluor_sites = nan(3,8,n_npc);
label_sites = nan(3,8,n_npc);

for i_npc = 1:n_npc
    n_sites = binornd(8,npc_params.prob_label);
    idx_sites = randperm(8,n_sites);
    npc_ang = (idx_sites + rand)*45;
    label_sites_ = npc_params.radial_dist .* [cosd(npc_ang); sind(npc_ang); zeros(size(npc_ang))];

    if any(npc_params.location_index == [2,3])
        fluor_sites_ = nan(3,8);
        for i_site = 1:n_sites
            fluor_sites_(:,i_site) = rotz(npc_ang(i_site)) * generateFluorophoreSite(n_sites,npc_params);
        end
    else
        fluor_sites_ = generateFluorophoreSite(n_sites,npc_params);
    end
    
    fluor_sites(:,1:n_sites,i_npc) = fluor_sites_;
    label_sites(:,1:n_sites,i_npc) = label_sites_;
end

end