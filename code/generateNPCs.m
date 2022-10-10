function tbl = generateNPCs(nup_species,membrane_width,npc_density,bg_density,visualize_data)
% This function generates NUPs on a planar nuclear membrane.
% 
% membrane_size is real-valued specifying a width in nm.
% npc_density is an integer from 1 to 9 NPCs/um2.
% bg_density is an integer that specifies the background noise.
%
%
% Created by Weihong Yeo, Northwestern University, 2022-08-19.
%

%% check inputs
if nargin < 1
    nup_species = "nup133";
end
if nargin < 2
    membrane_width = 10000;
end
if nargin < 3
    npc_density = 9;
end
if nargin < 4
    bg_density = 10;
end
if nargin < 5
    visualize_data = true;
end

%%
nup_parameters = getNupParameters(nup_species);
n_npc = npc_density * (membrane_width / 1000).^2;

centroids = generateCentroids(membrane_width,n_npc);
centroids = reshape(centroids,3,1,[]);

output = [];

[fluor_sites,label_sites] = generateLabelSites(nup_parameters,n_npc);

for i_npc = 1:n_npc
    n_sites = sum(~isnan(label_sites(1,:,i_npc)));
    
    for i_site = 1:n_sites
        localizations = generateLocalizations(nup_parameters,15);
        
        abs_xyz = centroids(:,:,i_npc) + fluor_sites(:,i_site,i_npc) + ...
            label_sites(:,i_site,i_npc) + localizations;
              
        length_append = size(localizations,2);
        output = [output; abs_xyz(1:2,:)',(1:length_append)',ones(length_append,2).*[i_site,i_npc]];
    end
end

if bg_density > 0
    bg_output = generateBackground(bg_density,membrane_width);
    output = [output; bg_output, zeros(size(bg_output,1),3)];
end

tbl = array2table(output,'VariableNames',{'x [nm]','y [nm]','frame','npc_site','npc_index'});
writetable(tbl,filename);

if visualize_data
    visualize(tbl,nup_species,membrane_width,centroids,centroids+label_sites);
end
