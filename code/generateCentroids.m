function npc_centroids = generateCentroids(membrane_width,n_npc,min_separation,keepout_width)
% This function generates the centroids of nuclear pore complexes on a
% planar membrane. This method uses the fast Poisson Disk Sampling method
% from Bridson, 2007.
% 
% The final output has size of (3,n_npc), whereby
% - 3: x, y, z
% - n_npc: number of NPCs
%
% Created by Weihong Yeo, Northwestern University, 2022-07-29.
% Last modified by Weihong Yeo, Northwestern University, 2022-08-19.
% 
% # Changelog
% ## 2022-08-19
% - added checks for inputs

%% Check inputs
assert(isreal(membrane_width));
if nargin <= 1  % if no. of npcs is not provided, assume density of 3/um2
    n_npc = 3 * (membrane_width / 1000).^2;
end
if nargin <= 2  % if min_separation is not provided, assume 250 um.
    min_separation = 250;
end
if nargin <= 3  % if keepout width is not provided, assume 100 um.
    keepout_width = 100;
end

%% Generate NPC centroid locations
% Use the fast Poisson Disk Sampling method from Bridson, 2007.
npc_centroids = fastpoissondisk(membrane_width,min_separation);

%% Remove points that are out of keepout width
npc_centroids = npc_centroids(all(abs(npc_centroids)<(membrane_width/2-keepout_width),2),:);
npc_centroids = npc_centroids+(membrane_width/2);

%% Randomly select available points
npc_centroids = npc_centroids(randperm(size(npc_centroids,1),n_npc),:);

%% Reshape the output to the prescribed size
npc_centroids = [npc_centroids';zeros(1,size(npc_centroids,1))];

end