function x = fastpoissondisk(domain_width,min_separation,n_points)
% Fast Poisson Disk Sampling method from Bridson, 2007
%
% Inputs:
%   domain_width: width of the domain (nm) 
%   min_separation: minimum separation between points (nm)
%   n_points: number of points to generate (optional)
%
% Outputs:
%   x: 2-column matrix with the x and y coordinates of the points
%
% 
% Written by Weihong Yeo, Northwestern University, 2022-08-19.
%
% ######################################################################### 
% Changelog
% #########################################################################
% 2022-08-19: Created file.
%

% Check inputs
grid_size = min_separation / sqrt(2);
if nargin < 3 || isempty(n_points)
    n_points = round(((domain_width / grid_size) ^ 2) / 2);
end

% Initialize output number of points
x = nan(n_points,2);
x(1,:) = (rand(1,2)-0.5).*grid_size;

% Generate points
kk = 1; jj = 2;
while (kk < n_points) && (jj < n_points)
    % Generate 30 points in a circle around x(kk,:)
    t = (rand(30,2).*[min_separation,2*pi])+[min_separation,0];
    y = [t(:,1).*cos(t(:,2)),t(:,1).*sin(t(:,2))]+x(kk,:);
    
    % Check if the points are within the domain and far enough from each other
    for ii = 1:30
        dist = pdist2(x,y(ii,:));
        if ~any(dist < min_separation,1) && all(abs(y(ii,:)) < (domain_width / 2))
            x(jj,:) = y(ii,:);
            jj = jj + 1;
        end
    end
    kk = kk + 1;
end

% Remove NaNs that were initialized
x = x(~isnan(x(:,1)),:);

end