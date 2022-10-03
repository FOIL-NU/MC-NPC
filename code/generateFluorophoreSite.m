function fluor_pos = generateFluorophoreSite(N,npc_params)
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

%% Set default parameters

if (nargin == 0) || isempty(N)
    N = 1;
end

if (nargin <= 1) || isempty(npc_params.location_index)
    npc_params.location_index = 1; % outer_membrane
end

del_phi = 22.5;
phi0 = [-1,1].* 45;

% values from https://pubs.acs.org/doi/10.1021/acsnano.1c03677
L = {[6.527, 6.669, 6.819, 7.157]; ...
     [2.520, 7.239]; ...
     [6.527, 6.669, 6.819, 7.157]; ...
     [3.127, 3.898, 5.523, 6.072]; ...
     [0.1, 0.63]};

phi = {phi0;
       [min(phi0-del_phi), max(phi0-del_phi)];
       [15,127.6];
       [min(90-del_phi), max(90+del_phi)];
       [15, 127.6];
       [10, 170]};

theta = {[-45, 45];
         [0, 360];
         [0, 360];
         [0, 360];
         [-80, 80]};

phi0_ = rand(N,1) * (phi0(2) - phi0(1)) + phi0(1);
y0 = 10.5 + 2.75 .* cosd(phi0_);
z0 = 2.75 .* sind(phi0_);

L_ = zeros(N,5);
phi_ = zeros(N,5);
theta_ = zeros(N,5);
x_ = zeros(N,5);
y_ = zeros(N,5);
z_ = zeros(N,5);

x = zeros(N,6);
y = zeros(N,6);
z = zeros(N,6);

for ii = 5:-1:1
    [L_(:,ii),phi_(:,ii),theta_(:,ii),x_(:,ii),y_(:,ii),z_(:,ii)] = get_xyz(L{ii},phi{ii},theta{ii},N);
    if ii < 5
        for jj = 1:N
            xyz = roty(theta_(jj,ii)) * rotx(phi_(jj,ii)) * [x(jj,ii+1:end)',y(jj,ii+1:end)',z(jj,ii+1:end)']';
            x(jj,ii+1:end) = xyz(1,:)' + x_(jj,ii);
            y(jj,ii+1:end) = xyz(2,:)' + y_(jj,ii);
            z(jj,ii+1:end) = xyz(3,:)' + z_(jj,ii);
        end
    else
        x(:,ii+1) = x_(:,ii);
        y(:,ii+1) = y_(:,ii);
        z(:,ii+1) = z_(:,ii);
        
    end
end

y = y + y0;
z = z + z0;

fluor_pos = [x(:,end),y(:,end),z(:,end)]';

%% orientate the primary and secondary antibodies to based on location of label
% we do not need to do anything if location is 'outer' or 'nuclear'.

if ~any(npc_params.location_index == [2,3])
    fluor_pos = rotx(90) * fluor_pos;
end

end

%% helper function
function [L_,phi_,theta_,x_,y_,z_] = get_xyz(L,phi,theta,N)

if nargin < 4
    N = 1;
end


if length(L) < 2
    L_ = repmat(L,N,1);
elseif length(L) == 2
    L_ = rand(N,1) .* (L(2) - L(1)) + L(1);
else
    L_ = datasample(L,N)';
end

if length(phi) < 2
    phi_ = repmat(phi,N,1);
else
    phi_ = rand(N,1) .* (phi(2) - phi(1)) + phi(1);
end

if length(phi) < 2
    theta_ = repmat(theta,N,1);
else
    theta_ = rand(N,1) .* (theta(2) - theta(1)) + theta(1);
end

x_ = L_ .* sind(phi_) .* sind(theta_);
y_ = L_ .* cosd(phi_);
z_ = L_ .* sind(phi_) .* cosd(theta_);

end
