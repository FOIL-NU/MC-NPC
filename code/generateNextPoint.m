function [next_point, next_quaternion] = generateNextPoint(link_length, ...
    local_polar, local_azimuth, local_quaternion)
% This function computes the global coordinates of the next point in a
% chain of links, given the local polar and azimuthal angles of the
% at the current point. The local polar angle is the angle away from the
% z-axis, and the local azimuthal angle is the angle away from the x-axis.
% The local angles are defined in the local coordinate system of the
% previous link. 
%
% Note that this function does not assume that the links are continuous, so
% this function may be used to generate for multiple next points in 
% independent chains.
%
% This function requires either the Robotics System Toolbox and/or 
% Aerospace Toolbox from MATLAB to use the functions `quaternion` and
% `quatrotate`.
%
% Inputs:
%   link_length: array of lengths of the next link
%   local_polar: array of local polar angles of the next link
%   local_azimuth: array of local azimuthal angles of the next link
%   local_quaternion: cell array of local quaternions of the current points
%
% Outputs:
%   next_point: A 3xN vector of the global coordinates of the next point 
%   next_quaternion: cell array of quaternions of the next point
%
%
% Created by Weihong Yeo, Northwestern University, 2023-04-13
% 
% #########################################################################
% Changelog
% #########################################################################
% 2023-04-13: Created file.
% 

% Check inputs
if nargin < 3
    error('Not enough inputs.')
end

N = length(link_length);

assert(length(local_polar) == N, ...
    'The number of local polar angles must be equal to the number of link lengths.')
assert(length(local_azimuth) == N, ...
    'The number of local azimuthal angles must be equal to the number of link lengths.')
assert(length(local_quaternion) == N, ...
    'The number of local quaternions must be equal to the number of links.')

% Initialize variables
next_point = zeros(3, N);
if nargin < 4 || isempty(local_quaternion)
    local_quaternion = cell(1, N);
    local_quaternion{1} = quaternion(1, 0, 0, 0);
end

next_quaternion = cell(1, N);

for ii = 1:N
    q_global = local_quaternion{ii};
    
    % Check if quaternion is valid
    if isempty(q_global)
        next_point(:, ii) = nan;
    else
        % Compute quaternions for local polar and azimuthal angles
        q_polar = quaternion(cosd(local_polar(ii)/2), sind(local_polar(ii)/2), 0, 0);
        q_azimuth = quaternion(cosd(local_azimuth(ii)/2), 0, 0, sind(local_azimuth(ii)/2));

        % Convert global angles to local angles
        q_polar_global = conj(q_global) * q_polar * q_global;
        q_azimuth_global = conj(q_global) * q_azimuth * q_global;

        % Update the global rotation quaternion
        q_global = q_global * q_polar_global * q_azimuth_global;
        next_quaternion{ii} = q_global;

        % Calculate the rotated vector
        rotated_vector = quatrotate(compact(q_global), [0, 0, 1]);

        % Add the rotated vector to the last point in the chain
        next_point(:, ii) = link_length(ii) * rotated_vector';
    end
end
end
