% Adapted from prefLoc.m by Simon Schultz
% This script finds the preferred location of a place cell by
%   phi_pref1 : calculating vector average place preference
%   phi_pref2 : finding peak of place field map (firing rate/occupancy)

function [ phi_pref1, phi_pref2 ] = prefLoc( pfMap )

Ncells = size( pfMap, 1 );
Nbins = size( pfMap, 2 );
phi = (1:Nbins)*360/Nbins;

phi_pref1 = zeros( Ncells, 1 );
phi_pref2 = zeros( Ncells, 1 );

for i = 1:Ncells
    
    % method 1: vector average place preference
    x = pfMap(i,:).*cosd(phi);
    y = pfMap(i,:).*sind(phi);
    nz = length(x>0);
    x_av = sum(x)/nz;
    y_av = sum(y)/nz;
    phi_pref1(i) = atan2d(y_av,x_av);
    
    % method 2: just find peak of firing rate profile
   [~,maxloc] = max( pfMap(i,:) );
   phi_pref2(i) = maxloc;
   
end

% method 1
phi_pref1(phi_pref1<0) = phi_pref1(phi_pref1<0) + 360;
phi_pref1 = discretize(phi_pref1,[0 phi]);

