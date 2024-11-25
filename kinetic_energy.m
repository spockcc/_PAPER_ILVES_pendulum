function y=kinetic_energy(dim, n, m, v)

% Calculates the kinetic energy of a collection of particles
%
% CALL SEQUNCE: y=kinetic_energy(dim, n, m, v)
%
% INPUT:
%   dim  the dimension of the vector space
%   n    the number the number of particles
%   m    m(i) is the mass of particle i
%   v    v(dim*(i-1)+1:dim*i) is the velocity of particle i
% 
% OUTPUT:
%   y   the total kinetic energy of the particles
%
% The kinetic energy is simply the sum
% 
%     y = 0.5*sum_{i=1}^n m(i)*norm(v(i),2)^2
% 
% where v(i) is the velocity of particle i.
%
% MINIMAL WORKING EXAMPLE: TODO

% PROGRAMMING by Carl Christian Kjelgaard Mikkelsen (spock@cs.umu.se)
%   2022-01-03 Initial programming and testing
%   2022-08-03 Added support for finite dimensional vector spaces

% Initialize the output
y=0;
% Loop over all particles
for i=1:n
    % Isolate the indices for particle i
    ii=dim*(i-1)+1:dim*i;
    % Isolate the velocity of particle i
    vi=v(ii);
    % Add the kinetic energy of particle i
    y=y+m(i)*norm(vi,2)^2;
end
% Perform the final scaling
y=0.5*y;