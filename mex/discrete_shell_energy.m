% [E,G,H] = discrete_shell_energy(V,F,...)
%
% Inputs:
%   V  #V by 3 
%   F  #F by 3 
%   Optional:
%     'U0' followed by current displacements
%     'Poissons' followed by the Poisson's ratio
%     'Youngs'  followed by Young's modulus
%     'SFF'  followed by one of 'sin', 'tan', 'average'
%     'Membrane'  followed by one of 'stvk' or 'neohookian'
%     'Thickness'  followed by thickness value in meters (todo: this could be
%       per element using libshell)
% Outputs:
%   E  scalar energy value
%   G  #V*3 by 1 gradient vector
%   H  3*#V by 3*#V Hessian sparse matrix
%
