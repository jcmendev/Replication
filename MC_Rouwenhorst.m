function [zgrid,Pz,dz] = MC_Rouwenhorst(nz,std_z,mu_z,rho_z)
% % Function MC_Rouwenhorst: Discretizes an AR(1) into a Markov Chain,
% following Rouwenhorst (1995) method
% 
% USAGE
%   Inputs: 
%           nz    : Number of nodes the exogenous states        
%           std_z : Standard deviation of the exogenous state
%           mu_z  : Mean of the exogenous state
%           rho_z : Autocorrelation of the exogenous state
%   Outputs:
%           zgrid  : Grid with exogenous states, nzx1
%           Pz     : Probability transition Matrix (PTM). nzxnz. Rows are current
%                    states, columns are future states
%           dz     : Ergodic distribution of exogenous states. nzx1
% Juan C. Mendez-Vizcaino, 2024
%% Technical parameters for the ergodic distribution
maxit = 1000;
tol   = 1e-10;
%% Build the grid of exogenous states
z_min = -sqrt(nz-1)*std_z/sqrt(1-rho_z.^2);
z_max = -z_min; 
zgrid = linspace(z_min,z_max,nz)';
zgrid = zgrid+mu_z;
%% Build the PTM
p  = (1+rho_z)/2;
Pz = [p,1-p;1-p,p];
for iz = 3:nz
    Pz = p*[Pz,zeros(iz-1,1);zeros(1,iz)     ]+...
    (1-p)*[zeros(iz-1,1) Pz ;zeros(1,iz)     ]+...
    (1-p)*[zeros(1,iz)      ;Pz,zeros(iz-1,1)  ]+...
        p*[zeros(1,iz)      ;zeros(iz-1,1) Pz  ];

        Pz(2:end-1,:) = Pz(2:end-1,:)/2;
end    
%% Build the ergodic distribution, starting from a uniform guess.
d0 = (1/nz)*ones(nz,1);
for it=1:maxit
dz = Pz'*d0;
if norm(d0-dz)<tol; break; end
d0 = dz;
end
