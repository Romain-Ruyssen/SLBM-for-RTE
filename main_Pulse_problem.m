%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  File to solve the Pulse problem from the paper                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
close all;
clc;

addpath(genpath('./Functions_SLBM'));
addpath(genpath('./Solvers'));
addpath(genpath('./Figures'));

% Simulation parameters
scale_Dx = 1/200;
omega = 0.9;   


Nb_speed_level = 2;
[speed_set,Nb_dir_tot] = construct_speed_set(Nb_speed_level);
[directions] = get_Directions(Nb_speed_level,speed_set);
[Weights_quad_MTL,angles_MTL] = get_Weights_trapez_HG(Nb_dir_tot,directions,0);
[Weights_HCR_MTL] = get_Weights_trapez_Half_Circle([1,0],Nb_dir_tot,directions);
[Weights_HCL_MTL] = get_Weights_trapez_Half_Circle([-1,0],Nb_dir_tot,directions);

[time_mesh, q_T, q_R] = ...
       solve_Pulse(scale_Dx,Nb_dir_tot,Nb_speed_level,speed_set,Weights_quad_MTL,Weights_HCR_MTL,Weights_HCL_MTL,omega);


return
