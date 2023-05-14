%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File to solve the Heterogeneous and Anisotropic problem from the        %
% paper                                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
close all;
clc;

addpath(genpath('./Functions_SLBM'));
addpath(genpath('./Solvers'));
addpath(genpath('./Figures'));

% Simulation parameters
scale_Dx = 1/200;
g = 0.7; 
choice_case = 1; % 0 - homogeneous media / 1 - heterogeneous media
Nb_speed_level = 2;

% Points of measure: backscattered flux
r_Nodes_meas = [ -2.5 , -0.5 ];
Nb_pt_meas = size(r_Nodes_meas,1);             

[speed_set,Nb_dir_tot] = construct_speed_set(Nb_speed_level);
[directions] = get_Directions(Nb_speed_level,speed_set);
[Weights_quad_MTL_iso,~] = get_Weights_trapez_HG(Nb_dir_tot,directions,0);
Weights_quad_MTL_iso = 2*pi*Weights_quad_MTL_iso;
[Weights_quad_MTL_aniso,~] = get_Weights_trapez_HG(Nb_dir_tot,directions,g);

[time_mesh_MTL,E_MTL,F_paper_MTL,time_instant] = ...
    solve_HetAni(scale_Dx,choice_case, ...
Nb_dir_tot,Nb_speed_level,speed_set, ...
Weights_quad_MTL_aniso,Weights_quad_MTL_iso, ...
Nb_pt_meas,r_Nodes_meas);

    

