%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File to solve the Laser problem from the paper                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
close all;
clc;

addpath(genpath('./Functions_SLBM'));
addpath(genpath('./Solvers'));
addpath(genpath('./Figures'));

% Simulation parameters
scale_Dx = 1/200;
Nb_speed_level = 2;
g = 0; 
case_choice = 0;        % Cases from the paper
% case_choice = 1;

% Points of measure
r_Nodes_meas_x = (0:scale_Dx*1:1)';
Nb_nodes = size(r_Nodes_meas_x,1);
r_Nodes_meas_y = 0.5*ones(Nb_nodes,1);
r_Nodes_meas = [ r_Nodes_meas_x , r_Nodes_meas_y ];
Nb_pt_meas = size(r_Nodes_meas,1);             


[speed_set,Nb_dir_tot] = construct_speed_set(Nb_speed_level);
[directions] = get_Directions(Nb_speed_level,speed_set);
[Weights_quad_MTL_iso,~] = get_Weights_trapez_HG(Nb_dir_tot,directions,0);
Weights_quad_MTL_iso = 2*pi*Weights_quad_MTL_iso;
[Weights_quad_MTL_aniso,~] = get_Weights_trapez_HG(Nb_dir_tot,directions,g);

[time_mesh,E,time_instant] = solve_Laser(case_choice, ...
            scale_Dx,Nb_dir_tot, ...
            Nb_speed_level,speed_set, ...
            Weights_quad_MTL_aniso,Weights_quad_MTL_iso, ...
            Nb_pt_meas,r_Nodes_meas);

   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Printing the energy in the section
if case_choice == 0
    sp = 150;
    lp = 2;
    figure
    hold on
    plot([0.49,0.5,0.5,0.52],[1.0,1.0,0.0,0.0],'LineWidth',2,'Color','b');
    plot([0.49,time_instant(172,1),time_instant(172,1),0.52],[1.0,1.0,0.0,0.0],'LineWidth',2,'Color','r');
    plot([0.49,0.505,0.505,0.52],[1.0,1.0,0.0,0.0],'LineWidth',2,'Color','g');
    scatter(r_Nodes_meas_x(99:105,1),E(99:105,171),sp,'o',"LineWidth",lp,'MarkerEdgeColor','b');
    scatter(r_Nodes_meas_x(99:105,1),E(99:105,172),sp,'*',"LineWidth",lp,'MarkerEdgeColor','r');
    scatter(r_Nodes_meas_x(99:105,1),E(99:105,173),sp,'^',"LineWidth",lp,'MarkerEdgeColor','g');
    Legend = ["exact $t = 0.5$" ; "exact $t = 0.502$" ; "exact $t = 0.505$" ; ...
              "SLBM $t = 0.5$" ; "SLBM $t = 0.502$" ; "SLBM $t = 0.505$" ];
    legend(Legend,'Interpreter','latex');
    hold off
    xlabel('$x$','Interpreter','latex');
    ylabel('$E$','Interpreter','latex');
    set(gcf,'Position',[200 100 900 600]);
    set(gca,'color','w'); 
    set(gcf,'color','w'); 
    set(gca,'FontSize',18);
    xlim([0.49,0.52]);
    ylim([0,1]);
end










