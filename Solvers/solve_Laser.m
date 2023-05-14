function [time_mesh,E_meas,time_instant] = solve_Laser(case_choice, ...
                scale_Dx,Nb_dir_tot, ...
                Nb_speed_level,speed_set, ...
                Weights_quad_aniso,Weights_quad_iso, ...
                Nb_pt_meas,r_Nodes_meas)
% Fonction pour résoudre le problème de la propagation du laser
% en utilisant l'approche MTL

% Définition des paramètres:
c = 1;
L_x = 1; % (m) length of the slab
L_y = 1;

T_period = (L_x)/c;

% Discrétisation du domaine spatial:
x_begin = 0;
y_begin = 0;

% Précision géométrique
e = 10^(-10);

% Lattice
Delta_x = scale_Dx*L_x;
Nb_nodes_x = L_x/Delta_x + 1;
Nb_nodes_y = L_y/Delta_x + 1;
Nb_nodes = Nb_nodes_x*Nb_nodes_y;
X_Nodes = (0:Delta_x:L_x) + x_begin;
Y_Nodes = (0:Delta_x:L_y) + y_begin;
r_Nodes = zeros(Nb_nodes,2);
Num_mat_Nodes = zeros(Nb_nodes,2);
X_grid_Nodes = zeros(Nb_nodes_x,Nb_nodes_y);
Y_grid_Nodes = zeros(Nb_nodes_x,Nb_nodes_y);
for j_nodes = 1:Nb_nodes_y
    for i_nodes = 1:Nb_nodes_x
        num_node_cons = (j_nodes - 1)*Nb_nodes_x + i_nodes;
        r_Nodes(num_node_cons,:) = [X_Nodes(i_nodes),Y_Nodes(j_nodes)];
        Num_mat_Nodes(num_node_cons,1) = i_nodes; 
        Num_mat_Nodes(num_node_cons,2) = j_nodes;
        X_grid_Nodes(i_nodes,j_nodes) = X_Nodes(i_nodes);
        Y_grid_Nodes(i_nodes,j_nodes) = Y_Nodes(j_nodes);
    end
end

% Construction de la lattice:
Lattice = get_lattice(Nb_speed_level,speed_set,Nb_nodes,Nb_nodes_x,Nb_nodes_y,0,0);

% Définition de la discrétisation temporelle
Delta_t = Delta_x/c;

time_steps = zeros(Nb_speed_level,1);
for i_speeds = 1:Nb_speed_level
    time_steps(i_speeds,1) = (Delta_x/c)*speed_set(i_speeds,1).norm;
end
% T_period = 100*max(time_steps);
% T_period = (L_x/2)/2;
time_mesh = construct_time_mesh(T_period,time_steps);
Nb_time_step = time_mesh.Nb_nodes;
sub_time_meshes = get_sub_time_meshes(Nb_speed_level,time_mesh);

disp(strcat("Nb of nodes: ",num2str(Nb_nodes)));
disp(strcat("Nb of time steps: ",num2str(Nb_time_step)));

% Obtention des noeuds appariable
Matchable = get_Matchable(Lattice,speed_set,Nb_speed_level);

% Caractérisation du media
S_a = 0;
S_s = 0;
S_t = S_a + S_s; % total scattering coefficient

% Caractérisation des frontières

% left border: x = 0
border_1.Num_nodes = find(r_Nodes(:,1) <= x_begin + 10^(-10));
border_1.Nb_nodes = size(border_1.Num_nodes,1);
border_1.Normal_in = [1 ; 0];
for i_levels = 1:Nb_speed_level
    border_1.Ps_dir_ode(i_levels,1).values = zeros(speed_set(i_levels,1).Nb_speed,1);
    for i_speeds = 1:speed_set(i_levels,1).Nb_speed
        border_1.Ps_dir_ode(i_levels,1).values(i_speeds,1) = (1/speed_set(i_levels,1).norm)*( (speed_set(i_levels,1).speeds(i_speeds,:)) * (border_1.Normal_in) );
    end
end

% bot border: y = 0
border_2.Num_nodes = find(r_Nodes(:,2) <= y_begin + 10^(-10));
border_2.Nb_nodes = size(border_2.Num_nodes,1);
border_2.Normal_in = [0 ; 1];
for i_levels = 1:Nb_speed_level
    border_2.Ps_dir_ode(i_levels,1).values = zeros(speed_set(i_levels,1).Nb_speed,1);
    for i_speeds = 1:speed_set(i_levels,1).Nb_speed
        border_2.Ps_dir_ode(i_levels,1).values(i_speeds,1) = (1/speed_set(i_levels,1).norm)*( (speed_set(i_levels,1).speeds(i_speeds,:)) * (border_2.Normal_in) );
    end
end

% right border: x = L_x
border_3.Num_nodes = find(r_Nodes(:,1) >= x_begin + L_x - 10^(-10));
border_3.Nb_nodes = size(border_3.Num_nodes,1);
border_3.Normal_in = [-1 ; 0];
for i_levels = 1:Nb_speed_level
    border_3.Ps_dir_ode(i_levels,1).values = zeros(speed_set(i_levels,1).Nb_speed,1);
    for i_speeds = 1:speed_set(i_levels,1).Nb_speed
        border_3.Ps_dir_ode(i_levels,1).values(i_speeds,1) = (1/speed_set(i_levels,1).norm)*( (speed_set(i_levels,1).speeds(i_speeds,:)) * (border_3.Normal_in) );
    end
end

% top border: y = L_y
border_4.Num_nodes = find(r_Nodes(:,2) >= y_begin + L_y - 10^(-10));
border_4.Nb_nodes = size(border_4.Num_nodes,1);
border_4.Normal_in = [0 ; -1];
for i_levels = 1:Nb_speed_level
    border_4.Ps_dir_ode(i_levels,1).values = zeros(speed_set(i_levels,1).Nb_speed,1);
    for i_speeds = 1:speed_set(i_levels,1).Nb_speed
        border_4.Ps_dir_ode(i_levels,1).values(i_speeds,1) = (1/speed_set(i_levels,1).norm)*( (speed_set(i_levels,1).speeds(i_speeds,:)) * (border_4.Normal_in) );
    end
end

% Recherche des noeuds de la frontière où l'on impose le flux du laser
e = 10^(-10);
if case_choice == 0
    
    check_1 = find(r_Nodes(:,1) < x_begin + e);
    check_2 = find( (r_Nodes(check_1,2) > 0.25 - e) & (r_Nodes(check_1,2) < 0.75 + e));
    Num_nodes_bd_left = check_1(check_2);
    
elseif case_choice == 1    
    
    check_1 = find(r_Nodes(:,1) < x_begin + e);
    check_2 = find( (r_Nodes(check_1,2) > 0.90 - e) );
    Num_nodes_bd_left = check_1(check_2);
    
    check_3 = find(r_Nodes(:,2) > y_begin + L_y - e);
    check_4 = find( (r_Nodes(check_3,1) < 0.10 + e) );
    
    Num_nodes_bd_top = check_3(check_4);
    
end

for i_levels = 1:Nb_speed_level
    Nb_dir = speed_set(i_levels,1).Nb_speed;
    I(i_levels,1).values = zeros(Nb_dir,Nb_nodes);
    I_old(i_levels,1).values = zeros(Nb_dir,Nb_nodes);
end

I_0 = 1;
% Initialisation des champs inconnus
if case_choice == 0
    
    num_dir_glob_cons = speed_set(1,1).num_speeds_glob(1,1); 
    W_cons = Weights_quad_iso(1,num_dir_glob_cons);
    
    I(1,1).values(1,Num_nodes_bd_left) = (1/W_cons)*I_0;
    I_old(1,1).values(1,Num_nodes_bd_left) = (1/W_cons)*I_0;
    
elseif case_choice == 1
    
    num_dir_glob_cons = speed_set(2,1).num_speeds_glob(4,1); 
    W_cons = Weights_quad_iso(1,num_dir_glob_cons);
    
    I(2,1).values(4,Num_nodes_bd_left) = (1/W_cons)*I_0;
    I_old(2,1).values(4,Num_nodes_bd_left) = (1/W_cons)*I_0;
    
    I(2,1).values(4,Num_nodes_bd_top) = (1/W_cons)*I_0;
    I_old(2,1).values(4,Num_nodes_bd_top) = (1/W_cons)*I_0;
    
    % Pour la condition initiale du volume
    check = [];
    for i_nodes = 1:Nb_nodes
        if ( r_Nodes(i_nodes,2) > ( r_Nodes(i_nodes,1) + 0.9) - e )
            check = [check ; i_nodes];
        end
    end
    
    I(2,1).values(4,check) = (1/W_cons)*I_0;
    I_old(2,1).values(4,check) = (1/W_cons)*I_0;
    
end

% Initialisation des données que nous relevons aux points de mesure
E_meas = zeros(Nb_pt_meas,Nb_time_step);

   
% Résolution du problème:
%%%%%%%%%%%%%%%%%%%%%%%%%

% Algorithme de résolution un peu particulier
Elems_unknown_levels = ones(Nb_speed_level,1); % Initialisation du numéro des éléments inconnus de chaque niveau temporel
state_levels = check_state_speed_levels(Nb_speed_level,sub_time_meshes,Elems_unknown_levels);
Num_calcule = 0;
% while state_levels
for i_time_step = 1:time_mesh.Nb_nodes-1    
    Num_calcule = Num_calcule + 1;
    
    if (i_time_step == floor(time_mesh.Nb_nodes/4)) || ...
       (i_time_step == floor(time_mesh.Nb_nodes/2)) || ...
       (i_time_step == floor(time_mesh.Nb_nodes*(3/4)))
        disp(strcat("Computing: ",num2str(Num_calcule),"/",num2str(time_mesh.Nb_nodes)));
    end
    
    % Détermination des niveaux de vitesse à calculer:
    Num_nodes_known = zeros(Nb_speed_level,1);
    for i_levels = 1:Nb_speed_level
        Nb_elems_level_cons = sub_time_meshes(i_levels,1).Nb_elems;
        if Elems_unknown_levels(i_levels,1) <= Nb_elems_level_cons
            Num_nodes_known(i_levels,1) = sub_time_meshes(i_levels,1).Elems(Elems_unknown_levels(i_levels,1),1);
        else
            Num_nodes_known(i_levels,1) = sub_time_meshes(i_levels,1).Elems(Nb_elems_level_cons,2);
        end
    end
    Num_node_begin_unknown = min(Num_nodes_known);
    Time_begin_unknown = time_mesh.Nodes(Num_node_begin_unknown,1);
    Levels_to_compute = find(Num_nodes_known == Num_node_begin_unknown);
    
    % Calcul des intensités interpolées (ou non)
    I_int_tot = zeros(Nb_dir_tot,Nb_nodes); 
    for i_levels = 1:Nb_speed_level
        if ismember(i_levels,Levels_to_compute) % Pas besoin d'interpolation
            Nb_speed_lvl_cons = speed_set(i_levels,1).Nb_speed;
            for i_speeds = 1:Nb_speed_lvl_cons
                num_speed_lvl_cons = speed_set(i_levels,1).num_speeds_glob(i_speeds,1);
                I_int_tot(num_speed_lvl_cons,:) = I(i_levels,1).values(i_speeds,:); 
            end
        else % Besoin d'une interpolation qui sera ici linéaire
            Last_Elem_known = Elems_unknown_levels(i_levels,1) - 1;
            Num_nodes_elem_cons = sub_time_meshes(i_levels,1).Elems(Last_Elem_known,:);
            Times_cons = time_mesh.Nodes(Num_nodes_elem_cons,1);
            
            I_interpolated = -((Time_begin_unknown - Times_cons(2))/(Times_cons(2) - Times_cons(1)))*I_old(i_levels,1).values + ...
                              ((Time_begin_unknown - Times_cons(1))/(Times_cons(2) - Times_cons(1)))*I(i_levels,1).values;
            Nb_speed_lvl_cons = speed_set(i_levels,1).Nb_speed;
            for i_speeds = 1:Nb_speed_lvl_cons
                num_speed_lvl_cons = speed_set(i_levels,1).num_speeds_glob(i_speeds,1);
                I_int_tot(num_speed_lvl_cons,:) = I_interpolated(i_speeds,:);  
            end              
        end
        
    end
      
    % Calcul du terme de scattering
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    S_scattered = Weights_quad_aniso*I_int_tot;
    
    % Calcul pour les niveaux qui le nécessite
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    Nb_levels_to_compute = size(Levels_to_compute,1);
    for i_levels = 1:Nb_levels_to_compute
        
        % Caractéristique du niveau considéré
        Num_level_cons = Levels_to_compute(i_levels,1);
        Nb_speed_level_cons = speed_set(Num_level_cons,1).Nb_speed;
        Num_speeds_glob_cons = speed_set(Num_level_cons,1).num_speeds_glob;
        
        % Etape de collision 
        Delta_t_cons = time_steps(Num_level_cons,1);
        I_tilde = I(Num_level_cons,1).values - Delta_t_cons*c*( (S_t').*I(Num_level_cons,1).values - (S_s').*S_scattered(Num_speeds_glob_cons,:) );
        
        % Etape de streaming
        I_streamed = zeros(Nb_speed_level_cons,Nb_nodes);
        for i_speeds = 1:Nb_speed_level_cons 
            Matchable_dir_cons = Matchable(Num_level_cons,1).speeds(i_speeds,1).num;
            Num_nodes_trans_dir_cons = Lattice(Num_level_cons,1).links(i_speeds,Matchable_dir_cons);
            I_streamed(i_speeds,Num_nodes_trans_dir_cons) = I_tilde(i_speeds,Matchable_dir_cons);
        end
        
        % Traitement des conditions limites qui deviennent ici "volumique
        Time_computed = Time_begin_unknown + Delta_t_cons;
        
        % Initialisation des champs inconnus
        if case_choice == 0
            
            % Pour la frontière 1
            if Num_level_cons == 1
                num_dir_glob_cons = speed_set(1,1).num_speeds_glob(1,1); 
                W_cons = Weights_quad_iso(1,num_dir_glob_cons);
                I_streamed(1,border_1.Num_nodes) = 0;
                I_streamed(1,Num_nodes_bd_left) = (1/W_cons)*I_0;
                for i_speeds = 2:Nb_speed_level_cons 
                    if border_1.Ps_dir_ode(Num_level_cons,1).values(i_speeds,1) > 0
                        I_streamed(i_speeds,border_1.Num_nodes) = 0;
                    end
                end 
            else
                for i_speeds = 1:Nb_speed_level_cons 
                    if border_1.Ps_dir_ode(Num_level_cons,1).values(i_speeds,1) > 0
                        I_streamed(i_speeds,border_1.Num_nodes) = 0;
                    end
                end
            end

            % Pour la frontière 2
            for i_speeds = 1:Nb_speed_level_cons 
                if border_2.Ps_dir_ode(Num_level_cons,1).values(i_speeds,1) > 0
                    I_streamed(i_speeds,border_2.Num_nodes) = 0;
                end
            end  

            % Pour la frontière 3
            for i_speeds = 1:Nb_speed_level_cons 
                if border_3.Ps_dir_ode(Num_level_cons,1).values(i_speeds,1) > 0
                    I_streamed(i_speeds,border_3.Num_nodes) = 0;
                end
            end 

            % Pour la frontière 4
            for i_speeds = 1:Nb_speed_level_cons 
                if border_4.Ps_dir_ode(Num_level_cons,1).values(i_speeds,1) > 0
                    I_streamed(i_speeds,border_4.Num_nodes) = 0;
                end
            end

        elseif case_choice == 1
            
            % Pour la frontière 1
            if Num_level_cons == 2
                num_dir_glob_cons = speed_set(2,1).num_speeds_glob(4,1); 
                W_cons = Weights_quad_iso(1,num_dir_glob_cons);
                I_streamed(4,border_1.Num_nodes) = 0;
                I_streamed(4,Num_nodes_bd_left) = (1/W_cons)*I_0;
                for i_speeds = 1:3 
                    if border_1.Ps_dir_ode(Num_level_cons,1).values(i_speeds,1) > 0
                        I_streamed(i_speeds,border_1.Num_nodes) = 0;
                    end
                end 
            else
                for i_speeds = 1:Nb_speed_level_cons 
                    if border_1.Ps_dir_ode(Num_level_cons,1).values(i_speeds,1) > 0
                        I_streamed(i_speeds,border_1.Num_nodes) = 0;
                    end
                end
            end

            % Pour la frontière 2
            for i_speeds = 1:Nb_speed_level_cons 
                if border_2.Ps_dir_ode(Num_level_cons,1).values(i_speeds,1) > 0
                    I_streamed(i_speeds,border_2.Num_nodes) = 0;
                end
            end  

            % Pour la frontière 3
            for i_speeds = 1:Nb_speed_level_cons 
                if border_3.Ps_dir_ode(Num_level_cons,1).values(i_speeds,1) > 0
                    I_streamed(i_speeds,border_3.Num_nodes) = 0;
                end
            end 

            % Pour la frontière 4
            if Num_level_cons == 2
                num_dir_glob_cons = speed_set(2,1).num_speeds_glob(4,1); 
                W_cons = Weights_quad_iso(1,num_dir_glob_cons);
                I_streamed(4,border_4.Num_nodes) = 0;
                I_streamed(4,Num_nodes_bd_top) = (1/W_cons)*I_0;
                for i_speeds = 1:3 
                    if border_4.Ps_dir_ode(Num_level_cons,1).values(i_speeds,1) > 0
                        I_streamed(i_speeds,border_4.Num_nodes) = 0;
                    end
                end    
            else
                for i_speeds = 1:Nb_speed_level_cons 
                    if border_4.Ps_dir_ode(Num_level_cons,1).values(i_speeds,1) > 0
                        I_streamed(i_speeds,border_4.Num_nodes) = 0;
                    end
                end
            end

        end
        
        
        %%%%%%%%%%
        % Mise à jour des données
        Elems_unknown_levels(Num_level_cons,1) = Elems_unknown_levels(Num_level_cons,1) + 1; 
        I_old(Num_level_cons,1).values = I(Num_level_cons,1).values;
        I(Num_level_cons,1).values = I_streamed; 
        
    end
    
    %%%%%%
    % test sur l'énergie
    E_field = Weights_quad_iso(1,:)*I_int_tot;
    %%%%%%
    
    % Vérification de l'état d'avancement
    state_levels = check_state_speed_levels(Nb_speed_level,sub_time_meshes,Elems_unknown_levels);
    
    
    % Calcul des grandeurs macro
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Energie totale:
    E_meas_tot = Weights_quad_iso(1,:)*I_int_tot;
    V_grid_Nodes = zeros(Nb_nodes_y,Nb_nodes_x);
    for i_nodes_grid = 1:Nb_nodes
        i_nd_cons = Num_mat_Nodes(i_nodes_grid,1);
        j_nd_cons = Num_mat_Nodes(i_nodes_grid,2);
        V_grid_Nodes(i_nd_cons,j_nd_cons) = E_meas_tot(1,i_nodes_grid);
    end
    F = griddedInterpolant(X_grid_Nodes,Y_grid_Nodes,V_grid_Nodes);
    Vq = F(r_Nodes_meas(:,1),r_Nodes_meas(:,2));
    E_meas(:,i_time_step) = Vq ;

    % Output pour vérification
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
  
    time_instant(i_time_step,1) = Time_begin_unknown;
    
    if case_choice == 0
        t_1 = floor( (1/2)*(time_mesh.Nb_nodes) ) ;
        t_2 = t_1 + 1;
        t_3 = t_1 + 2;
        if (i_time_step == t_1) || (i_time_step == t_2) || (i_time_step == t_3)
            fig = figure;
            hold on
            scatter(r_Nodes(:,1),r_Nodes(:,2),20,E_field,'filled');
            axis equal
            colormap(jet);
            col = colorbar;
            col.Label.Interpreter = 'latex';
            col.Label.String = "E";
            xlb = xlabel('$x$','Interpreter','latex');
            ylb = ylabel('$y$','Interpreter','latex');
            set(gcf,'Position',[200 100 900 600]);
            set(gca,'color','w'); 
            set(gcf,'color','w'); 
            set(gca,'FontSize',18);
            xlim([0,1]);
            ylim([0,1]);
            hold off
            folder = "./Figures/Laser/";
            name_fig = strcat(folder,"E_tot_",num2str(case_choice),"_",num2str(i_time_step),".eps");
            saveas(fig,name_fig);
            
            close all
            
            fig = figure;
            hold on
            scatter(r_Nodes(:,1),r_Nodes(:,2),200,E_field,'filled');
            axis equal
            colormap(jet);
            col = colorbar;
            col.Label.Interpreter = 'latex';
            col.Label.String = "E";
            xlb = xlabel('$x$','Interpreter','latex');
            ylb = ylabel('$y$','Interpreter','latex');
            set(gcf,'Position',[200 100 900 600]);
            set(gca,'color','w'); 
            set(gcf,'color','w'); 
            set(gca,'FontSize',18);
            xlim([0.49,0.51]);
            ylim([0.74,0.76]);
            hold off
            folder = "./Figures/Laser/";
            name_fig = strcat(folder,"E_tot_",num2str(case_choice),"_",num2str(i_time_step),"_zoom.eps");
            saveas(fig,name_fig);
            
            close all

        end
    else 
        
        t_1 = floor( (1/2)*(time_mesh.Nb_nodes) ) - 1;
        t_2 = time_mesh.Nb_nodes-1;
        if (i_time_step == t_1) || (i_time_step == t_2)
            fig = figure;
            hold on
            scatter(r_Nodes(:,1),r_Nodes(:,2),20,E_field,'filled');
            x_t = [0.0,c*(1/sqrt(2))*Time_begin_unknown,0.1+c*(1/sqrt(2))*Time_begin_unknown,0.1];
            y_t = [0.9,0.9-c*(1/sqrt(2))*Time_begin_unknown,1.0-c*(1/sqrt(2))*Time_begin_unknown,1.0];
            plot(x_t,y_t,'LineWidth',2,'Color','w');
            axis equal
            colormap(jet);
            col = colorbar;
            col.Label.Interpreter = 'latex';
            col.Label.String = "E";
            xlb = xlabel('$x$','Interpreter','latex');
            ylb = ylabel('$y$','Interpreter','latex');
            set(gcf,'Position',[200 100 900 600]);
            set(gca,'color','w'); 
            set(gcf,'color','w'); 
            set(gca,'FontSize',18);
            xlim([0,1]);
            ylim([0,1]);
            hold off
            folder = "./Figures/Laser/";
            name_fig = strcat(folder,"E_tot_",num2str(case_choice),"_",num2str(i_time_step),".eps");
            saveas(fig,name_fig);
            close 
            
        end
        
    end
    
end   


end