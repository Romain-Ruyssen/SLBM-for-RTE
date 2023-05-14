function [ time_mesh, q_T, q_R] = ...
     solve_Pulse(scale_Dx,Nb_dir_tot,Nb_speed_level,speed_set,Weights_quad,Weights_HCR,Weights_HCL,omega)
% Fonction pour résoudre le problème du pulse avec l'approche
% Multi Time Levels
 

% Définition des paramètres:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

L_x = 1; % (m) length of the slab
L_y = 0.2;

c = 1;

S_t = 1;       % (m^(-1)) extinction coefficient
S_s = omega*S_t;
S_a = S_t - S_s;

T_period = 6*(L_x/c);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solution numérique

% Discrétisation du domaine spatial:

x_begin = 0;
y_begin = 0;

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

% Construction sous forme de maillage pour pouvoir regarder la conservation
Mesh.nb_nodes = Nb_nodes;
Mesh.nodes = r_Nodes;
Mesh.nb_elems = (Nb_nodes_x-1)*(Nb_nodes_y-1);
for j_nodes = 1:Nb_nodes_y-1
    for i_nodes = 1:Nb_nodes_x-1
        num_elem_cons = (j_nodes - 1)*(Nb_nodes_x-1) + i_nodes;
        num_node_cons = (j_nodes - 1)*Nb_nodes_x + i_nodes;
        Mesh.elems(num_elem_cons,:) = [ num_node_cons , num_node_cons + 1 , num_node_cons + Nb_nodes_x + 1 , num_node_cons + Nb_nodes_x];
    end
end

% Définition de la Lattice
Lattice = get_lattice(Nb_speed_level,speed_set,Nb_nodes,Nb_nodes_x,Nb_nodes_y,0,1);

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

%%%%%%%%%%%%%%%%%%%%%%%%
% Caractérisation des frontières volumiques nécessaires pour imposer les conditions de symétrie
% Nous allons pour chaque niveau de vitesse définir les correspondances
Delta_x_speed_lvl = zeros(Nb_speed_level,1);
for i_speed = 1:Nb_speed_level
    Delta_x_speed_lvl(i_speed,1) = speed_set(i_speed,1).norm*Delta_x;
end
bd_thickness = max(Delta_x_speed_lvl);

% bot border: y = 0
border_volumic_2.Num_nodes = find(r_Nodes(:,2) <= y_begin + bd_thickness + 10^(-10));


border_2.Num_nodes = find(r_Nodes(:,2) <= y_begin + 10^(-10));
border_2.Nb_nodes = size(border_2.Num_nodes,1);
border_2.Normal_in = [0 ; 1];
for i_levels = 1:Nb_speed_level
    border_2.Ps_dir_ode(i_levels,1).values = zeros(speed_set(i_levels,1).Nb_speed,1);
    for i_speeds = 1:speed_set(i_levels,1).Nb_speed
        border_2.Ps_dir_ode(i_levels,1).values(i_speeds,1) = (1/speed_set(i_levels,1).norm)*( (speed_set(i_levels,1).speeds(i_speeds,:)) * (border_2.Normal_in) );
    end
end

% top border: y = L_y

%%%%%%%%%%%%%%%%%%%%%%%%

% Initialisation des champs d'énergie avec une strucutre en niveau
for i_levels = 1:Nb_speed_level
    Nb_dir = speed_set(i_levels,1).Nb_speed;
    I(i_levels,1).values = zeros(Nb_dir,Nb_nodes);
    I_old(i_levels,1).values = zeros(Nb_dir,Nb_nodes);
end


% Résolution du problème:
%%%%%%%%%%%%%%%%%%%%%%%%%

% Définition des pulses d'intensité
I_p = ones(1,Nb_time_step);
% I_0 = (1/Weights_quad(1,1))*1;
I_0 = 1;
Nb_iter_time_stop = max(find(time_mesh.Nodes <= 1 + 10^(-10)));
for i_time = Nb_iter_time_stop:Nb_time_step
    I_p(1,i_time) = 0;
end

% Pour la condition limite
I(1,1).values(1,border_1.Num_nodes) = I_0;

% Initialisation de la Réflectance et de la Transmittance
q_R = zeros(1,Nb_time_step);
q_T = zeros(1,Nb_time_step);
check_1 = find(r_Nodes(:,1) <= 10^(-10));
check_2 = find((r_Nodes(check_1,2) <= L_y/2 + 10^(-10)) & (r_Nodes(check_1,2) >= L_y/2 - 10^(-10)));
Num_node_q_R = check_1(check_2);
clear check_1 check_2
check_1 = find(r_Nodes(:,1) >= L_x - 10^(-10));
check_2 = find((r_Nodes(check_1,2) <= L_y/2 + 10^(-10)) & (r_Nodes(check_1,2) >= L_y/2 - 10^(-10)));
Num_node_q_T = check_1(check_2);
clear check_1 check_2

% Initialisation des coefficients pour le calculr de la Réflectance et de la Transmittance
Ps_dir_ode_q_R = zeros(Nb_dir_tot,1);
Ps_dir_ode_q_T = zeros(Nb_dir_tot,1);
for i_levels = 1:Nb_speed_level
    for i_speeds = 1:speed_set(i_levels,1).Nb_speed
        num_speed_glob_cons = speed_set(i_levels,1).num_speeds_glob(i_speeds,1);
        if border_1.Ps_dir_ode(i_levels,1).values(i_speeds,1) < 0 - 10^(-10)
            Ps_dir_ode_q_R(num_speed_glob_cons,1) = abs(border_1.Ps_dir_ode(i_levels,1).values(i_speeds,1));
        elseif border_3.Ps_dir_ode(i_levels,1).values(i_speeds,1) < 0 - 10^(-10)
            Ps_dir_ode_q_T(num_speed_glob_cons,1) = abs(border_3.Ps_dir_ode(i_levels,1).values(i_speeds,1));              
        end
    end
end

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
    
    % Calcul du terme de scattering
    I_int_tot = zeros(Nb_dir_tot,Nb_nodes); % Intensité interpolée  
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
    E = Weights_quad*I_int_tot;
    
    % Calcul pour les niveaux qui le nécessite
    Nb_levels_to_compute = size(Levels_to_compute,1);
    for i_levels = 1:Nb_levels_to_compute
        
        % Caractéristique du niveau considéré
        Num_level_cons = Levels_to_compute(i_levels,1);
        Nb_speed_level_cons = speed_set(Num_level_cons,1).Nb_speed;
        Num_speeds_glob_cons = speed_set(Num_level_cons,1).num_speeds_glob;
        
        % Etape de collision 
        
        E_cons = E(Num_speeds_glob_cons,:);
        
        % Calcul des densités d'énergie collidées
        Delta_t_cons = time_steps(Num_level_cons,1);
        I_tilde = I(Num_level_cons,1).values - Delta_t_cons*c*(S_t*I(Num_level_cons,1).values - (S_s/(2*pi))*E_cons);
        
        % Etape de streaming
        I_streamed = zeros(Nb_speed_level_cons,Nb_nodes);
        for i_speeds = 1:Nb_speed_level_cons 
            Matchable_dir_cons = Matchable(Num_level_cons,1).speeds(i_speeds,1).num;
            Num_nodes_trans_dir_cons = Lattice(Num_level_cons,1).links(i_speeds,Matchable_dir_cons);
            I_streamed(i_speeds,Num_nodes_trans_dir_cons) = I_tilde(i_speeds,Matchable_dir_cons);
        end
        
        % Traitement des conditions limites

        % Pour la frontière 1
        if Num_level_cons == 1
            I_streamed(1,border_1.Num_nodes) = I_p(1,i_time_step);
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
        
        
        % Pour la frontière 3
        for i_speeds = 1:Nb_speed_level_cons 
            if border_3.Ps_dir_ode(Num_level_cons,1).values(i_speeds,1) > 0
                I_streamed(i_speeds,border_3.Num_nodes) = 0;
            end
        end  
   
        
        % Mise à jour des données
        Elems_unknown_levels(Num_level_cons,1) = Elems_unknown_levels(Num_level_cons,1) + 1; 
        I_old(Num_level_cons,1).values = I(Num_level_cons,1).values;
        I(Num_level_cons,1).values = I_streamed; 
        
    end
    
    % Calcul de la transmittance et de la reflectance
    q_R_cons = (1/I_0)*Weights_HCL* ...
                      (Ps_dir_ode_q_R.*I_int_tot(:,Num_node_q_R));
    q_T_cons = (1/I_0)*Weights_HCR* ...
                      (Ps_dir_ode_q_T.*I_int_tot(:,Num_node_q_T));
    
    q_T(1,i_time_step) = q_T_cons;
    q_R(1,i_time_step) = q_R_cons;
    
    % Vérification de l'état d'avancement
    state_levels = check_state_speed_levels(Nb_speed_level,sub_time_meshes,Elems_unknown_levels);
    
end

end