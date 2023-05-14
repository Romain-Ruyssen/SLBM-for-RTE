function [time_mesh,E_meas,F_paper,time_instant] = ... 
              solve_HetAni(scale_Dx,choice_case,Nb_dir_tot, ...
                Nb_speed_level,speed_set, ...
                Weights_quad_aniso,Weights_quad_iso, ...
                Nb_pt_meas,r_Nodes_meas)
% Fonction pour résoudre le problème hétérogène avec scattering anisotrope
% en utilisant l'approche MTL
% choice_case: 0 - milieu homogène / 1 - milieu hétérogène

% Définition des paramètres:
c = 1;
L_x = 5; % (m) length of the slab
L_y = 5;

T_period = 3*(L_x)/c;

% Discrétisation du domaine spatial:
x_begin = -2.5;
y_begin = -2.5;

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
S_a = 0.01;
S_s_out = 1; % valeur hors du disque
S_s_in  = 5;  % valeur dans le disque
S_s = S_s_out*ones(Nb_nodes,1); % scattering cross section

if choice_case == 1 % Cas hétérogène
    x_c = 1;
    y_c = 1.5;
    hete_rad = 0.3;
    rad_Nodes = sqrt((r_Nodes(:,1) - x_c).^2 + (r_Nodes(:,2) - y_c).^2);
    Num_nodes_in = find(rad_Nodes <= hete_rad + e);
    S_s(Num_nodes_in,1) = S_s_in; % scattering cross section
end
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

% Initialisation des champs d'énergie avec une strucutre en niveau
for i_levels = 1:Nb_speed_level
    Nb_dir = speed_set(i_levels,1).Nb_speed;
    I(i_levels,1).values = zeros(Nb_dir,Nb_nodes);
    I_old(i_levels,1).values = zeros(Nb_dir,Nb_nodes);
end

% Recherche des noeuds de la frontière où l'on impose le flux du laser
% f_t = zeros(Nb_time_step+1,1); % Fonction temporelle
t_0 = 1; % Attention ?
w_t = 1/sqrt(2);
% for i_time = 1:(Nb_time_step+1)
%     time_cons = (i_time-1)*Delta_t;
%     f_t(i_time,1) = (2/pi)*exp( -((time_cons - t_0)/w_t)^2 );
% end
f_t = @(var_time) (2/pi)*exp( -((var_time - t_0)/w_t).^2 );

f_x = zeros(border_1.Nb_nodes,1);
w_s = 1/sqrt(10);
x_s = -2.5;
y_s = 1.5;
for i_nodes = 1:border_1.Nb_nodes
    num_node_cons = border_1.Num_nodes(i_nodes,1);
    y_cons = r_Nodes(num_node_cons,2);
    f_x(i_nodes,1) = exp(-((y_cons - y_s)/w_s)^2);
end
f_x = f_x/Weights_quad_iso(1,1);

% Initialisation des champs inconnus
I(1,1).values(1,border_1.Num_nodes) = f_x*f_t(0); % f_x*f_t(1,1);
I_old(1,1).values(1,border_1.Num_nodes) = f_x*f_t(0); % f_x*f_t(1,1);

% Définition du point où l'on mesure le backscattered flux
Num_Nodes_bf = zeros(Nb_pt_meas,1);
for i_nodes = 1:Nb_pt_meas
    dist_nodes_bsf = sqrt( (r_Nodes(:,1) - r_Nodes_meas(i_nodes,1)).^2 + (r_Nodes(:,2) - r_Nodes_meas(i_nodes,2)).^2 );
    [~,P_bf] = min(dist_nodes_bsf);
    Num_Nodes_bf(i_nodes,1) = P_bf;
end

% Initialisation des données que nous relevons aux points de mesure
E_meas = zeros(Nb_pt_meas,Nb_time_step);
F_paper = zeros(Nb_pt_meas,Nb_time_step);

% Pour calculer facilement les fluxs on définit les directions
directions_for_flux = [];
for i_levels = 1:Nb_speed_level        
    directions_for_flux = [directions_for_flux ; 
       (1/speed_set(i_levels,1).norm)*speed_set(i_levels,1).speeds];
end

% Pour faire les mesures d'énergies
E_injected = zeros(Nb_time_step,1);
% E_injected(1,1) = sum(Weights_quad_iso(1,1)*f_x*f_t(0),1);

E_out.left = zeros(Nb_time_step,1);
E_out.bot = zeros(Nb_time_step,1);
E_out.right = zeros(Nb_time_step,1);
E_out.top = zeros(Nb_time_step,1);

time_instant = zeros(Nb_time_step,1);

for i_levels = 1:Nb_speed_level
    Num_dir_out_left_tot{i_levels,1}.values = find(border_1.Ps_dir_ode(i_levels,1).values < 0);
    Num_dir_out_bot_tot{i_levels,1}.values = find(border_2.Ps_dir_ode(i_levels,1).values < 0);
    Num_dir_out_right_tot{i_levels,1}.values = find(border_3.Ps_dir_ode(i_levels,1).values < 0);
    Num_dir_out_top_tot{i_levels,1}.values = find(border_4.Ps_dir_ode(i_levels,1).values < 0);
end
    
E_tot = zeros(Nb_time_step,1);
% E_tot(1,1) = E_injected(1,1); %Weights_quad_iso(1,:)*I;

% Ce que j'ai ajouté juste pour les 3 images
instant_images = [2.0 ; 4.0 ; 6.0];
nb_times = 3;
[~,t_1] = min(abs(time_mesh.Nodes - instant_images(1,1))); 
[~,t_2] = min(abs(time_mesh.Nodes - instant_images(2,1)));    
[~,t_3] = min(abs(time_mesh.Nodes - instant_images(3,1)));

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
        
        % Pour la frontière 1
        if Num_level_cons == 1
            I_streamed(1,border_1.Num_nodes) = f_x*f_t(Time_computed);
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
        
        %%%%%%
        % test sur l'énergie
        if Num_level_cons == 1
            E_injected(i_time_step,1) = E_injected(i_time_step,1) + ...
                                    Delta_t_cons*sum(Weights_quad_iso(1,1)*f_x*f_t(Time_begin_unknown),1);
        end
%         E_field = Weights_quad_iso(1,Num_speeds_glob_cons)*I_tilde;
%         E_tot(i_time_step,1) = E_tot(i_time_step,1) + sum(E_field);
    
        Num_dir_out_left = Num_dir_out_left_tot{Num_level_cons,1}.values;
        Num_dir_out_bot = Num_dir_out_bot_tot{Num_level_cons,1}.values; 
        Num_dir_out_right = Num_dir_out_right_tot{Num_level_cons,1}.values;
        Num_dir_out_top = Num_dir_out_top_tot{Num_level_cons,1}.values;
        E_out.left(i_time_step,1) = E_out.left(i_time_step,1) + Delta_t_cons*sum(Weights_quad_iso(1,Num_dir_out_left)*((directions_for_flux(Num_dir_out_left,:)*(-border_1.Normal_in)).*I_tilde(Num_dir_out_left,border_1.Num_nodes)),2);
        E_out.bot(i_time_step,1) = E_out.bot(i_time_step,1) + Delta_t_cons*sum(Weights_quad_iso(1,Num_dir_out_bot)*((directions_for_flux(Num_dir_out_bot,:)*(-border_2.Normal_in)).*I_tilde(Num_dir_out_bot,border_2.Num_nodes)),2);
        E_out.right(i_time_step,1) = E_out.right(i_time_step,1) + Delta_t_cons*sum(Weights_quad_iso(1,Num_dir_out_right)*((directions_for_flux(Num_dir_out_right,:)*(-border_3.Normal_in)).*I_tilde(Num_dir_out_right,border_3.Num_nodes)),2);
        E_out.top(i_time_step,1) = E_out.top(i_time_step,1) + Delta_t_cons*sum(Weights_quad_iso(1,Num_dir_out_top)*((directions_for_flux(Num_dir_out_top,:)*(-border_4.Normal_in)).*I_tilde(Num_dir_out_top,border_4.Num_nodes)),2);
        %%%%%%
        
        %%%%%%%%%%
        % Mise à jour des données
        Elems_unknown_levels(Num_level_cons,1) = Elems_unknown_levels(Num_level_cons,1) + 1; 
        I_old(Num_level_cons,1).values = I(Num_level_cons,1).values;
        I(Num_level_cons,1).values = I_streamed; 
        
    end
    
    %%%%%%
    % test sur l'énergie
    E_field = Weights_quad_iso(1,:)*I_int_tot;
    E_tot(i_time_step,1) = sum(E_field);
    %%%%%%
    
    % Vérification de l'état d'avancement
    state_levels = check_state_speed_levels(Nb_speed_level,sub_time_meshes,Elems_unknown_levels);
    
    
    % Calcul des grandeurs macro
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Energie totale:
    E_meas(:,i_time_step) = Weights_quad_iso(1,:)*I_int_tot(:,Num_Nodes_bf);
    
    % Flux:
    F_paper(:,i_time_step) = (Weights_quad_iso(1,:).*(abs(directions_for_flux(:,1)')))*I_int_tot(:,Num_Nodes_bf);
    
    
end


end
