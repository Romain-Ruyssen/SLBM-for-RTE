function [speed_set,Nb_dir_tot] = construct_speed_set(Nb_speed_level)
% Fonction qui construit l'ensemble des niveaux de vitesse à partir d'un
% nombre de niveau donné

e = 10^(-10); % Précision numérique

% En gros on construit au fur et à mesure tel que la norme discrète l1 
% (petit l et pas grand L) de couple de N^2 (entiers naturels)

Ray_max = 30; % Lui quelque part nous donne le Dt_max/Dt_min maximal 
Size_Mesh_grid_max = floor(Ray_max) + 1;

% Construction des vecteurs liens:
[X,Y] = meshgrid( (-Size_Mesh_grid_max:1:Size_Mesh_grid_max) , (-Size_Mesh_grid_max:1:Size_Mesh_grid_max) );
xy_Nodes = [X(:) Y(:)];
check_1 = find(xy_Nodes(:,1) == 0);
check_2 = find(xy_Nodes(check_1,2) == 0);
xy_Nodes(check_1(check_2),:) = []; % On retire [0 0]

Nb_links = size(xy_Nodes,1);

% Calcul des normes des liens:
Norm_links = sqrt(xy_Nodes(:,1).^2 + xy_Nodes(:,2).^2);
Num_links = (1:Nb_links)';

% Rangement des directions par normes
Nb_dir_tot = 0;
norm_cons = 0;
Num_lvl_cons = 0;
dir_already_cons = [];
while (norm_cons < Ray_max) && (Num_lvl_cons < Nb_speed_level) 
    
    norm_cons = min(Norm_links);
    
    [R,~] = find( (Norm_links < norm_cons + e) & (Norm_links > norm_cons - e) );
    
    nb_speeds = size(R,1);
    num_links_cons = Num_links(R,1);
    dir_links_cons = xy_Nodes(num_links_cons,:);
    
    % Rangement dans l'ordre croissant des angles
    % Calcul des angles de chaque direction (de 0 à 2*pi)
    angles = zeros(nb_speeds,1);
    for i_angle = 1:nb_speeds
        theta = get_Angle(dir_links_cons(i_angle,1),dir_links_cons(i_angle,2));
        if theta < 0
            theta = theta + 2*pi;
        end
        angles(i_angle,1) = theta;
    end

    % Arrangement des angles dans l'ordre croissant pour simplifier
    [~,I_permut] = sort(angles);
    dir_links_cons = dir_links_cons(I_permut,:);
    
    % Vérification que les directions à ajouter ne sont pas déjà
    % considérées
    
    if isempty(dir_already_cons) % Si l'ensemble total est vide on ajoute
        
        check_isempty_dir_add = 0;
        nb_speeds_add = nb_speeds;
        dir_links_cons_add = dir_links_cons;
        
    else
        
        dir_kept = [];
        for i_verif = 1:nb_speeds
            
            vect_prod_all = abs(dir_links_cons(i_verif,1)*dir_already_cons(:,2) ...
                          - dir_links_cons(i_verif,2)*dir_already_cons(:,1));
                      
            if min(vect_prod_all) > e % Il n'y a pas de zéro 
                dir_kept = [dir_kept ; i_verif];
            end

        end
        
        if ~isempty(dir_kept)
            check_isempty_dir_add = 0;    
            
             nb_speeds_add = size(dir_kept,1);
            dir_links_cons_add = dir_links_cons(dir_kept,:);
               
        else
            check_isempty_dir_add = 1;
        end
        
    end
    
    if ~check_isempty_dir_add % Si l'ensemble des directions à ajouter n'est pas vide
        
        Num_lvl_cons = Num_lvl_cons + 1;
        
        speed_set(Num_lvl_cons,1).Nb_speed = nb_speeds_add;
        speed_set(Num_lvl_cons,1).speeds = dir_links_cons_add;
        speed_set(Num_lvl_cons,1).norm = norm_cons;

        % Construction d'une numérotation globale:
        speed_set(Num_lvl_cons,1).num_speeds_glob = (Nb_dir_tot+1:1:Nb_dir_tot+speed_set(Num_lvl_cons,1).Nb_speed)';

        % Incrémentation du nombre de vitesse:
        Nb_dir_tot = Nb_dir_tot + speed_set(Num_lvl_cons,1).Nb_speed;
        
        % Ajout des directions à l'ensemble global
        dir_already_cons = [dir_already_cons ; dir_links_cons_add];
        
    end
    
    % Suppession des directions considérées
    Norm_links(R) = [];
    Num_links(R) = [];
    
    
end

end