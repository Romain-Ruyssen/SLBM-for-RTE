function sub_time_meshes = get_sub_time_meshes(Nb_speed_level,time_mesh)
% Fonction pour obtenir les sous maillages temporels de chaque niveau de
% vitesse à partir du maillage global


    for i_levels = 1:Nb_speed_level
        
        [~,C] = find(time_mesh.who_levels(i_levels,:));
        
        sub_time_meshes(i_levels,1).Nb_nodes = size(C,2);                                   % Nombre de noeuds
        sub_time_meshes(i_levels,1).Num_node = C;                                           % Numérotation globale de ses noeuds
        sub_time_meshes(i_levels,1).Elems = [ C(1:1:end-1)' , C(2:1:end)' ];                % Eléments avec la numérotation globale des noeuds
        sub_time_meshes(i_levels,1).Nb_elems = size(sub_time_meshes(i_levels,1).Elems,1);   % Nombre d'éléments
    
    end


end