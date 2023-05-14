function state_levels = check_state_speed_levels(Nb_speed_level,sub_time_meshes,Elems_unknown_levels)
% Fonction pour vérifier si tous les éléments de chaque niveau de vitesse
% sont connus.

state_levels = true;

check_all = zeros(Nb_speed_level,1);
for i_levels = 1:Nb_speed_level
    Nb_elems_level_cons = sub_time_meshes(i_levels,1).Nb_elems;
    if Elems_unknown_levels(i_levels,1) >= (Nb_elems_level_cons+1)
        check_all(i_levels,1) = 1;
    else
        break
    end
end

if sum(check_all,1) == Nb_speed_level
    state_levels = false;
end

end