function Matchable = get_Matchable(Lattice,speed_set,Nb_speed_level)
% Fonction pour obtenir les noeuds appariables

    for i_levels = 1:Nb_speed_level
        Nb_speeds_cons = speed_set(i_levels,1).Nb_speed;
        for i_speeds = 1:Nb_speeds_cons
            Matchable(i_levels,1).speeds(i_speeds,1).num = find(Lattice(i_levels,1).links(i_speeds,:)); % Appariable
        end 
    end

end