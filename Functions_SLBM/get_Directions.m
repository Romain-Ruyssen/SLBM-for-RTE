function [directions] = get_Directions(Nb_speed_level,speed_set)
% Une fonction simplement pour concaténer toutes les directions
% dans une seule matrice

directions = [];
for i_dir = 1:Nb_speed_level
    directions = [directions ; speed_set(i_dir,1).speeds];
end

end