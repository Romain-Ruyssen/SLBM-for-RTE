function Lattice = get_lattice_New(Nb_speed_level,speed_set,Nb_nodes,Nb_nodes_x,Nb_nodes_y,sym_left_right,sym_bot_top)

% Fonction pour construire, en fonction d'un niveau souhaité, les réseaux.

for i_levels = 1:Nb_speed_level
    Nb_speeds_cons = speed_set(i_levels,1).Nb_speed;
    Lattice(i_levels,1).links = zeros(Nb_speeds_cons,Nb_nodes);
    for i_nodes = 1:Nb_nodes_x
        for j_nodes = 1:Nb_nodes_y
            Num_node_cons = (j_nodes-1)*Nb_nodes_x + i_nodes;
            for i_speeds = 1:Nb_speeds_cons
                speed_cons = speed_set(i_levels,1).speeds(i_speeds,:);
                i_node_link = i_nodes + speed_cons(1,1);
                j_node_link = j_nodes + speed_cons(1,2);
                
                % En fonction de si on veut une symétrie ou non
                Num_node_link = 1;
                
                if (i_node_link <= 0)
                    if sym_left_right
                        i_node_link = i_node_link + Nb_nodes_x;
                    else
                        Num_node_link = 0;
                    end  
                elseif (i_node_link >= Nb_nodes_x+1)
                    if sym_left_right
                        i_node_link = i_node_link - Nb_nodes_x;
                    else
                        Num_node_link = 0;
                    end 
                end
                
                if (j_node_link <= 0)
                    if sym_bot_top
                        j_node_link = j_node_link + Nb_nodes_y;
                    else
                        Num_node_link = 0;
                    end 
                elseif (j_node_link >= Nb_nodes_y+1)
                    if sym_bot_top
                        j_node_link = j_node_link - Nb_nodes_y;
                    else
                        Num_node_link = 0;
                    end  
                end

                if Num_node_link ~= 0 % Si après tous les tests il n'est toujours pas nul
                    Num_node_link = (j_node_link-1)*Nb_nodes_x + i_node_link;
                end
                
                Lattice(i_levels,1).links(i_speeds,Num_node_cons) = Num_node_link;
            end
        end    
    end    
end


end