function [Int_operators] = get_Int_operators(Nb_nodes,speed_off_grid,Interpolation,Matchable)
% Fonction pour construire pour chaque direction qui tombe off grid un
% opérateur d'interpolation

Nb_speed_off = size(speed_off_grid,1);
    
    for i_speed_off = 1:Nb_speed_off
        
        dir_cons = speed_off_grid(i_speed_off,1);
        Nb_nodes_int = size(Matchable.off(i_speed_off,1).num,1);
%         Int_operators(i_speed_off,1).op = zeros(Nb_nodes_int,Nb_nodes);
        Int_operators(i_speed_off,1).op = sparse(Nb_nodes_int,Nb_nodes);
        
        for i_nodes = 1:size(Matchable.off(i_speed_off,1).num,1)
 
            Num_node_cons = Matchable.off(i_speed_off,1).num(i_nodes,1);
            num_nodes_int = Interpolation.directions(i_speed_off,Num_node_cons).num;
            weight_int = Interpolation.directions(i_speed_off,Num_node_cons).weight;
            Int_operators(i_speed_off,1).op(i_nodes,num_nodes_int) = weight_int;
               
        end
        
    end  

end