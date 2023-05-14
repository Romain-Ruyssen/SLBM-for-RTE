function [Weights] = get_Weights_trapez_Half_Circle(dir_HC,Nb_dir,directions)
% Fonction qui sort les poids de la quadrature en utilisant la loi des
% trapèzes pour un scattering avec une fonction de phase de Henyey–Greenstein
% dir_HC est la direction caractérisant la demi cercle
% [1 , 0] : demi cercle vertical de droite
% [0 , 1] : demi cercle horizontal supérieur

% Initialisation de l'output
Weights = zeros(1,Nb_dir);

% Calcul des angles de chaque direction (de 0 à 2*pi)
angles = zeros(Nb_dir,1);
for i_angle = 1:Nb_dir
    theta = get_Angle(directions(i_angle,1),directions(i_angle,2));
    if theta < 0
        theta = theta + 2*pi;
    end
    angles(i_angle,1) = theta;
end

% Arrangement des angles dans l'ordre croissant pour simplifier
[angles_tot_sort,I_permut] = sort(angles);


% Boucle sur les éléments du maillage du cercle
for i_m = 1:Nb_dir
    if i_m ~= Nb_dir
        num_nodes = [i_m , i_m+1];
        theta_nodes = angles_tot_sort(num_nodes,1);
        Delta_theta = theta_nodes(2,1) - theta_nodes(1,1); 
    else
        num_nodes = [i_m , 1];
        theta_nodes = angles_tot_sort(num_nodes,1);
        Delta_theta = 2*pi - theta_nodes(1,1) + theta_nodes(2,1); 
    end
    scalar_prods = dir_HC*[ [cos(theta_nodes(1,1)) , cos(theta_nodes(2,1))] ;
                            [sin(theta_nodes(1,1)) , sin(theta_nodes(2,1))] ];             
    if ( (scalar_prods(1,1) > -10^(-10)) && (scalar_prods(1,2) > -10^(-10)) )
        % Les deux ps positifs ou nuls c'est qu'on est du bon coté du
        % cercle
        num_nodes_real = I_permut(num_nodes);
        Weights(1,num_nodes_real) = Weights(1,num_nodes_real) + Delta_theta/2;
    end
    
end

end
