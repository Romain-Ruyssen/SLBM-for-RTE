function [Weights,angles] = get_Weights_trapez_HG(Nb_dir,directions,g)
% Fonction qui sort les poids de la quadrature en utilisant la loi des
% trapèzes pour un scattering avec une fonction de phase de Henyey–Greenstein

% Initialisation de l'output
Weights = zeros(Nb_dir);

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

% Calcul des poids: on parcourt les angles dans leur organisation
% croissante
for i_m = 1:Nb_dir
    thetha_m = angles_tot_sort(i_m);
    i_m_perm = I_permut(i_m,1); % Numéro dans la numérotation initiale de direction
    for i_mp = 1:Nb_dir
        thetha_mp = angles_tot_sort(i_mp);
        i_mp_perm = I_permut(i_mp,1);
        if i_mp == 1
            theta_mp_m1 = angles_tot_sort(end,1) - 2*pi;
            theta_mp_p1 = angles_tot_sort(2,1);
        elseif i_mp == Nb_dir
            theta_mp_m1 = angles_tot_sort(Nb_dir-1,1);
            theta_mp_p1 = 2*pi + angles_tot_sort(1,1);
        else
            theta_mp_m1 = angles_tot_sort(i_mp-1,1);
            theta_mp_p1 = angles_tot_sort(i_mp+1,1);
        end
        Delta_theta_m1 = rem(abs(thetha_mp - theta_mp_m1),pi);
        Delta_theta_p1 = rem(abs(thetha_mp - theta_mp_p1),pi);
        Delta_theta_tot = (Delta_theta_m1/2) + (Delta_theta_p1/2);
        theta_int_a = (theta_mp_m1 + thetha_mp)/2;
        theta_int_b = (thetha_mp + theta_mp_p1)/2;
        f_HG_cons = @(theta_p) (1/(2*pi))*((1-g^2)./(1+g^2-2*g*cos(thetha_m-theta_p)));
        Weights(i_m_perm,i_mp_perm) = Delta_theta_tot*(0.5*f_HG_cons(theta_int_a) + 0.5*f_HG_cons(theta_int_b));
    end
end

end
