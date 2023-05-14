function time_mesh = construct_time_mesh(Period,time_steps)
% Fonction pour construire le maillage temporel

% ATTENTION: il faut que les pas de temps soient rangés dans l'ordre
% croissant.

e = 10^(-10); % Précision numérique

Nb_times_levels = size(time_steps,1); % Nombre de couche temporelle


% Initialisation du maillage
x_Nodes = [];
who_time = ones(Nb_times_levels,1);

% Boucle sur tous les instants (ensemble formé de tous les instants de tous
% les niveaux temporels)
x_Nodes = [x_Nodes ; 0];
counters = ones(Nb_times_levels,1);
keep_going = 1;
while keep_going
    
    times_cons = counters.*time_steps; % Le temps maximal de chaque niveau temporel
    R = find(times_cons <= min(times_cons) + e); % On trouve les emplacements du temp minimum 
                                                 % et en même temps les niveaux de temps égaux
    if min(times_cons) <= Period                                              
        x_Nodes = [x_Nodes ; min(times_cons)];
        who_time_cons = zeros(Nb_times_levels,1);
        who_time_cons(R,1) = 1;
        who_time = [who_time , who_time_cons];
        counters(R,1) = counters(R,1) + 1; % Incrémentation des compteurs nécessaires
    else
        keep_going = 0;
    end
    
end

% construction du maillage temporel
time_mesh.Nodes = x_Nodes; % Coordonnées des noeuds
time_mesh.Nb_nodes = size(time_mesh.Nodes,1); % Nombre de noeuds
time_mesh.Elems = [(1:1:time_mesh.Nb_nodes-1)' , (2:1:time_mesh.Nb_nodes)'];
time_mesh.Nb_elems = size(time_mesh.Elems,1);
time_mesh.Nb_times_level = Nb_times_levels;
time_mesh.times_steps = time_steps;
time_mesh.who_levels = who_time;


end