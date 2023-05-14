function angle = get_Angle(coor_x,coor_y)
% Fonction qui retourne l'angle polaire compris entre ]-pi,pi].

e = 10^(-10);
norm = sqrt( coor_x^2 + coor_y^2 ); 

if coor_y > e
    angle = acos(coor_x/norm);
elseif coor_y < e
    angle = -acos(coor_x/norm);
else % Dans ce cas coor_y est nul
    if coor_x > e
        angle = 0;
    elseif coor_x < e
        angle = pi; % On aurait pu mettre -pi
    else
        angle = 0; % On va mettre ca comme base
    end
end


end