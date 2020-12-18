% Subroutine that generates a discretized circle in finite elements
%   -Each element has a linear geometry and constant value over their
%   domain
% Author : Edward H-Hannan
% Date of Creation : 2019-07-15
% Last update : 2020-12-17

function [cont, elem, nb_cont, nb_elem] = generate_circle(cont, elem, nb_cont, nb_elem, ordering, r, origin, nel_circle, type, val)

    if ordering == 0 %0=clockwise 
        angles_geo = flip(linspace(0,2*pi,nel_circle+1));
    elseif ordering == 1 %1=anticlockwise
        angles_geo = linspace(0,2*pi,nel_circle+1);
    end

    for i = 1 : nel_circle
        index = nb_elem + i; 
        elem(index).p1 = [r*cos(angles_geo(i))+origin(1) , r*sin(angles_geo(i))+origin(2)]; 
        elem(index).p2 = [r*cos(angles_geo(i+1))+origin(1) , r*sin(angles_geo(i+1))+origin(2)];
        elem(index).l = elem(index).p2 - elem(index).p1; %array of elements belonging to the contour
        ln = elem(index).l/norm(elem(index).l); %normalized vector
        elem(index).ln = ln;
        n = cross([ln 0],[0 0 1]); %normal vector
        elem(index).n = n(1:2); 
        elem(index).type = type; % 0=Neumann, 1=Dirichlet
        elem(index).val = val;
    end
    nb_cont = nb_cont + 1; %number of contour
    cont(nb_cont).nb_elem = nel_circle; %number of element in the contour
    cont(nb_cont).elem = [nb_elem+1 : nb_elem+nel_circle]; %array of elements belonging to the contour
    cont(nb_cont).rotation = ordering; %0=clockwise, 1=anticlockwise
    nb_elem=nb_elem + nel_circle;

end