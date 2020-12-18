% Subroutine that generates a discretized squared in finite elements
%   -Each element has a linear geometry and constant value over their
%   domain
% Author : Edward H-Hannan
% Date of Creation : 2019-07-15
% Last update : 2020-12-17

function [cont, elem, nb_cont, nb_elem] = generate_rectangle(cont, elem, nb_cont, nb_elem, ordering, d1, d2, translation, nb_elem_rec, type_rec, val_rec)
   
    if ordering == 0
        c1_x = ones(1,nb_elem_rec(1)+1) * d2/2;          
        c1_y = flip(-d1/2 : d1/nb_elem_rec(1) : d1/2);
        
        c2_x = (d2/2 : -d2/nb_elem_rec(2) : -d2/2);      
        c2_y = ones(1, nb_elem_rec(2)+1) * -d1/2;
        
        c3_x = ones(1,nb_elem_rec(3)+1) * -d2/2;         
        c3_y = (-d1/2 : d1/nb_elem_rec(3) : d1/2); 
        
        c4_x = flip(d2/2 : -d2/nb_elem_rec(4) : -d2/2);  
        c4_y = ones(1,nb_elem_rec(4)+1) * d1/2;
    elseif ordering == 1
        c1_x = ones(1,nb_elem_rec(1)+1) * d2/2;          
        c1_y = -d1/2 : d1/nb_elem_rec(1) : d1/2;
        
        c2_x = (d2/2 : -d2/nb_elem_rec(2) : -d2/2);     
        c2_y = ones(1,nb_elem_rec(2)+1) * d1/2;
        
        c3_x = ones(1,nb_elem_rec(3)+1) * -d2/2;        
        c3_y = flip(-d1/2 : d1/nb_elem_rec(3) : d1/2); 
        
        c4_x = flip(d2/2 : -d2/nb_elem_rec(4) : -d2/2);  
        c4_y = ones(1, nb_elem_rec(4)+1) * -d1/2;
    end

    % Construction of the points of the rectangle 
    point_rec_x = [c1_x  , c2_x(2:end) , c3_x(2:end) , c4_x(2:end)];
    point_rec_y = [c1_y  , c2_y(2:end) , c3_y(2:end) , c4_y(2:end)];

    % Data structure for contour
    for i = 1 : sum(nb_elem_rec)
        index = nb_elem + i; 
        elem(index).p1 = [point_rec_x(i)+translation(1)   , point_rec_y(i)+translation(2) ];  
        elem(index).p2 = [point_rec_x(i+1)+translation(1) , point_rec_y(i+1)+translation(2)];  
        elem(index).l = elem(index).p2 - elem(index).p1; %vector going from p1 à p2
        ln = elem(index).l/norm(elem(index).l); %normalized vector
        elem(index).ln = ln;
        n = cross([ln 0],[0 0 1]); %normal vector of the element
        elem(index).n = n(1:2); 

        %Type of condition on the element
        if any( i == (1 : nb_elem_rec(1)))
                elem(index).type = type_rec(1); % 0=Neumann, 1=Dirichlet
                elem(index).val = val_rec(1);

        elseif any( i == (nb_elem_rec(1)+1 : sum(nb_elem_rec(1:2))))
                elem(index).type = type_rec(2); % 0=Neumann, 1=Dirichlet
                elem(index).val = val_rec(2);

        elseif any( i == (nb_elem_rec(2)+1 : sum(nb_elem_rec(1:3))))    
                elem(index).type = type_rec(3); % 0=Neumann, 1=Dirichlet
                elem(index).val = val_rec(3);

        elseif any( i == (nb_elem_rec(3)+1 : sum(nb_elem_rec(1:4))))      
                elem(index).type = type_rec(4); % 0=Neumann, 1=Dirichlet
                elem(index).val = val_rec(4);
        end

    end
    nb_cont = nb_cont + 1; %number of contour
    cont(nb_cont).nb_elem = sum(nb_elem_rec); %number of element in the contour
    cont(nb_cont).elem = [nb_elem+1 : nb_elem + sum(nb_elem_rec)];  %array of elements belonging to the contour
    cont(nb_cont).rotation = ordering; %0=clockwise, 1=anticlockwise
    nb_elem = nb_elem + sum(nb_elem_rec);
end 