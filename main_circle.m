%% Cleaning
clc;
clear all;

%% Initialization 
cont = struct(); % structure containing the information of the contour 
elem = struct(); % structure containing the information of the finite elements
nb_cont = 0; 
nb_elem = 0; 
cpt_fig = 0;

%% Parametrization and generation of the frontier of a circle

% Parameters for generating the circular frontier
ordering = 1;             % 0=clockwise, 1=anticlockwise
r = 2;                    % radius
origin = [0,0];
nb_elem_circle = 10;          % nb of elements                     
type = 1;                 % 0=Neumann, 1=Dirichlet                         
val = 1;                  % value associated with the Neumann/Dirichlet condition        

% Generation of the circular frontier
[cont, elem, nb_cont, nb_elem] = generate_circle(cont, elem, nb_cont, nb_elem, ordering, r, origin, nb_elem_circle, type, val);

%% Data structure for degrees of liberty
[ddl, nb_ddl, cont] = generate_ddl(cont, elem, nb_cont);

%% Display of degrees of liberty
cpt_fig = display_ddl(cont, elem, ddl, nb_cont, cpt_fig);

%% Matrix H and G
% The matrices H and G are construct line-by-line. For a given source
% point, ddl(i).coord (line), we compute all his interaction with every
% other elements j (columns) from de data structure "elem"

% Note : the for-loop on the elements are in the function "calc_GH0."

H = zeros(nb_ddl);
G = zeros(nb_ddl);
for i = 1:nb_ddl
    [G(i,:), H(i,:)] = calc_GH0(ddl(i).coord, elem);
    H(i,i) = H(i,i)-ddl(i).angle/(2*pi);
end


%% Solver
[u, u_prime] = solver(nb_ddl, elem, H, G);

%% Display Dirichlet and Neuman quantitites
cpt_fig = display_neumann_dirichlet(cpt_fig, ddl, nb_ddl, nb_cont, cont, u , u_prime);
    
%% Compute the artificial angle with points between elements
% We want to compute the artifical angle that are created at the
% intersection between the "j" and "j+1" element (note: 1 and n+1 are
% paired)

x_corner = zeros(1,nb_ddl);
y_corner = zeros(1, nb_ddl);
angles_corner = zeros(1, nb_ddl);

cpt = 0; 
for i = 1:nb_cont
    for j = cont(i).elem
        cpt = cpt+1;
        x_corner(cpt) = elem(j).p1(1); 
        y_corner(cpt) = elem(j).p1(2);
        
        if j == 1 % the first element
            prev_elem = cont(i).elem(end);    
        elseif i >1 && (j == 1 + cont(i-1).n_elem)
            prev_elem = cont(i).elem(end);     
        else
            prev_elem = j-1;
        end
        l1 = elem(prev_elem).ln;
        l2 = elem(j).ln;
        zz = cross([l1 0],[l2 0]);
        alpha = asin(zz(3));
        
        if cont(i).rotation == 0 %anti-horaire
             angles_corner(cpt) = pi-alpha;
        
        elseif cont(i).rotation == 1 %horaire
            angles_corner(cpt) = pi+alpha;
        end
            
    end
end 

%% Reconstruction of the solution in the domain
nb_points = 100; 
[X,Y] = meshgrid(linspace(-r , r, nb_points) , linspace(-r , r, nb_points) ); 

sol = zeros( length(X));
artifical_corner = false; %bool that helps us verify if the point is an artifical corner
for i = 1 : size(sol,1)
    for j = 1 : size(sol,2)
        artifical_corner = false;
        
        % Verify if the point is an artifical corner
        for k = 1 : length(x_corner)
            if (round(x_corner(k),8) == round(X(i,j),8)) && (round(y_corner(k),8) == round(Y(i,j),8))  
                artifical_corner = true;
                coin_angle = angles_corner(k);
            end 
        end
        
        % artifical corner
        if artifical_corner == true
            [G_domaine, H_domaine] = calc_GH0( [X(i,j), Y(i,j)], elem );
            sol(i,j) = (coin_angle/(2*pi))^(-1) * (dot(H_domaine,u) - dot(G_domaine,u_prime));
        
        % outside the domain
        elseif  sqrt( (X(i,j)).^2 + (Y(i,j)).^2 ) > (r) 
            sol(i,j) = 0;
            
        % frontier of the domain
        elseif  sqrt( (X(i,j)).^2 + (Y(i,j)).^2 ) == (r) 
            [G_domaine, H_domaine] = calc_GH0( [X(i,j), Y(i,j)], elem );
            sol(i,j) = 2* (dot(H_domaine,u) - dot(G_domaine,u_prime)); 
            
        % interior of the domain 
        else                                          
            [G_domaine, H_domaine] = calc_GH0( [X(i,j), Y(i,j)], elem );
            sol(i,j) = dot(H_domaine,u) - dot(G_domaine,u_prime) ; 
        end
         
    end
end
    
%% Solution in domain
cpt_fig = cpt_fig + 1;
figure(cpt_fig)
h = surf(X,Y,sol);
set(h,'edgecolor','none')
title('Domaine');
xlim([-r, r])
ylim([-r r])
colorbar;
grid off;