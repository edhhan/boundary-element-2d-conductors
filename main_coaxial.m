%% Cleaning
clc;
clear all;

%% Initialization 
cont = struct(); %structure containing the information of the contour 
elem = struct(); %structure containing the information of the finite elements
nb_cont = 0; 
nb_elem = 0; 
cpt_fig = 0;

nb_elem_circle = 200; %number of per circle

%% Parametrization and generation of the outer circle

% Parameters for generating the circular frontier
ordering1 = 1;
r1 = 12;                    % radius
origin1 = [0,0];            % origin                                          
nb_elem_circle1 = nb_elem_circle;   % nb of elements for circle 1                     
type1 = 1;                  % 0=Neumann, 1=Dirichlet                         
val1 = 0;                   % value associated with the Neumann/Dirichlet condition          

% Generation of the circular frontier
[cont, elem, nb_cont, nb_elem] = generate_circle(cont, elem, nb_cont, nb_elem, ordering1, r1, origin1, nb_elem_circle1, type1, val1);


%% Parametrization and generation of the inner circle

% Parameters for generating the circular frontier
ordering2 = 0;
r2 = 4;                     % radius
origin2 = [0,0];            % origin
nel_circle2 = nb_elem_circle;   % nb of elements for circle 2                  
type2 = 1;                  % 0=Neumann, 1=Dirichlet                         
val2 = 1;                   % value associated with the Neumann/Dirichlet condition        

% Generation of the circular frontier
[cont, elem, nb_cont, nb_elem] = generate_circle(cont, elem, nb_cont, nb_elem, ordering2, r2, origin2, nel_circle2, type2, val2);

%% Data structure for degrees of liberty
[ddl, nb_ddl, cont] = generate_ddl(cont, elem, nb_cont);

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

%% Display of degrees of liberty
cpt_fig = display_ddl(cont, elem, ddl, nb_cont, cpt_fig);

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

cpt = 0; %compteur
for i = 1:nb_cont
    for j = cont(i).elem
        cpt = cpt+1;
        x_corner(cpt) = elem(j).p1(1); 
        y_corner(cpt) = elem(j).p1(2);
        
        if j == 1 % the first element
            prev_elem = cont(i).elem(end);    
        elseif i >1 && (j == 1 + cont(i-1).nb_elem)  
            prev_elem = cont(i).elem(end);     
        else
            prev_elem = j-1;
        end
        l1 = elem(prev_elem).ln;
        l2 = elem(j).ln;
        zz = cross([l1 0],[l2 0]);
        alpha = asin(zz(3));
        
        % clockwise=0, anticlockwise=1
        if cont(i).rotation == 1
             angles_corner(cpt) = pi-alpha;
        
        elseif cont(i).rotation == 0 
            angles_corner(cpt) = pi+alpha;
        end
            
    end
end 

%% Reconstruction of the solution in the domain
nb_points = 251; 
d = 2*r1; %domain displaying
[X,Y] = meshgrid(linspace(-(d/2) , (d/2), nb_points) , linspace(-(d/2) , (d/2), nb_points) ); 

sol = zeros( length(X));
for i = 1 : size(sol,1)
    for j = 1 : size(sol,2)
        
        %r1 : outer circle
        %r2 : inner circle
        
        artifical_corner = false; 
        for k = 1 : nb_ddl
            if (round(x_corner(k),10) == round(X(i,j),10)) && (round(y_corner(k),10) == round(Y(i,j),10))  
                artifical_corner = true;
                coin_angle = angles_corner(k);
            end 
        end
        
        
        % frontier of a circle
        if sqrt( (X(i,j).^2 + Y(i,j).^2) ) == r1 || sqrt( (X(i,j).^2 + Y(i,j).^2) ) == r2
            [G_domaine, H_domaine] = calc_GH0( [X(i,j), Y(i,j)], elem );
            sol(i,j) = 2*(dot(H_domaine,u) - dot(G_domaine,u_prime));
        
        % artifical corner
        elseif artifical_corner == true
            [G_domaine, H_domaine] = calc_GH0( [X(i,j), Y(i,j)], elem );
            sol(i,j) = (coin_angle/(2*pi))^(-1) * (dot(H_domaine,u) - dot(G_domaine,u_prime));
            
        % outside domain
        elseif sqrt( (X(i,j).^2 + Y(i,j).^2) ) < r2 || sqrt( (X(i,j).^2 + Y(i,j).^2) ) > r1
            sol(i,j) = NaN;
            
        % interior of domain
        else
            [G_domaine, H_domaine] = calc_GH0( [X(i,j), Y(i,j)], elem );
            sol(i,j) = (dot(H_domaine,u) - dot(G_domaine,u_prime));
        end
        
    end
end
    
%% Solution in domain
cpt_fig = cpt_fig + 1;
figure(cpt_fig)
h = surf(X,Y,sol);
set(h,'edgecolor','none')
%title('Domaine');
colorbar;
hold on
contour3(X,Y,sol,11,'w');
hold off
grid off;



