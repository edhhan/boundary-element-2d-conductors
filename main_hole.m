%% Cleaning
clc;
clear all;

%% Initialization 
cont = struct(); %structure containing the information of the contour 
elem = struct(); %structure containing the information of the finite elements
nb_cont = 0; 
nb_elem = 0; 
cpt_fig = 0;

%% Parametrization and generation of the frontier of a outer rectangle

% Parameters for generating the rectangular frontier
ordering = 1; % 0=clockwise, 1=anticlockwise
d1 = 10; % length for pair-wise opposite side : right and left  
d2 = 10; % length for pair-wise opposite side : top and bottom
translation = [0,0];
nb_elem_d = 25; 
nb_elem_rec =  [nb_elem_d , nb_elem_d , nb_elem_d, nb_elem_d];   
type_rec = [1     ,   1   ,   1  ,     1];   % 0=Neumann, 1=Dirichlet
val_rec =  [0     ,   0   ,   0  ,     0];   % value associated with the Neumann/Dirihlet condition per side

% Generation of the rectangular frontier
[cont, elem, nb_cont, nb_elem] = generate_rectangle(cont, elem, nb_cont, nb_elem, ordering, d1, d2, translation, nb_elem_rec, type_rec, val_rec);

%% Parametrization and generation of the frontier of the inner circle

% Parameters for generating the rectangular frontier
ordering = 0; %0=clockwise, 1=anticlockwise
r = 2;                     % radius 
origin = [0,0];
nb_elem_circle = 200;       % nb of element                    
type = 1;                  % 0=Neumann, 1=Dirichlet                         
val = 1;                   % value associated with the Neumann/Dirichlet condition        

% Generation of the circle
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

%% Reconstruction of the solution in the domain
nb_points = 100;
[X,Y] = meshgrid(linspace(-d2/2 , d2/2 , nb_points), linspace(-d1/2 , d1/2 , nb_points)); 
sol = zeros(size(X,1), size(X,2));

for i = 1 : size(sol,1)
    for j = 1 : size(sol,2)

        % corner of domain
        if abs(X(i,j)) == d2/2 && abs(Y(i,j)) == d1/2
            [G_domaine, H_domaine] = calc_GH0( [X(i,j), Y(i,j)], elem );
            sol(i,j) = 4* (dot(H_domaine,u) - dot(G_domaine,u_prime));

        % frontier of domain
        elseif abs(X(i,j)) == d2/2 || abs(Y(i,j)) == d1/2
            [G_domaine, H_domaine] = calc_GH0( [X(i,j), Y(i,j)], elem );
            sol(i,j) = 2* (dot(H_domaine,u) - dot(G_domaine,u_prime));
         
        % interior of circle (undefined domain)
        elseif sqrt( (X(i,j)^2 + Y(i,j)^2) )< r
            sol(i,j) = NaN;

        % interior of domain : between circle and rectangle (defined domain)
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
contour3(X,Y,sol,10,'w');
hold off
grid off;
