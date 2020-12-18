function [cpt_fig] = display_neumann_dirichlet(cpt_fig, ddl, nb_ddl, nb_cont, cont, u , u_prime)

    x_array = zeros(1,nb_ddl);
    y_array = zeros(1, nb_ddl);
    for i = 1 : nb_ddl
       x_array(i) = ddl(i).coord(1);
       y_array(i) = ddl(i).coord(2); 
    end

    cpt_fig = cpt_fig + 1;
    figure(cpt_fig);
    for i = 1:nb_cont
        plot3(x_array(cont(i).ddl) , y_array(cont(i).ddl), u(cont(i).ddl) );
        hold on
    end
    title('Dirichlet');
    hold off

    cpt_fig = cpt_fig + 1;
    figure(cpt_fig);
    for i = 1:nb_cont
        plot3(x_array(cont(i).ddl) , y_array(cont(i).ddl), u_prime(cont(i).ddl) );
        hold on
    end
    title('Neumann');
    hold off


end