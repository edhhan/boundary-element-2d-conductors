function [ddl, nb_ddl, cont] = generate_ddl(cont, elem, nb_cont)

    ddl = struct(); 
    nb_ddl = 0;
    for i = 1 : nb_cont
        ddl_list = [];
        for j = cont(i).elem
            nb_ddl = nb_ddl+1;
            ddl_list = [ddl_list nb_ddl];
            ddl(nb_ddl).coord = (elem(j).p1 + elem(j).p2)/2;
            ddl(nb_ddl).angle = pi; %mid-point
            elem(j).ddl = nb_ddl;       
        end
        cont(i).ddl = ddl_list;
        cont(i).nb_ddl = length(ddl_list); 
    end
    
end