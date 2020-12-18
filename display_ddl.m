function cpt_fig = display_ddl(cont, elem, ddl, nb_cont, cpt_fig)

    cpt_fig = cpt_fig +1;
    figure(cpt_fig)
    hold on
    for i = 1 : nb_cont
        for j = 1:cont(i).nb_elem
            xy = [elem(cont(i).elem(j)).p1 ; elem(cont(i).elem(j)).p2];
            line(xy(:,1),xy(:,2));
        end
        for j = 1:cont(i).nb_ddl
            xy = ddl(cont(i).ddl(j)).coord;
            line(xy(1),xy(2),'Marker','o');
        end
    end
    axis equal
end