function [yd,yu] = find_cell_y(position_y,yd,yu,ym,y_domain)
if(position_y <= y_domain(1,ym))
    yd=yd; yu=ym; ym=floor((yd+yu)/2);
    if ((yu-yd) ==1)
        yd=yd; yu=yu;
        return;
    else
         [yd,yu] = find_cell_x(position_y,yd,yu,ym,y_domain);
    end
else
    yd=ym; yu=yu; ym=floor((yd+yu)/2);
    if((yu-yd) ==1)
        yd=yd;yu=yu;
        return;
    else
        [yd,yu]= find_cell_x(position_y,yd,yu,ym,y_domain);
    end
end
end