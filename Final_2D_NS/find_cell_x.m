function [xl,xr] = find_cell_x(position_x,xl,xr,xm,x_domain)
if(position_x <= x_domain(1,xm))
    xl=xl; xr=xm; xm=floor((xl+xr)/2);
    if ((xr-xl) ==1)
        xl=xl; xr=xr;
        return;
    else
         [xl,xr] = find_cell_x(position_x,xl,xr,xm,x_domain);
    end
else
    xl=xm; xr=xr; xm=floor((xl+xr)/2);
    if((xr-xl) ==1)
        xl=xl;xr=xr;
        return;
    else
        [xl,xr]= find_cell_x(position_x,xl,xr,xm,x_domain);
    end
end
end
