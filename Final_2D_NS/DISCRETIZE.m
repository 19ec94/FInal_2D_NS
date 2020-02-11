function [total_nb_cells,x_domain,y_domain,dx,dy]=DISCRETIZE (xStart, ...
    xEnd,yStart,yEnd,nb_cells_in_x,nb_cells_in_y)
%total number of cells in our domain
total_nb_cells = nb_cells_in_x* nb_cells_in_y;
%discretize the domain into nb_cells+1 points
x_domain= linspace(xStart,xEnd,nb_cells_in_x+1);
y_domain= linspace(yStart,yEnd,nb_cells_in_y+1);
dx = x_domain(2)-x_domain(1);
dy = y_domain(2)-y_domain(1);

end