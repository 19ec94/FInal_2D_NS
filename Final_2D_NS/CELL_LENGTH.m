function [x_cell,y_cell] = CELL_LENGTH(nb_cells_in_x,nb_cells_in_y,x_domain, ...
    y_domain)
x_cell = zeros(nb_cells_in_x);
y_cell = zeros(nb_cells_in_y);
%We store each cell's absolute length from orgin in "x" direction
for index=1:nb_cells_in_x
    x_cell(index)=x_domain(index+1);
end
%We store each cell's absolute length from orgin in "y" direction
for index=1:nb_cells_in_y
    y_cell(index)=y_domain(index+1);
end
end