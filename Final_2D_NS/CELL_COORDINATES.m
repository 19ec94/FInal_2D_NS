function [cell_coord, cell_centre_coord] = ...
    CELL_COORDINATES(x_domain,y_domain,nb_cells_in_x,...
    nb_cells_in_y,total_nb_cells)

    cell_coord = zeros(5,total_nb_cells);
    cell_centre_coord = zeros(5,total_nb_cells);
for i=1:nb_cells_in_x %pointer is at individual column and progress through row
    for j=i:nb_cells_in_x:total_nb_cells
        cell_coord(1,j)=x_domain(1,i);
        cell_coord(2,j)=x_domain(1,i+1);
        cell_centre_coord(1,j)=x_domain(1,i)+ (x_domain(1,i+1)-x_domain(1,i))/2;
        cell_coord(5,j)=j;
    end
end

for i=1:nb_cells_in_y %pointer is at row and progress through row
    for j= (((i-1)*nb_cells_in_x)+1):(i*nb_cells_in_x)
        cell_coord(3,j)=y_domain(1,i);
        cell_coord(4,j)=y_domain(1,i+1);
        cell_centre_coord(2,j)=y_domain(1,i)+ (y_domain(1,i+1)-y_domain(1,i))/2;
    end
end
%for i=1:total_nb_cells
%    cell_coord(5,i)=i;
%end

end
