function [pp,current_cell,new_cell_curl]= NEW_CELL_CURL (current_cell,total_nb_particles,total_nb_cells,nb_cells_in_x, ...
    nb_cells_in_y,x_domain,y_domain,par_new,new_cell_curl,pp)
num_par_in_cell= zeros(1,total_nb_cells);
new_cell_curl = 0.*new_cell_curl;
position_x =0; position_y=0;
for par_number=1:total_nb_particles
    %There is a scope to paralalize this loop, as it runs thorugh number of
    %particles.For each particle  it calculates the cell number and
    %accumulates in the new_cell_curl and accumulates count in num_par_in_a_cell.
    %Attention should be paid to the "new_cell_curl","num_par_in_a_cell".
    %Crtical or atomic would be fine.
    xl=1; xr=nb_cells_in_x+1; xm=floor((xl+xr)/2);
    position_x = par_new(1,par_number);
    [xl,xr] = find_cell_x(position_x,xl,xr,xm,x_domain);
    position_y = par_new(2,par_number);
    yd=1; yu=nb_cells_in_y+1; ym=floor((yd+yu)/2);
    [yd,yu] = find_cell_y(position_y,yd,yu,ym,y_domain);
    current_cell(1,par_number)=(xl+(yd-1)*nb_cells_in_x);
    %Note xl and yd are the actual new cell identification
    %I have to indentify the actual new_curl_cell array to store
    %the curl value
    new_cell_curl( (xl+(yd-1)*nb_cells_in_x) ) =  ...
        new_cell_curl((xl+(yd-1)*nb_cells_in_x)) + par_new(5,par_number);
    %I also store the count of the number of particles in the
    %"xl+(yd-1)*nb_cells_in_x" cell. 
    num_par_in_cell((xl+(yd-1)*nb_cells_in_x))= ...
        num_par_in_cell((xl+(yd-1)*nb_cells_in_x))+1;
end
%New average cell curl
%I divide each cell by the current number of particles I get the average
%value of curl of the cell
new_cell_curl = new_cell_curl ./ num_par_in_cell;
pp = num_par_in_cell;
end