function par_old = INIT_POS_PAR(total_nb_cells, nb_of_particles_in_a_cell, ...
    cell_centre_coord,total_nb_particles,dx,dy)
rng('default')
par_old = zeros(1,total_nb_particles);
for i=1:total_nb_cells
    for j=(((i-1)*nb_of_particles_in_a_cell)+1):(i*nb_of_particles_in_a_cell)
        par_old(1,j)=cell_centre_coord(1,i) + dx * rand(1,1) ; %Uniform 
    end
end
for i=1:total_nb_cells
    for j=(((i-1)*nb_of_particles_in_a_cell)+1):(i*nb_of_particles_in_a_cell)
        par_old(2,j)=cell_centre_coord(2,i) + dy * rand(1,1); % uniform random
    end
end
end