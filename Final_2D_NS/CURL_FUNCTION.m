function cell_curl = CURL_FUNCTION(total_nb_cells,cell_centre_coord)
cell_curl = zeros(1,total_nb_cells);
syms vel_vec(x,y,z) %creation of symbolic variables
vel_vec(x,y,z) = [sin(2*pi*y) cos(2*pi*x) 0]; %velocity vector
sm_curl = curl(vel_vec,[x y z]); % curl of the velocity vector and it is 
%still a symbolic function. so we use "matlabFunction" to convert it to 
% a oridinary function func_sm_curl. And then we pass arguments to this
%function which gives us normal numerical values.
func_sm_curl = matlabFunction(sm_curl,'vars', [x y z],'Optimize',false);
for i=1:total_nb_cells
    %here "numerical_sm_curl" is a 3D vector with first two elements zero
    %and only the third component has a numerical value.
    numerical_sm_curl = func_sm_curl(cell_centre_coord(1,i), ...
         cell_centre_coord(2,i),0); 
     %so we store only the 3rd component in our cell curl array
    cell_curl(i) = numerical_sm_curl(3);
end
end