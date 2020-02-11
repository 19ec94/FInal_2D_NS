%% Domain
%We want to solve NS equation in a 2D domain. For that we define our domain.
%we specify "x and y" domain.
%Start and End point in x and y direction
xStart = 0; xEnd = 1;  yStart = 0; yEnd= 1;
%We want to discretize the "x" and "y" domain into "N" and "N" number of
% sub-domians. Importantly we can have different number of cells in "x" and
%"y" directions. So we need two variables.
nb_cells_in_x = 20; nb_cells_in_y = 20;
[total_nb_cells,x_domain,y_domain,dx,dy]=DISCRETIZE (xStart, ...
    xEnd,yStart,yEnd,nb_cells_in_x,nb_cells_in_y);
%We store each cell's absolute length in "x" and "y" directions from the
%respective origin
[x_cell,y_cell] = CELL_LENGTH(nb_cells_in_x,nb_cells_in_y,x_domain, ...
    y_domain);
%Number of particle in a cell
nb_of_particles_in_a_cell = 100;
%Total_number_of_particles
total_nb_particles= nb_of_particles_in_a_cell*nb_cells_in_x*nb_cells_in_y;
%populate the cells with particles. It is  5*number_of_particles array
%that stores the whole data about the particle's  x,y,
%velocity_x,velocity_y,vorticity  of the cell it comes from
%par_old = zeros(5,total_nb_particles);
par_new = zeros(5,total_nb_particles);
%diffusion coefficient
diff_co_eff=0.01;
%Each cell has a 2D velocity vector. And particles have the velocity of the
%cell it comes from.  We define two variables to store the values at time
%"t-1" and "t"
cell_vel = zeros(2,total_nb_cells);
new_cell_vel = zeros(2,total_nb_cells);
%cell_curl = zeros(1,total_nb_cells);
new_cell_curl = zeros(1,total_nb_cells);
current_cell = zeros(1,total_nb_particles);
%time_step
num_par_in_a_cell=0;
pp = zeros(1,total_nb_cells);

%We store each cells xStart,xEnd, yStart and yEnd and it's cell number.
%We can use it later ?
%cell_coord = zeros(5,total_nb_cells);
%We store the centre point of each cell. This is used in initializing
%particles at this position and calculating velocity and curl at this point.
%cell_centre_coord = zeros(2,total_nb_cells);
[cell_coord, cell_centre_coord] = ...
    CELL_COORDINATES(x_domain,y_domain,nb_cells_in_x,...
    nb_cells_in_y,total_nb_cells);
% Initialize each particle's position in the domain
%We set each particle's initial position as centre of the respective cell
par_old = INIT_POS_PAR(total_nb_cells, nb_of_particles_in_a_cell, ...
    cell_centre_coord,total_nb_particles,dx,dy);
% initial velocity
%we calculate cell velocity at it's centre position
%cell_curl = zeros(1,total_nb_cells);
for i=1:total_nb_cells
        cell_vel(1,i)=sin(2*pi* cell_centre_coord(2,i));
        cell_vel(2,i)=cos(2 *pi*cell_centre_coord(1,i));
      % cell_curl(1,i)= -( 2*pi*cos(2 *pi*cell_centre_coord(2,i)) )-( 2*pi*sin(2*pi*cell_centre_coord(1,i)) ) ;
end
%**********TIME STEP**************
CFL=0.9;
final_time=0.0005;
Xmax=max(cell_vel(1,:));
Ymax=max(cell_vel(2,:));
if Xmax>Ymax
    dt = (CFL*dx) / Xmax;
elseif Ymax>Xmax
    dt = (CFL*dx) / Ymax;
else
    dt = (CFL*dx) /Xmax;
end
%***********TIME STEP*************
% Initial cell curls
%we calculate cell curl at it's centre position
%*************************CURL**************************
cell_curl1 = CURL_FUNCTION(total_nb_cells,cell_centre_coord);
cell_curl =cell_curl1;
X =zeros(nb_cells_in_x,nb_cells_in_x);
Y =zeros(nb_cells_in_x,nb_cells_in_x);
U =zeros(nb_cells_in_x,nb_cells_in_x);
V =zeros(nb_cells_in_x,nb_cells_in_x);
wo =zeros(nb_cells_in_x,nb_cells_in_x);
for i=1:nb_cells_in_x
    for j=1:nb_cells_in_x
        jj = (i-1)*nb_cells_in_x+j;
        X(i,j) = cell_centre_coord(1,jj);
        Y(i,j) = cell_centre_coord(2,jj);
        U(i,j) = cell_vel(1,jj);
        V(i,j) = cell_vel(2,jj);
     end
end
[wr,ca]=curl(X,Y,U,V);
for i=1:nb_cells_in_x
    for j=1:nb_cells_in_x
        jj = (i-1)*nb_cells_in_x+j;
        wo(i,j)=cell_curl1(1,jj);
    end
end
for i=1:nb_cells_in_x
    for j=1:nb_cells_in_x
        jj = (i-1)*nb_cells_in_x+j;
          cell_curl(1,jj)=wr(i,j);
    end
end
cell_curl=cell_curl1;
%figure(1)
%surf(X,Y,wo)
%figure(2)
%surf(X,Y,wr)

%**********************CURL*******************************

% Initialize each particle's velocity and vorticity
for i = 1:total_nb_cells
    for xInd=(((i-1)*nb_of_particles_in_a_cell)+1):(i*nb_of_particles_in_a_cell)
        par_old(3,xInd)=cell_vel(1,i);
        par_old(4,xInd)=cell_vel(2,i);
        par_old(5,xInd)=cell_curl(1,i);
    end
end
% copy
par_new(1,:) = par_old(1,:);
par_new(2,:) = par_old(2,:);
par_new(5,:) = par_old(5,:);
t=0;
%%
for time=0:dt:final_time
    t=t+1;
    % update the particle position
    [par_new] = POS_UPDATE(total_nb_particles,par_old,par_new,diff_co_eff,dt,...
        xStart,xEnd,yStart,yEnd);
    
    % update cell vorticity
    [pp,current_cell,new_cell_curl]= NEW_CELL_CURL (current_cell,total_nb_particles,total_nb_cells,nb_cells_in_x, ...
        nb_cells_in_y,x_domain,y_domain,par_new,new_cell_curl,pp);
    
    %cell velocity after celculating only from the respective curl
    new_cell_vel = 0 .* new_cell_vel;
    for j=1:total_nb_cells
        for i=1:total_nb_particles
            dist_x = cell_centre_coord(1,j)-par_new(1,i);
            dist_y = cell_centre_coord(2,j)-par_new(2,i);
            modulo = dist_x^2+dist_y^2;
            new_cell_vel(1,j)= new_cell_vel(1,j) + par_new(5,i)* Kx(modulo,dist_y);
            new_cell_vel(2,j)= new_cell_vel(2,j) + par_new(5,i)* Ky(modulo,dist_x);
        end
        new_cell_vel(1,j) = (1/total_nb_particles) * new_cell_vel(1,j);
        new_cell_vel(2,j) = (1/total_nb_particles) * new_cell_vel(2,j); %% check for paranthesis (2*pi)
    end
    %new_cell_vel = new_cell_vel ./ total_nb_particles;
    %********************New REconstructed CELL CURL**********************
    
    
    
    
    %********************************************************
    %Update newposition to old position and new average velocity to old par
    %and new par
    for i=1:total_nb_particles
        r = current_cell(1,i);
        par_old(3,i)= new_cell_vel(1,r);
        par_old(4,i)= new_cell_vel(2,r);
        par_new(5,i)= new_cell_curl(1,r);
    end
end
%%

x = cell_centre_coord(1,:);
y = cell_centre_coord(2,:);
%%


u1 = cell_vel(1,:);
v1=cell_vel(2,:);
%subplot(1,2,1)
quiver(x,y,u1,v1,'linewidth',2);
title('t=0');
axis([0 1 0 1],'square');
legend('u(x,y)');
xlabel('x');
ylabel('y');
u = new_cell_vel(1,:);
v=new_cell_vel(2,:);
%subplot(1,2,2)
%{
quiver(x,y,u,v,'linewidth',2);
title('t=0.0005');
axis([0 1 0 1],'square');
legend('u(x,y)');
xlabel('x');
ylabel('y');

%}
%%
%{
fileID1=fopen('U10.txt','w');
fileID2=fopen('U20.txt','w');
fileID3=fopen('U40.txt','w');
fileID4=fopen('U60.txt','w');
fprintf(fileID1,'%16f, ',n1);
fprintf(fileID2,'%16f, ',n2);
fprintf(fileID3,'%16f, ',n3);
fprintf(fileID4,'%16f, ',n4);
fclose(fileID1);
fclose(fileID2);
fclose(fileID3);
fclose(fileID4);
%}
%%
%{
e1=zeros(2,100);
e2=zeros(2,400);
e3=zeros(2,1600);
e4=zeros(2,3600);

v1=zeros(1,100);
v2=zeros(1,400);
v3=zeros(1,1600);
v4=zeros(1,3600);

for i=1:100
    e1(1,i)=n1(1,i)-n5(1,64*i);
    e1(2,i)=n1(2,i)-n5(2,64*i);
    v1(1,i)=n5(1,36*i);
end
no1=norm(e1)/norm(v1)

for i=1:400
    e2(1,i)=n2(1,i)-n5(1,16*i);
    e2(2,i)=n2(2,i)-n5(2,16*i);
    v2(1,i)=n5(1,16*i);
end
no2=norm(e2)/norm(v2)

for i=1:1600
    e3(1,i)=n3(1,i)-n5(1,4*i);
    e3(2,i)=n3(2,i)-n5(2,4*i);
    v3(1,i)=n5(1,4*i);
end
no3=norm(e3)/norm(v3)

%}
%%
%{
grid = [10,20,40];
error=[1.4044,1.3672,1.3525];
loglog(grid,error);
title('Grid convergence test');
xlabel('log(N)');
ylabel('log(\mid\mide\mid\mid_{2})');
legend({'Error'})
box off
%}
%%
%matlab2tikz('c_ns_ec.tex','width','\fwidth','extraAxisOptions','enlargelimits=false')
%%
function vector = Kx(modulo,value)
%modulo = sqrt(modulo);
vector = (1/(2*pi*modulo))*(-(value));
end
function vector = Ky(modulo,value)
%modulo = sqrt(modulo);
vector = (1/(2*pi*modulo))* ((value));
end
