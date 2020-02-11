function [par_new] = POS_UPDATE (total_nb_particles,par_old,par_new,diff_co_eff,dt,...
    xStart,xEnd,yStart,yEnd)

for i=1:total_nb_particles
    par_new(1,i)=par_old(1,i)+dt * par_old(3,i) + sqrt(2*diff_co_eff)* sqrt(dt)*randn;
    %periodic boundary condition
    if par_new(1,i) > xEnd
        par_new(1,i) = par_new(1,i)-xEnd;
    elseif par_new(1,i) < xStart
        par_new(1,i) = par_new(1,i)+xEnd;
    end
    par_new(2,i)=par_old(2,i)+dt * par_old(4,i) + sqrt(2*diff_co_eff)* sqrt(dt)*randn;
    if par_new(2,i) > yEnd
        par_new(2,i) = par_new(2,i)-yEnd;
    elseif par_new(2,i) < yStart
        par_new(2,i) = par_new(1,i)+yEnd;
    end
end

end