function z = f_mean_traj(x,area_w)

y = (area_w')*x(:,1:end)/sum(area_w);

z= nan(size(y));
wd = 1; 
for i=(1+wd):(size(y,2)-wd)
    
    z(1,i)=mean(y(i-wd:i+wd));
    
end

