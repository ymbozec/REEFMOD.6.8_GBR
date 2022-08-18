%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Y.-M. Bozec, MSEL, created Aug 2015.
%
% Track individuals for John Hedley's visualization tool
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [environ_list,colony_list] = f_track_populations(environ_list,colony_list,algal,coral,t,META,grazable_cell)


% FIRST RECORD SAND
[X_sand, Y_sand] = ind2sub([META.grid_x_count,META.grid_y_count],find(grazable_cell==0));
col_zero = zeros(size(X_sand,1),1);
sand_temp = [col_zero+t,X_sand, Y_sand,col_zero+META.cell_area_cm2,col_zero];
environ_list = [environ_list ; sand_temp] ;

% THEN RECORD ALGAE
algal_cm2 = full(([algal(:).cover_cm2])) ;
[cells,J]=size(algal_cm2(grazable_cell~=0,:)) ;

for a = 1:META.nb_algal_types
    
    list_temp = zeros(cells,5);
    list_temp(:,1) = t;
    [list_temp(:,2), list_temp(:,3)] = ind2sub([META.grid_x_count,META.grid_y_count],find(grazable_cell~=0));
    list_temp(:,4) = algal_cm2(grazable_cell~=0,a) ;
    list_temp(:,5) = a ;
    
    environ_list = [environ_list ; list_temp] ;
    
end

% NOW RECORD CORAL COLONIES
for s = 1:META.nb_coral_types
    
    id1=zeros(size(coral(s).cover_cm2));
    id1(coral(s).cover_cm2>0) = 1;
    
    nb_colonies = nnz(coral(s).cover_cm2(id1==1)) ;
    
    list_temp = zeros(nb_colonies, 6);
    
    list_temp(:,1) = t ;
    [cell,J]=ind2sub(size(coral(s).cover_cm2),find(id1==1));
    
    [list_temp(:,2),list_temp(:,3)] = ind2sub([META.grid_x_count,META.grid_y_count],cell);
    list_temp(:,4)=coral(s).cover_cm2(id1==1);
    list_temp(:,5)=coral(s).colony_ID(id1==1);
    
    list_temp(:,6)=s;
    
    colony_list = [colony_list ; list_temp] ;
end