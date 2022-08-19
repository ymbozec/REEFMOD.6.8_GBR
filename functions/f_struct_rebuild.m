% -------------------------------------------------------------------------
% Y.-M. Bozec, MSEL, created Aug 2015.
% For optimization
% -------------------------------------------------------------------------

function [coral] = f_struct_rebuild (coral_cm2, surface_cm2, volume_cm3, colony_ID, clade, species_ID, nb_coral_types, doing_clades, doing_3D)

coral(nb_coral_types).cover_cm2 =[]; % initialisation

if sum(species_ID)~=size(coral_cm2,2)

    size(species_ID)
    size(coral_cm2)
    error('species_ID not consistent with storage in coral_cm2')

end
  
col_start = 1;
col_stop = 0;

% First thing is to update with recent mortality
id0=find(coral_cm2==0 & colony_ID~=0);
colony_ID(id0)=0;


% then proceed to the reconstruction
for s = 1:nb_coral_types
      
    col_stop = col_stop + species_ID(s) ;
    
    coral(s).cover_cm2 = sparse(coral_cm2(:,col_start:col_stop));
    coral(s).colony_ID = sparse(colony_ID(:,col_start:col_stop));
    coral(s).surface_cm2 = 0 ; % need to populate surface_cm2 even if not doing 3D
    coral(s).volume_cm3 = 0 ; % need to populate volume_cm3 even if not doing 3D
    coral(s).clade = 0 ; % need to populate clade even if not doing it
    col_start = col_start + species_ID(s) ;
    
end

if doing_clades == 1
    
    clade(id0)=0;
    col_start = 1;
    col_stop = 0;
    
    for s = 1:nb_coral_types
        col_stop = col_stop + species_ID(s) ;
        coral(s).clade = sparse(clade(:,col_start:col_stop)) ;
        col_start = col_start + species_ID(s) ;
    end
    
end


if doing_3D == 1 % then re-arrange the other matrices accordingly
    
    surface_cm2(id0)=0;
    volume_cm3(id0)=0;
    
    for s = 1:nb_coral_types
        
        col_stop = col_stop + species_ID(s) ;
        
        coral(s).surface_cm2 = sparse(surface_cm2(:,col_start:col_stop));
        coral(s).volume_cm3  = sparse(volume_cm3(:,col_start:col_stop));
        
        col_start = col_start + species_ID(s) ;
        
    end
end