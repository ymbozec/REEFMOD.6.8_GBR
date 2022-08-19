% -------------------------------------------------------------------------
% Y.-M. Bozec, MSEL, created Aug 2015.
% For optimization
% -------------------------------------------------------------------------

function [coral, genes] = f_struct_arrange(coral, genes, META)

% Re-arrange the matrix coral_cm2 for optimization
% We try to limit the number of columns by moving colonies back to first colums 
% if space is available, then delete the last column if all colonies have been displaced
% Note a displaced colony stay in its cell (same row)

% Probably not worth doing this at every time step

for s = 1:META.nb_coral_types
    
    [rows, cols, ~] = size(coral(s).cover_cm2);
    a = [1:rows]' ;
    R = a(:,ones(cols,1)) ;
    [temp_coral_cm2, I] = sort(coral(s).cover_cm2,2,'descend') ;
    nIdx = R + (I-1)*rows ;
    % then re-organize the matrix of ID and clades accordingly
    temp_colony_ID = coral(s).colony_ID(nIdx);
    % check if the last column(s) has zeros
    id_col = spones(temp_coral_cm2);
    id_sum=sum(id_col,1);
    % then delete last column(s)
    temp_coral_cm2(:,id_sum==0)=[] ;
    temp_colony_ID(:,id_sum==0)=[] ;
    
    
    % Store the new (optimized) coral cover
    coral(s).cover_cm2 = sparse(temp_coral_cm2) ;
    coral(s).colony_ID = sparse(temp_colony_ID) ;
    
    if META.doing_clades == 1
        temp_clade = coral(s).clade(nIdx);
        temp_clade(:,id_sum==0)=[] ;
        coral(s).clade = sparse(temp_clade) ;
    end
    
    if META.doing_3D == 1 % then re-arrange the other matrices accordingly
        
        temp_surface_cm2 = coral(s).surface_cm2(nIdx);
        temp_volume_cm3 = coral(s).volume_cm3(nIdx);
        
        temp_surface_cm2(:,id_sum==0)=[] ;
        temp_volume_cm3(:,id_sum==0)=[] ;
        
        coral(s).surface_cm2 = sparse(temp_surface_cm2) ;
        coral(s).volume_cm3 = sparse(temp_volume_cm3) ;
        
    end
    
    % then re-arrange the matrices of TQL accordingly
    if META.doing_genetics == 1 && META.genetics.group(s)==1
        
        % First update the list of QTLs to only keep the survivors from previous mortality 
        list_old = genes(s).list_coral_ID ;
        list_new = coral(s).colony_ID(coral(s).colony_ID~=0) ;
        
        check = ismember(list_old,list_new) ;
        
        genes(s).QTLs(check==0,:,:)=[];
        genes(s).list_coral_ID(check==0,:)=[];
        genes(s).phenotypes(check==0,:)=[];
        
        % Now re-order following list_new (re-arranged matrices)
        [~,y]=sort(genes(s).list_coral_ID,'ascend');
        [~,z]=sort(list_new,'ascend');

        genes(s).QTLs(z,:,:) = genes(s).QTLs(y,:,:);
        genes(s).list_coral_ID(z,:) = genes(s).list_coral_ID(y,:);
        genes(s).phenotypes(z,:) = genes(s).phenotypes(y,:);
        
    end
end