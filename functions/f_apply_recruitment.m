% -------------------------------------------------------------------------
% Y.-M. Bozec, MSEL, created Nov 2011.
% Last modified: 08/2021
%
% Coral recruitment per grid cell: 'unlimited' nb of recruits, determined
% by Poisson distribution, reduced by available space (EAM only)
% -------------------------------------------------------------------------

function [coral, genes, total_settled, ID_colony_tracking, algal] = ...
    f_apply_recruitment(coral, algal, genes, META, REEF, max_density_settlers, clade_prop,ID_colony_tracking)

total_settled = zeros(1,META.nb_coral_types) ;

% random_ct = randperm(META.nb_coral_types) ; % list coral types in random order
random_ct=[randperm(6) 7:META.nb_coral_types];

EAM_cm2 = algal(1).cover_cm2 ;

for s = random_ct
    
    if max_density_settlers(s)>0
        
        free_space_cm2 = EAM_cm2;
        
        [l,c] = size(coral(s).cover_cm2);
        id1 = zeros(l,c);
        id1(coral(s).cover_cm2>0)=1;
        
        colony_count = sum(id1,2) ; % number of colonies of a given species for each cell
        
        lambda = full(max_density_settlers(s) * free_space_cm2./ REEF.floor_SA_cm2); % mean parameter of the poisson distribution of number of recruits
        
        % Adjust to available space (prop of EAM)
        new_settler_count = poissrnd(lambda, l, 1).*REEF.grazable_cell ;  %cannot settle on sand
        new_settler_count(colony_count + new_settler_count > META.max_colonies) = META.max_colonies - colony_count(colony_count + new_settler_count > META.max_colonies);
        MAX = max(new_settler_count); % max nb of candidates for settlement in a cell over the entire grid
        check_max_allowed = META.max_colonies - MAX - c;
        
        if check_max_allowed<0
            MAX = MAX + check_max_allowed; %force max number of recruits to match maximum number of corals allowed in cells
            new_settler_count(new_settler_count>MAX)=MAX;
        end
        
        if MAX==0
            continue % No recruitment for species s, so jump to the next
        end
        
        add_settler = zeros(l,MAX) ;
        
        for j = 1:MAX
            add_settler(:,j)=1; % add a settler in all cell
            add_settler(sum(add_settler,2)>new_settler_count,j)=0; % delete settlers in excess
        end
        
        
        coral(s).cover_cm2(:,c+(1:MAX))= add_settler ;
        
        total_settled(1,s) = sum(new_settler_count); % record the total number of settlers for this species
        
        settler_ID = zeros(size(add_settler));
        settler_ID(add_settler==1)=[1:1:total_settled(1,s)]+ID_colony_tracking(s);
        
        list_new = coral(s).colony_ID(coral(s).colony_ID~=0); % do this before the list of ID is extended with recruits
        
        coral(s).colony_ID(:,c+(1:MAX)) = settler_ID;% assigns ID to the new colonies, starting from ID max
        ID_colony_tracking(1,s) = ID_colony_tracking(1,s)+total_settled(1,s);% Record the new ID max
        
        % Assign TQLs
        if META.doing_genetics == 1 && META.genetics.group(s) == 1
            
            % First update the list of QTLs to only keep the survivors from previous mortality
            % (needs to be repeated here after fecundation in f_runmodel since not all the reefs were visited
            list_old = genes(s).list_coral_ID;
            check = ismember(list_old,list_new) ;
            
            genes(s).QTLs(check==0,:,:)=[];
            genes(s).list_coral_ID(check==0,:)=[];
            genes(s).phenotypes(check==0,:)=[];
            
            idx = randi(size(REEF.coral(s).QTL_pool_IN,1),total_settled(1,s),1);
            TEMP_QTL = REEF.coral(s).QTL_pool_IN(idx,:,:) ;
            
            %         idx = randi(size(REEF.coral(s).QTL_pool,1),total_settled(1,s),1);
            %         TEMP_QTL = REEF.coral(s).QTL_pool(idx,:,:) ;
            
            genes(s).QTLs = [genes(s).QTLs ; TEMP_QTL] ;
            genes(s).list_coral_ID = [genes(s).list_coral_ID ; settler_ID(settler_ID~=0)];
            
            % Compute genotypes with heritability
            % First calculate breeding value
            BREED = sum(sum(TEMP_QTL,3),2) ;
            % Then determine phenotype following heritability (esd=0 implies perfect heritability)
            Env_effect = normrnd(0, META.genetics.esd(s), size(BREED));
            %         Env_effect(Env_effect>4) = 4; % Limit the environmental effect (add no more that +6 deg C to BREED)
            
            tmp_phenotypes = BREED + Env_effect ;
            
            max_phenotype = floor(max(META.List_Topt)-REEF.Topt_baseline);
            tmp_phenotypes(tmp_phenotypes > max_phenotype) = max_phenotype ;
            
            genes(s).phenotypes = [genes(s).phenotypes ; tmp_phenotypes] ;
            
        end
        
        % Assign clade
        if META.doing_clades == 1
            rand_clade = rand(size(add_settler)).*add_settler ; % by default, all recruits have the sensitive clade (clade = 1)
            add_settler(rand_clade > clade_prop) = 2 ; % recruits have the same proportion of clade C over D than at initial step
            coral(s).clade(:,c+(1:MAX)) = add_settler ;
        end
        
        EAM_cm2 = EAM_cm2 - new_settler_count; % Update turf (only one recruit per species per cell)
        
        if META.doing_3D == 1
            coral(s).surface_cm2(:,c+(1:MAX))= add_settler;
            coral(s).volume_cm3(:,c+(1:MAX))= add_settler;
        end
        
    end
end

algal(1).cover_cm2 = sparse(EAM_cm2); %needs to be sparsed
