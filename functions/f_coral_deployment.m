% -------------------------------------------------------------------------
% Y.-M. Bozec, MSEL, created Aug 2018.
% Deployment of corals for restoration (with genetics)
% -------------------------------------------------------------------------

function   [coral, algal, genes, ID_colony_tracking, total_deployed] = ...
    f_coral_deployment(coral, algal, genes, REEF, ID_colony_tracking, META, Density_to_deploy, coral_diameter_mean, coral_diameter_sd, restored_cells)

% Feb 11 2022: we now have pre-determined cells for deployment (varies for each run)
% Tried to minimise code changes to keep track of previous version
select_species = find(Density_to_deploy>0);
total_deployed = zeros(1, META.nb_coral_types);

for s = select_species
    
    avail_space_cm2 = algal(1).cover_cm2 + algal(4).cover_cm2 ; % can be only outplanted on EAM and TURF

    [l,c] = size(coral(s).cover_cm2);   
    id1 = zeros(l,c);
    id1(coral(s).cover_cm2>0)=1;
    colony_count = sum(id1,2) ; % number of colonies of a given species for each cell
    
    transplant_layer = zeros(l, META.max_colonies-c);
    nb_outplants = zeros(l,1);      
            
    if META.restore_random_density == 1
        
        nb_outplants(restored_cells,1) = poissrnd(Density_to_deploy(s),1,length(restored_cells)) ;
%         nb_outplants(nb_outplants > 5) = 5;
        
        check = META.max_colonies - nb_outplants - colony_count ;
        nb_outplants(check<0) = META.max_colonies ;
    else
        
        nb_outplants(restored_cells,1) = ceil(Density_to_deploy(s))*ones(length(restored_cells),1) ;
        % should first determine total number of transplants over the grid then distribute     
    end
    % A revoir: seulement select restored cells, ajouter une ligne de
    % random diameters correspondant aux nombre d'outplants
    % puis ajuster les size avec available space
    % update available space for the next species
    
    if sum(nb_outplants)>0
        
        for k=1:length(restored_cells)
            n = nb_outplants(restored_cells(k),1);
            transplant_layer(restored_cells(k), 1:n) = normrnd(coral_diameter_mean(s), coral_diameter_sd(s),n,1);
        end
    end
    %     for j = 1:size(transplant_layer,2)
    %
    %         whereto_outplant = find(nb_outplants>0);
    %
    %         if sum(nb_outplants)>0
    %             % Need to implement cap of coral cover regarding available space
    %             do_transplant = zeros(l,1);
    %             do_transplant(nb_outplants>0)=1;
    %             transplant_layer(whereto_outplant,j)=normrnd(coral_diameter_mean(s), coral_diameter_sd(s),length(whereto_outplant),1);
    %             nb_outplants = nb_outplants - do_transplant ;
    %         end
    %     end
    
    cols = sum(transplant_layer,1);
    transplant_layer = transplant_layer(:,cols~=0) ;
    %     I = spones(transplant_layer);
    transplant_layer_cm2 = zeros(size(transplant_layer));
    transplant_layer_cm2(restored_cells,:) = ceil(pi*(transplant_layer(restored_cells,:)/2).^2); % Note the rounding can give 0 cm2 (if generated diameter <0.8)
    %     J = spones(transplant_layer_cm2);
    %     transplant_layer_cm2((I-J)==1)=1; % replace the 0 cm2 by a minimum 1

    check_space = avail_space_cm2(restored_cells) - sum(transplant_layer_cm2(restored_cells,:),2);
    K = find(check_space<0);
    
    if isempty(K)==0 % if not enough space (EAM) is available in some cells
        
        transplant_layer_cm2(restored_cells(K),:) = 0*transplant_layer_cm2(restored_cells(K),:); % just cancel deployment in those cells
    end
    
    % Check if this doesn't overpass the max number of colonies per species/cell
    check_max_allowed = META.max_colonies - size(transplant_layer_cm2,2) - c;
    
    % If it does, then delete layers of transplants
    if check_max_allowed < 0 % note this works with any negative
        transplant_layer_cm2 = transplant_layer_cm2(:,1:(size(transplant_layer_cm2, 2)+check_max_allowed));
    end
    
    total_deployed(1,s) = sum(sum(spones(transplant_layer_cm2)));
    
    outplant_ID = zeros(size(transplant_layer_cm2));
    outplant_ID(transplant_layer_cm2>0)=[1:1:total_deployed(1,s)]+ID_colony_tracking(s);

    coral(s).cover_cm2 = [coral(s).cover_cm2 transplant_layer_cm2];
    coral(s).colony_ID = [coral(s).colony_ID outplant_ID]; % marked with -1
    
%     if META.doing_clades == 1        
%         coral(s).clade = [coral(s).clade transplant_layer/outplant_size_cm2]; %transplant clade 1        
%     end
    
    ID_colony_tracking(1,s) = ID_colony_tracking(1,s)+total_deployed(1,s);% Record the new ID max

%     algal(4).cover_cm2 = algal(4).cover_cm2 - sum(transplant_layer,2); % should overgrow EAM as well
    algal(1).cover_cm2 = algal(1).cover_cm2 - sum(transplant_layer_cm2,2); % should overgrow EAM as well
    
    % Assign TQLs
    if META.doing_genetics == 1
         
%         enhance_tolerance = META.genetics.enhanced_tolerance(s)/(size(META.genetics.QTL_pool,2)*size(META.genetics.QTL_pool,3));
        enhance_tolerance = META.genetics.enhanced_tolerance(s)/(size(REEF.coral(s).QTL_pool_IN,2)*size(REEF.coral(s).QTL_pool_IN,3));
        
%         idx = randi(size(META.genetics.QTL_pool,1),total_transplanted(1,s),1);
%         TEMP_QTL = META.genetics.QTL_pool(idx,:,:) + enhance_tolerance ;
%         
%         genes(s).QTLs = [genes(s).QTLs ; TEMP_QTL] ;
%         genes(s).list_coral_ID = [genes(s).list_coral_ID ; transplant_ID(transplant_ID~=0)];
%         
%         % Compute genotypes with heritability
%         % First calculate breeding value
%         BREED = sum(sum(TEMP_QTL,3),2) ;
%         % Then deeermine phenotype following heritability (esd=0 implies perfect heritability)
%         tmp_phenotypes = BREED + normrnd(0, META.genetics.esd(s), size(BREED)) ;
%         genes(s).phenotypes = [genes(s).phenotypes ; tmp_phenotypes] ;
%         
        
        %%%%%% NEW WAY
        idx = randi(size(REEF.coral(s).QTL_pool_IN,1),total_deployed(1,s),1);
        TEMP_QTL = REEF.coral(s).QTL_pool_IN(idx,:,:) + enhance_tolerance ;
        
        genes(s).QTLs = [genes(s).QTLs ; TEMP_QTL] ;
        genes(s).list_coral_ID = [genes(s).list_coral_ID ; outplant_ID(outplant_ID~=0)];
        
        % Compute genotypes with heritability
        % First calculate breeding value
        BREED = sum(sum(TEMP_QTL,3),2) ;
        % Then determine phenotype following heritability (esd=0 implies perfect heritability)
        Env_effect = normrnd(0, META.genetics.esd(s), size(BREED));
%         Env_effect(Env_effect>4) = 4; % Limit the environmental effect (add no more that +6 deg C to BREED)
              
        % Then deeermine phenotype following heritability (esd=0 implies perfect heritability)
        tmp_phenotypes = BREED + Env_effect ;
        
        max_phenotype = floor(max(META.List_Topt)-REEF.Topt_baseline);
        tmp_phenotypes(tmp_phenotypes > max_phenotype) = max_phenotype ;        
       
        genes(s).phenotypes = [genes(s).phenotypes ; tmp_phenotypes] ;
        
    end
    
end

algal(1).cover_cm2(algal(1).cover_cm2<0)=0; % final adjustment of EAM after deployment. In rare ases EAM might be negative so forced to 0.