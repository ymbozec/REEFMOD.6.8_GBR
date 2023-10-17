% -------------------------------------------------------------------------
% Y.-M. Bozec, MSEL, created Nov 2011.
%
% Last modified: 07/2020 (track coral loss per species)
%
% Modified (08/2018) to integrate genetic adaptation to warming
% -------------------------------------------------------------------------


function [coral, genes, algal, total_coral_loss, total_mortality] = f_bleaching_new(coral, genes, algal, bleaching_whole_mortality,...
    CORAL, doing_3D, nb_coral_types, doing_clades, doing_genetics, bleaching_whole_offset, bleaching_partial_offset,Topt_baseline, Topt2index)

% NEED as inputs
% - matrix of DHW for the given reef for each MMM (Topt) at time step

% Extract data from the structures (need to be filled again when leaving)
algal_cm2 = [algal.cover_cm2] ;
[coral_cm2, surface_cm2, volume_cm3, clade, colony_ID, species_ID] = f_struct_deploy (coral);

% Locate brooders and spawners
% id_brooders = zeros(size(coral_cm2)) ; % initialize brooders id
% id_spawners = zeros(size(coral_cm2)) ; % initialize spawners id

%%%% This is new stuff (August 2013) for implementing species-specific bleaching mortalities
id0 = zeros(size(coral_cm2)) ;

sensitivity_bleaching = id0 ;
extent_bleaching = id0 ;
proba_switching = id0 ;

id1 = id0+1;
id1(coral_cm2 <= 0) = 0 ; % Remove the already dead ones

col_start = 1;
col_stop = 0;

%% Added thermal resistance from genetics (August 2018)
if doing_genetics == 1
    
    MORT = id0 ;
    
    for s = 1:nb_coral_types
        
        if species_ID(s)>0
            
            col_stop = col_stop + species_ID(s) ;
            
            id1_tmp = id1(:,col_start:col_stop);
            s_mort = zeros(size(id1_tmp)) ;
            colony_ID_tmp = colony_ID(:,col_start:col_stop);
            
            % First update the list of QTLs to only keep the survivors from previous mortality
            list_old = genes(s).list_coral_ID;
            list_new = colony_ID_tmp(id1_tmp==1);
            
            check = ismember(list_old,list_new) ;
            
            genes(s).QTLs(check==0,:,:)=[];
            genes(s).list_coral_ID(check==0,:)=[];
            genes(s).phenotypes(check==0,:)=[];
            
            % Phenotype is Topt so need to find out the corresponding DHW
            % Then convert into whole colony mortality using Hughes relationship
            % Output must be a matrix of proba of whole colony-mortality
            All_Topts = round(10*(Topt_baseline + genes(s).phenotypes))/10 ;
            
            if isempty(All_Topts)==0
                
%                 max(All_Topts)               
                s_mort(id1_tmp==1) = bleaching_whole_mortality(10*All_Topts - Topt2index) ;
                
            end
            
            sensitivity_bleaching(:,col_start:col_stop) = CORAL.sensitivity_bleaching(s);
            extent_bleaching(:,col_start:col_stop) = CORAL.bleaching_partial_extent(s);
            proba_switching(:,col_start:col_stop)= CORAL.proba_switching(s) ;
            MORT(:,col_start:col_stop) = s_mort ;
            
            col_start = col_start + species_ID(s) ;
            
        end
    end
    
else
    
    MORT = bleaching_whole_mortality.*id1 ;
    
    for s = 1:nb_coral_types
        
        col_stop = col_stop + species_ID(s) ;
        sensitivity_bleaching(:,col_start:col_stop)= CORAL.sensitivity_bleaching(s);
        extent_bleaching(:,col_start:col_stop)= CORAL.bleaching_partial_extent(s);
        proba_switching(:,col_start:col_stop)= CORAL.proba_switching(s) ;
        
        col_start = col_start + species_ID(s) ;
        
    end
    
end

if doing_clades == 1
    % Clade-induced tolerance to thermal stress
    sensitivity_bleaching(clade==2) = sensitivity_bleaching(clade==2) * CORAL.bleaching_tolerance_clade ;
end

%________________________________
%
% Whole-colony mortality
%________________________________

% Generate random mortalities for adol + adults
id1 = ones(size(coral_cm2)) ; % assigns 1 to every colony
% id1(coral_cm2 < CORAL.adol_size) = 0 ; % Only keep adults + adol (note this also excludes negative/dead colonies for speed)
id1(coral_cm2 < CORAL.size_threshold_wcm) = 0 ; % Only keep adults + adol (note this also excludes negative/dead colonies for speed)

rand_mort1 = rand(size(id1)) ; % Generates random probability for adol + adults

prob_initial_mortality = sensitivity_bleaching*CORAL.bleaching_depth.*MORT;
prob_initial_mortality(prob_initial_mortality>1)=1; % cap to 1 before extrapolating to 6 month
prob_whole_mortality = 1-(1-prob_initial_mortality).^bleaching_whole_offset;

% prob_whole_mortality = id1.*(1-(1-sensitivity_bleaching.*MORT).^bleaching_whole_offset);
% prob_whole_mortality = id1.*sensitivity_bleaching .* MORT * bleaching_whole_offset;

id_dead = id1 ;
id_dead(rand_mort1 > prob_whole_mortality) = 0 ; % exclude the survivors
coral_loss = coral_cm2.*id_dead ;

algal_cm2(:,1) = algal_cm2(:,1) + sum(coral_loss,2) ;
% coral_cm2(id_dead==1) = - coral_cm2(id_dead==1);  % now dead (negatives)
coral_cm2(id_dead==1) = 0; 
total_mortality = sum(sum(id_dead))/sum(sum(id1));

%________________________________
%
% Partial mortality
%________________________________

% These are used to update the status of each coral. If a coral has never been bleached
% and then survives a bleaching event it's status changes from 0 to 1.
% NOTE from JH: the record of previous bleaching doesn't  have any effect on likelihood
% of partial mortality whereas it does on total mortality, doesn't make sense.
% NOTE from YM: we now (01/2015) apply the same reduction to the probability of partial mortality 

id2 = id1 - id_dead ;
rand_mort2 = rand(size(id2)) ; % Generates random probability for adol + adults

prob_partial_mortality = id2.*(1-(1-prob_initial_mortality).^bleaching_partial_offset); % May 2023 -> use initial mortality which is already capped to 1
% prob_partial_mortality = id2.*(1-(1-sensitivity_bleaching.*MORT).^bleaching_partial_offset);
% prob_partial_mortality = id2.* sensitivity_bleaching .* MORT * bleaching_partial_offset;

id_part = id2 ;
id_part(rand_mort2 > prob_partial_mortality) = 0 ; 
bleach_extent = floor(extent_bleaching.* coral_cm2 .* id_part) ; 
algal_cm2(:,1) = algal_cm2(:,1) + sum(bleach_extent, 2) ;
coral_cm2 = coral_cm2 - bleach_extent ;

%________________________________
%
% Clade switching
%________________________________
if doing_clades == 1
    rand_switch = rand(size(id2)) ; % Generates random probability of switching to the thermally-tolerant clade (clade 2)
    clade(rand_switch < proba_switching & coral_cm2>0) = 2 ; % NOTE THIS INDEPENDENT OF BLEACHING MORTALITY
end

%%%%%%%% Before leaving, store the new covers into 'coral' and 'algal'%%%%%%%%%%%%%%%%%%%%
[coral] = f_struct_rebuild (coral_cm2, surface_cm2, volume_cm3, colony_ID, clade, species_ID, nb_coral_types, doing_clades, doing_3D);

for a=1:size(algal_cm2, 2) 
    algal(a).cover_cm2(:,1) = algal_cm2(:,a) ;
end

% total_coral_loss = sum(sum(coral_loss,2)) + sum(sum(bleach_extent, 2)) ;

% NEW Jul 2020: keep track of losses per species
count = 1;
total_coral_loss = zeros(1,nb_coral_types);

for s = 1:nb_coral_types
    select_sp = count:(count+species_ID(s)-1);
    total_coral_loss(s) = sum(sum(coral_loss(:,select_sp)))+sum(sum(bleach_extent(:,select_sp)));
    count = count+species_ID(s);
end    
    
