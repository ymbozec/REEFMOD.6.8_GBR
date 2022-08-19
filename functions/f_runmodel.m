% -------------------------------------------------------------------------
% Y.-M. Bozec, MSEL, created Nov 2011.
% Last modified: 08/2022
%
% REEFMOD RUN FILE
% -------------------------------------------------------------------------

% This function runs 1 simulation of ReefMod, either for a single reef or multiple reefs
% The number of reefs is specified in META.nb_reefs (default is 1).
%
% Coral and algal covers for reef #n are stored respectively in the struct arrays:
% - metapop(n).coral
% - metapop(n).algal
% where ~.coral is itself a struct array that gives for every coral species 's':
%       - ~.coral(s).cover_cm2 = the 2D planimetric area (coral size) in cm2 of every colony in every cell
%               -> positive values for live corals
%               -> negative values for dead corals
%       - ~.coral(s).surface_cm2 = the 3D surface area in cm2 (paraboloid) of every colony in every cell
%       - ~.coral(s).volume_cm3 = the volume in cm3 (paraboloid) of every colony in every cell
%       - ~.coral(s).clade = clade 1 (thermally sensitive) or 2 (tolerant)
% and where ~.algal is itself a struct array that gives for every algal type 'a':
%       - ~.algal(a).cover_cm2 = the size (planimetric area in cm2) of their patch in every cell
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [RESULT, RECORD] = f_runmodel(META, REEF, CORAL, ALGAL, CONNECT, REEF_POP, REEF_COTS)

%___________________________________________________________________________________________________
%%
%       SPACE ALLOCATION
%___________________________________________________________________________________________________

metapop(META.nb_reefs).coral = {} ;
metapop(META.nb_reefs).algal = {} ;
metapop(META.nb_reefs).genes = {} ;

if META.doing_bleaching == 1   
    load('BleachingModelHughes.mat')       
end

% 'RESULT' stores the outcome of variables over the grid at the end of each time step (ie, t+1) for every reef
RESULT.coral_pct2D = zeros(META.nb_reefs, META.nb_time_steps+1, META.nb_coral_types);
RESULT.algal_pct = zeros(META.nb_reefs, META.nb_time_steps+1, META.nb_algal_types);
RESULT.coral_settler_count = zeros(META.nb_reefs, META.nb_time_steps+1, META.nb_coral_types);
RESULT.coral_total_fecundity = zeros(META.nb_reefs, META.nb_time_steps+1, META.nb_coral_types);
RESULT.coral_larval_supply = zeros(META.nb_reefs, META.nb_time_steps+1, META.nb_coral_types);

if META.doing_clades == 1
    RESULT.clade_prop = zeros(META.nb_reefs, META.nb_time_steps+1, META.nb_coral_types);
end

if META.doing_3D==1
    RESULT.live_CaCO3 = zeros(META.nb_reefs, META.nb_time_steps+1, META.nb_coral_types);
    RESULT.dead_CaCO3 = zeros(META.nb_reefs, META.nb_time_steps+1, META.nb_coral_types);
    RESULT.substrate_cm2 = zeros(META.nb_reefs, META.nb_time_steps+1);
end

if META.doing_size_frequency==1
    RESULT.coral_juv_count = zeros(META.nb_reefs, META.nb_time_steps+1, META.nb_coral_types, length(0:CORAL.diam_bins(1):CORAL.juv_max_diam)-1);
    RESULT.coral_adol_count = zeros(META.nb_reefs, META.nb_time_steps+1, META.nb_coral_types, length(CORAL.juv_max_diam:CORAL.diam_bins(2):CORAL.adol_max_diam)-1);
    RESULT.coral_adult_count = zeros(META.nb_reefs, META.nb_time_steps+1, META.nb_coral_types, length(CORAL.adol_max_diam:CORAL.diam_bins(3):(CORAL.adult_max_diam+CORAL.diam_bins(3)))-1);
end

if META.doing_COTS==1
    RESULT.COTS_all_densities = zeros(META.nb_reefs, META.nb_time_steps+1, META.COTS_maximum_age);
    RESULT.COTS_total_perceived_density = zeros(META.nb_reefs, META.nb_time_steps+1);
    RESULT.COTS_settler_densities = zeros(META.nb_reefs, META.nb_time_steps+1); % density of settlers before mortality (then become 6mo old recruits)
    RESULT.COTS_larval_supply = zeros(META.nb_reefs, META.nb_time_steps+1);
    RESULT.COTS_larval_output = zeros(META.nb_reefs, META.nb_time_steps+1);  % amount of larvae produced per reef (fecundity*area*survival)
    RESULT.COTS_fecundity = zeros(META.nb_reefs, META.nb_time_steps+1); % amount of larvae per grid (direct output of adult density)
    RESULT.COTS_adult_densities = zeros(META.nb_reefs, META.nb_time_steps+1);   
    % Note: outputs of CoTS control are all in RESULT.COTS_records (created in the CoTS control module)
end

if META.tracking_rubble ==1
    RESULT.rubble_cover_pct2D = zeros(META.nb_reefs, META.nb_time_steps+1) ;
    RESULT.rubble_cover_pct2D(:,1)= cat(1,REEF(1:META.nb_reefs).initial_rubble_pct) ; 
end


if META.doing_genetics == 1
    REEF(META.nb_reefs).scale_mutation_rate = zeros(1,1);
    RESULT.genetic_diversity = zeros(META.nb_reefs, META.nb_time_steps+1, META.nb_coral_types);
    RESULT.relative_fitness = zeros(META.nb_reefs, META.nb_time_steps+1, META.nb_coral_types);
    Topt_bins = [22:0.2:32] ;
    RESULT.Topt_distri = zeros(META.nb_reefs, META.nb_time_steps+1, META.nb_coral_types, length(Topt_bins));
    RESULT.Topt_mean = zeros(META.nb_reefs, META.nb_time_steps+1, META.nb_coral_types);
else
    META.Topt2index = 0;
end

% 'RECORD' stores events occuring at each time step (i.e., t)
RECORD.applied_DHWs = zeros(META.nb_reefs, META.nb_time_steps,'single') ;
bleaching_mortalities = zeros(META.nb_reefs, META.nb_time_steps,'single') ; % internal variable to pre-determine bleaching mortality
RECORD.applied_bleaching_mortality = zeros(META.nb_reefs, META.nb_time_steps,'single') ; % actual mortality given the DHW trehsold and cyclone occurence
RECORD.hurricane_events = uint8(zeros(META.nb_reefs, META.nb_time_steps)) ;
RECORD.hurricane_parm_k = zeros(META.nb_reefs, META.nb_time_steps);

if META.doing_restoration == 1   
    RECORD.total_outplanted = zeros(META.nb_reefs, META.nb_time_steps, length(META.outplant_species)) ;
    RECORD.outplanted_reefs = uint16(zeros(META.nb_reefs, META.nb_time_steps+1)) ;
    RECORD.total_enriched = zeros(META.nb_reefs, META.nb_time_steps, length(META.enriched_species)) ;
    RECORD.enriched_reefs = uint16(zeros(META.nb_reefs, META.nb_time_steps+1)) ;
    RECORD.rubble_cover_pct2D_stabilised = zeros(META.nb_reefs, META.nb_time_steps) ;
    RECORD.stabilised_reefs = uint16(zeros(META.nb_reefs, META.nb_time_steps+1)) ; 
    RECORD.fogged_reefs = uint16(zeros(META.nb_reefs, META.nb_time_steps+1)) ; 
    RESULT.coral_pct2D_restored_sites = RESULT.coral_pct2D;
    never_deployed = ones(META.nb_reefs,1);
end

% Track percentage cover loss per species for ALL stressors (03 Feb 2020)
RECORD.coral_pct2D_lost_bleaching = zeros(META.nb_reefs, META.nb_time_steps, META.nb_coral_types) ;
RECORD.coral_pct2D_lost_cyclones = zeros(META.nb_reefs, META.nb_time_steps, META.nb_coral_types) ;
RECORD.coral_pct2D_lost_COTS = zeros(META.nb_reefs, META.nb_time_steps, META.nb_coral_types) ;

% Internal variables
ID_colony_tracking = ones(META.nb_reefs,META.nb_coral_types);
colony_list = zeros(1,6);
environ_list = zeros(1,5);
convert_area_ha = (10^8)/META.total_area_cm2;

if META.doing_coral_connectivity == 0   
    META.area_habitat = ones(META.nb_reefs,1); % needs to be defined if connectivity is OFF
end

%______________________________________________________________________________________________________________________________________________
%%
% INITIALISATION
%______________________________________________________________________________________________________________________________________________
% disp('Initialising reef states...')

nb_cells = META.grid_x_count * META.grid_y_count ;
list_cell = uint32(1:nb_cells) ;

t = 0 ;

for n = 1:META.nb_reefs % This must be done for every reef before time simulations
    
    
    %%%%%% RE-INITIALISE GROWTH RATE IN CASE IT PREVIOUSLY CHANGED %%%%%%
    if META.doing_calcification == 1
        % Re-initialise coral growth rates because they may have changed
        CORAL.growth_rate = init_coral_growth_rate ;
    end
    
    %%%%%% DEFINE UNGRAZABLE CELLS FIRST %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    REEF(n).grazable_cell = ones(nb_cells,1) ; % create a list of grazable cells (grazable==1)
    id = randsample(list_cell, uint32(REEF(n).nongrazable_substratum * nb_cells)) ; % select randomly cell indices
    REEF(n).grazable_cell(id) = 0 ; % switch to 0 the randomly selected cells (nongrazable==0)
    
    % Define the surface area (SA) of the reef substrate of each grazable cell
    % By default, every cell is flat. If doing 3D, the substrate SA will change over time,
    % otherwise will stay flat. Non-grazable cells (sand) are flat (SA = cell area)
    REEF(n).substrate_SA_cm2 = ones(nb_cells,1)*META.cell_area_cm2;
    REEF(n).floor_SA_cm2 = ones(nb_cells,1)*META.cell_area_cm2;
    
    %%%%%% INITIALISE THE GRID WITH CORALS AND ALGAE %%%%%%%%%%%%%%%%%%%%%%
    [metapop(n).coral, metapop(n).algal] = f_initialise_population(META, REEF(n), CORAL) ;
        
    %%%%%% SET UP 3D %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if META.doing_3D == 1
        % Estimate surface, volume of living and dead (after erosion) colonies.
        % Produces also different surface areas for each cell in REEF that will vary over time:
        % - surface area of the substrate underneath live corals (REEF.substrate_SA_cm2)
        % - surface area of the reef floor, ie including living colonies (REEF.floor_SA_cm2)
        [metapop(n).coral, metapop(n).algal, REEF] = f_initialise_rugosity(META, REEF(n), CORAL, metapop(n).coral, metapop(n).algal, ALGAL(n).initial_cover);
        
        RESULT.volume_eroded(n,t+1) = 0 ;
        RESULT.substrate_cm2(n,t+1) = sum(REEF(n).substrate_SA_cm2); % total 3D area of the substrate
        RESULT.floor_cm2(n,t+1) = sum(REEF(n).floor_SA_cm2); % total 3D area of the sea floor (including coral 3D surface)
    end
    
    %%%%%% SET UP COTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Reminder: t=0 at initialization
    if META.doing_COTS == 1
        
%         if REEF_COTS.densities_M(n,t+1) > -1 % if empirical estimate available
        if isnan(REEF_COTS.densities_M(n,t+1))==0  % if empirical observation is available
            
            switch META.randomize_initial_COTS_densities
                
                case 1 % Gaussian generated from specified mean and SD (use SD=0 to force to the mean)
                    
                    COTS_total_density = normrnd(REEF_COTS.densities_M(n,1),REEF_COTS.densities_SD(n,1));
                    COTS_total_density(COTS_total_density<0)=0;
                    
                case 2 % Poisson generated to avoid creating low densities (relative to SD)
                    
                    COTS_total_density = poissrnd(REEF_COTS.densities_M(n,1));
                    
            end
            
        else % otherwise we initialise with 0 CoTS
            
            COTS_total_density = 0 ;
            
        end
        
        % Initialise with reference pop structure in winter (already corrected for imperfect detection)
        % Density is number of CoTS per reef grid (400m2)
        COTS_densities = COTS_total_density * META.COTS_init_age_distri_OUTBREAK(2,:) ;
        mature_COTS_density = sum(COTS_densities(:,META.COTS_fecundity~=0),2);
        
        % Fertilization success (0-1) from number of mature CoTS per hectare (Babcock et al. 2014)
        fertilization_success = 0.14 * (convert_area_ha * mature_COTS_density).^0.61 ;
        fertilization_success(fertilization_success>0.9)=0.9;  %fertilization cap (max obtained by Babcock)
        
        % Population estimates at initialization (t=0)
        RESULT.COTS_total_perceived_density(n,t+1) = sum(COTS_densities(META.COTS_adult_min_age:end).*META.COTS_detectability(META.COTS_adult_min_age:end));
        RESULT.COTS_all_densities(n,t+1,:) = COTS_densities ;
        RESULT.COTS_fecundity(n,t+1) = sum(COTS_densities.*META.COTS_fecundity).*fertilization_success ;
        
    end
    
    %%%%%% SET UP CYCLONES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % If randomized with the chosen regime (makes the chronology starting at different dates)
    if META.doing_hurricanes == 1
        
        if META.randomize_hurricane_chronology == 1
            
            [c]=size(REEF(n).hurricane_chronology,2); % c is the number of time steps (seasons) of the chronology
            date = randsample(c,1);% Pick up a date randomly for starting the new chronology
            % Cut the chronology in two pieces at the chosen date, and collate the left piece after the right piece
            REEF(n).hurricane_chronology = [REEF(n).hurricane_chronology(:,date:c) REEF(n).hurricane_chronology(:,1:c-1)] ;
        end
        
        % March 2022: pre-determination of the constant k of the mortality rate with or without noise.
        % Was previously in f_hurricane_effect.m. See script for definition of parameters
        a = -3e-007;
        k_cat5 = 0.0551-a*((0.0007/(-2*a))^2);
        relHurrMort = [0.0461 0.1183 0.2504 0.5677 1]; %e.g. a cat 4 has 57% of the impact as a cat 5
        
        for t0=1:2:META.nb_time_steps
            
            hurricaneCat = REEF(n).hurricane_chronology(t0);
            
            if hurricaneCat>0
                
                if META.deterministic_hurricane_mortality == 1
                    
                    RECORD.hurricane_parm_k(n,t0) = k_cat5*relHurrMort(hurricaneCat);
                else
                    % NEW (Jul 2019): add Gaussian noise into the calculation of mortality)
                    RECORD.hurricane_parm_k(n,t0) = k_cat5*(normrnd(relHurrMort(hurricaneCat),0.1));
                end
            end
        end
    end
   
    %%%%%% SET UP BLEACHING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if META.doing_bleaching == 1
        
        bleaching_mortalities(n,:) = f_generate_bleaching_mortalities(BleachingLinearModel,REEF(n).predicted_DHWs, META.deterministic_bleaching) ;
        bleaching_mortalities(n,2:2:end) = 0; % no bleaching in winter
    end
    
    %%%%%% STORE INITIAL COVERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Reminder: t=0 at initialization

    % 1) Algal cover: total cover (%) of each algal type over the grid at initial step
    RESULT.algal_pct(n,t+1,1:META.nb_algal_types) = 100*full(sum([metapop(n).algal(1:META.nb_algal_types).cover_cm2]))./sum(REEF(n).substrate_SA_cm2) ;
    
    % 2) Coral cover: total cover (%) of each coral type at initial step
    for s=1:META.nb_coral_types
        
        living_planar_coral_cover_cm2 = metapop(n).coral(s).cover_cm2(metapop(n).coral(s).cover_cm2>0) ;
        RESULT.coral_pct2D(n,t+1,s) = 100*sum(living_planar_coral_cover_cm2)./sum(REEF(n).substrate_SA_cm2) ;
        
        if META.doing_clades == 1
            living_coral_clades = metapop(n).coral(s).clade(metapop(n).coral(s).cover_cm2>0) ;
            RESULT.clade_prop(n,t+1,s) = nnz(living_coral_clades==1)/nnz(living_coral_clades);
        end
        
        % Estimaton of living and dead carbonate volumes (per species)
        if META.doing_3D == 1
            living_coral_volume_cm3 = metapop(n).coral(s).volume_cm3(metapop(n).coral(s).cover_cm2>0) ;
            dead_coral_volume_cm3 = metapop(n).coral(s).volume_cm3(metapop(n).coral(s).cover_cm2<0) ;
            RESULT.live_CaCO3(n,t+1,s)= sum(living_coral_volume_cm3);
            RESULT.dead_CaCO3(n,t+1,s)= sum(dead_coral_volume_cm3);
        end
        
        % Coral size structure at initial step
        if META.doing_size_frequency == 1
            [count, ~] = f_count_sizefreq(living_planar_coral_cover_cm2, CORAL) ;
            RESULT.coral_juv_count(n,t+1,s,:)=count.juv;
            RESULT.coral_adol_count(n,t+1,s,:)=count.adol;
            RESULT.coral_adult_count(n,t+1,s,:)=count.adult;
        end
        
        % Record the ID max
        if sum(sum(metapop(n).coral(s).cover_cm2))~=0
            ID_colony_tracking(n,s) = max(max(metapop(n).coral(s).colony_ID));
        end
    end
    
    % 3) Estimate reef rugosity
    if META.doing_3D == 1
        coral_sp_pct2D = squeeze(squeeze(RESULT.coral_pct2D(n,t+1,:)));
        [rugosity, SI_reef] = f_estimate_rugosity(REEF(n), META, CORAL, coral_sp_pct2D);
        RESULT.rugosity(n,t+1)=rugosity;
        RESULT.SI_reef(n,t+1)=SI_reef;
    end
    
    if META.track_populations == 1
        [environ_list,colony_list] = f_track_populations(environ_list,colony_list,metapop(n).algal,metapop(n).coral,t,META,REEF(n).grazable_cell);
    end
    
    % 4) Create genotypes
    if META.doing_genetics == 1
        
        REEF(n).scale_mutation_rate = META.area_habitat(n)/(1e-10 * META.total_area_cm2);
        
        for s = 1:META.nb_coral_types
            
            metapop(META.nb_reefs).genes(s).locus = {} ;
            
            if META.genetics.group(s)==1
                            
                list_coral_ID = metapop(n).coral(s).colony_ID(metapop(n).coral(s).colony_ID~=0);
                idx = randi(META.genetics_pop_size,length(list_coral_ID),1);
                
                metapop(n).genes(s).QTLs = REEF(n).coral(s).QTL_pool_IN(idx,:,:);
                metapop(n).genes(s).list_coral_ID = list_coral_ID ;
                
                % Compute genotypes with heritability
                % First calculate breeding value
                BREED = sum(sum(metapop(n).genes(s).QTLs,3),2) ;
                % Then determine phenotype following heritability (esd=0 implies perfect heritability)
                Env_effect = normrnd(0, META.genetics.esd(s), size(BREED));
%                 Env_effect(Env_effect>4) = 4; % Limit the environmental effect (add no more that +6 deg C to BREED)
                
                metapop(n).genes(s).phenotypes = BREED + Env_effect ;
                
                max_phenotype = floor(max(META.List_Topt)-REEF(n).Topt_baseline);
                
                metapop(n).genes(s).phenotypes(metapop(n).genes(s).phenotypes > max_phenotype) = max_phenotype ;
            
                All_Topts = round(10*(REEF(n).Topt_baseline + metapop(n).genes(s).phenotypes))/10 ;
                    
                [count_Topt, ~] = hist(All_Topts,Topt_bins) ;
                RESULT.Topt_distri(n,t+1,s,:) = count_Topt;
                RESULT.Topt_mean(n,t+1,s) = mean(All_Topts);
                
            end
        end
    else
        REEF(n).coral.QTL_pool_IN=[];
    end
    
    
    % Estimate coral fecundity of the reef with associated genetics
    for s=1:META.nb_coral_types
        
        find_colonies =  find(metapop(n).coral(s).cover_cm2>0);
        
        if isempty(find_colonies)==0
            
            F_list = ones(length(metapop(n).coral(s).cover_cm2(find_colonies)),1);
            % (everyone is fit at initial step)
            
            % Number of larvae produced by the reef (per 400m2 of reef)
            [RESULT.coral_total_fecundity(n,t+1,s)] = f_estimate_fecundity(metapop(n).coral(s).cover_cm2(find_colonies), F_list,...
                CORAL.fecund_min_size(s), CORAL.fecund_a(s), CORAL.fecund_b(s));
                        
            if META.doing_genetics == 1
                
                % List the potential parents
                list_coral_ID = metapop(n).genes(s).list_coral_ID ;
                  
                % Genotype of parents
                QTLs = metapop(n).genes(s).QTLs ;
                
                select1 = randsample(1:length(list_coral_ID), META.genetics_pop_size, 'true', F_list')'; % pick up parent 1 for each larva
                select2 = randsample(1:length(list_coral_ID), META.genetics_pop_size, 'true', F_list')'; % pick up parent 2 for each larva (note this does not avoid selfing)
                
                % Look for hermaphrodites
                test_self =  select1 - select2;
                select2(test_self==0) = randi(length(list_coral_ID), length(test_self(test_self==0)),1);
                % Note this won't completely exclude self-fertilization, which would become increasingly prevalent
                % when local population is getting reduced to few individuals (in which case the number of larvae
                % should be relatively low so that the impact of self-fertiization on the network shoud be negligible)
                
                % Create a representative pool of larval genotypes produced by the reef
                REEF(n).coral(s).QTL_pool_OUT(1:META.genetics_pop_size,:,1) = squeeze(QTLs(select1,:,1));
                REEF(n).coral(s).QTL_pool_OUT(1:META.genetics_pop_size,:,2) = squeeze(QTLs(select2,:,1));
                
                % Record the mutation rate on each source reef
                % Mutation rate is scaled to the estimated population size, which is nb colonies * habitat area
                LocalMuteRate = REEF(n).scale_mutation_rate*META.genetics.MuteRate(s)*ID_colony_tracking(n,s)/META.genetics_pop_size;
                rand_MuteIncidence = rand(META.genetics_pop_size, META.genetics.nb_loci, 2)<LocalMuteRate; % rand_MuteIncidence is logical!
                mutations = zeros(size(rand_MuteIncidence));
                
                mutations(rand_MuteIncidence==1) = normrnd(0, META.genetics.MuteEffect(s),sum(sum(sum(rand_MuteIncidence))),1);
                REEF(n).coral(s).QTL_pool_OUT = REEF(n).coral(s).QTL_pool_OUT + mutations ;
                
                % Update the pool of genotypes coming in = due to retention
                % Note this pool will be updated at the next step with larvae coming from other reefs
                REEF(n).coral(s).QTL_pool_IN = REEF(n).coral(s).QTL_pool_OUT ;
                
            end            
        end
    end   
end

META.genetics.QTL_pool =[]; %don't need this anymore


%______________________________________________________________________________________________________________________________________________
%%
% RUN TIME SIMULATIONS
%______________________________________________________________________________________________________________________________________________

arrange_step = zeros(1,META.nb_time_steps);
arrange_step(1:4:META.nb_time_steps)=1;
arrange_step(1)=0;

if META.randomize_WQ_chronology == 1
    % Randomize the pick-up of WQ and connectivity layers
    RECORD.WQ_chronology = randi(6,1,META.nb_time_steps);
else
% OLD HINDCAST:
%     WQ_chronology_tmp = repmat([6 6 7 7 8 8 1 1 2 2 3 3 4 4 5 5 6 6],1,META.nb_time_steps); % too long but doesn't matter here
% NEW FOR HINDCAST & FORECAST:
% Cycle the 2011-2018 water years twice but with 2011 (extreme wet year) applied once, then repeat in the future. This ensures that
% the extreme wet season is applied once every 15 years (2011, 2026, 2041, ...) to match historical outbreak patterns (Pratchett et al. 2014)
    WQ_chronology_tmp = repmat([6 6 7 7 8 8 1 1 2 2 3 3 4 4 5 5 6 6 7 7 8 8 2 2 3 3 4 4 5 5],1,20); % too long but doesn't matter here

    RECORD.WQ_chronology = WQ_chronology_tmp(1:META.nb_time_steps);   
end


for t = 1:META.nb_time_steps
    
    select_layer = RECORD.WQ_chronology(t) ;

    season = iseven(t); % 1 is winter, 0 is summer (first time step is summer)
    
    if season == 0 % reproductive season is summer
        
        %% %%%%%%% ESTIMATE CORAL LARVAL SUPPLY FOR EVERY REEF %%%%%%%%%%%%%%%%%
        output_larvae = squeeze(RESULT.coral_total_fecundity(:,t,:).*META.area_habitat); %larvae produced per 400m2 'scaled up' to reef area
        % Formally this should be multiplied by 2500 but simplified to avoid carrying on large numbers (would be divided by 2500 afterwards any way) 
        
        if META.doing_water_quality == 1
            % Suspended sediments reduce fertilization and the number of competent larvae
            output_larvae = output_larvae.*REEF_POP(select_layer).CORAL_larvae_production(:,1) ;
        end
        
        if META.doing_coral_connectivity == 1
            
            coral_conn_matrix = CONNECT(select_layer).ACROPORA; % Only one matrix for all species
            
            % If only a subset of reefs is simulated, need to account for larvae coming from surrounding reefs
            if isempty(META.outside_reef_ID)==0
                in_degree_larval_prop = CONNECT(select_layer).ACROPORA_ext';
                % Assume the pool of larvae in the greater region can be approximated by the average pool of larvae of the focal region (per reef),
                % then distribute according to in-degree edges
                BOUNDARY_CONNECTIVITY = mean(output_larvae,1) .* in_degree_larval_prop(:,ones(1,META.nb_coral_types));
            else
                BOUNDARY_CONNECTIVITY = 0;
            end
            
            if META.nb_reefs>1
                
                % This is the number of larvae supplied per species per unit reef area (400m2)
                RESULT.coral_larval_supply(:,t,:) = (((output_larvae')*coral_conn_matrix)' + BOUNDARY_CONNECTIVITY)./META.area_habitat ;              
            else     
                RESULT.coral_larval_supply(:,t,:) = ((output_larvae*coral_conn_matrix)' + BOUNDARY_CONNECTIVITY)./META.area_habitat ;             
            end
            
        else
            
            % Fixed regime of external larval supply
            RESULT.coral_larval_supply(:,t,:) = (META.coral_min_selfseed*output_larvae./META.area_habitat) + squeeze(META.coral_immigration(:,t,:));
        end
        
        %% Restoration: larval enrichment
%         %% INACTIVATED FOR NOW (MARCH 2022): moving larvae simulated as a deployment of coral outplants (1 yr old)
%         if META.doing_restoration == 1
%             
%             if META.nb_reefs_enriched > 0
%                 
%                 if META.doing_larval_enrichment(t)==1
%                     
%                     All_reef_states = squeeze(sum(RESULT.coral_pct2D(:,t,:),3));
%                     
%                     Priority_reef_ID = META.priority_list_LarvalEnrich;
%                     
%                     % Exclude from priority list those reefs with large amount of corals (no need to enrich larval supply)
%                     Priority_reef_ID(All_reef_states(META.priority_list_LarvalEnrich)>META.threshold_for_larval_enrichment)=[];
%                     
%                     if isempty(Priority_reef_ID)==0
%          
%                         nb_restored_reefs = META.nb_reefs_enriched ;
%                         nb_restored_reefs(nb_restored_reefs>length(Priority_reef_ID)) = length(Priority_reef_ID);
%                         
%                         extra_larvae = META.enriched_larvae(ones(1,nb_restored_reefs),:);
%                         X=reshape(extra_larvae,size(extra_larvae,1),1,size(extra_larvae,2));
%                         
%                         RESULT.coral_larval_supply(Priority_reef_ID(1:nb_restored_reefs),t,:) = ...
%                             RESULT.coral_larval_supply(Priority_reef_ID(1:nb_restored_reefs),t,:) + X ;
%                         
%                         RECORD.enriched_reefs(Priority_reef_ID(1:nb_restored_reefs),t)=1;
%                         
%                     end
%                 end
%             end
%         end
            
        
        %% %%%%%%% ESTIMATE COTS LARVAL SUPPLY %%%%%%%%%%%%%%%%%
        if META.doing_COTS == 1
            % First weigh output larvae by suitable habitat for corals
            COTS_larval_output = RESULT.COTS_fecundity(:,t).*META.area_habitat;
            
            if META.doing_water_quality == 1 && META.doing_Chl_forcing == 1
                % mean Chlorophyll concentration (sum of Chl a) increases larval survivorship
                RESULT.COTS_larval_output(:,t) = COTS_larval_output.*REEF_POP(select_layer).COTS_larvae_survival(:,1) ;
            else
                RESULT.COTS_larval_output(:,t) = COTS_larval_output.*META.COTS_min_larval_survival ; % without WQ default survival
            end
            
            if META.doing_COTS_connectivity == 1
                
                cots_conn_matrix = CONNECT(select_layer).COTS;
                RESULT.COTS_larval_supply(:,t) = (RESULT.COTS_larval_output(:,t)')*cots_conn_matrix./META.area_habitat' ;
                
            else
                % Self recruitment with fixed regime of COTS immigration (= no connectivity)
                RESULT.COTS_larval_supply(:,t) = (META.COTS_min_selfseed*RESULT.COTS_larval_output(:,t)./META.area_habitat) + META.COTS_immigration(1,t);
            end            
        end
    end
    
    %% %%%%%%% PROCESS EVERY REEF ONE AFTER THE OTHER  %%%%%%%%%%%%%%%%%%%%%
%     seed = rng;

    for n = 1:META.nb_reefs % This must be done for every reef before time simulations
             
        %%%% --------------------------------------------------------------
        %%%% Water quality
        %%%% --------------------------------------------------------------
        if META.doing_water_quality == 1

            REEF(n).CORAL_recruit_survival = REEF_POP(select_layer).CORAL_recruit_survival(n,1).*ones(1,META.nb_coral_types);
            REEF(n).CORAL_recruit_survival(1,4:6) = 1; % Full survival for non-Acropora groups
            REEF(n).CORAL_juvenile_growth = REEF_POP(select_layer).CORAL_juvenile_growth(n,season+1)*CORAL.juvenile_growth_rate;
        else            
            REEF(n).CORAL_recruit_survival = ones(1,META.nb_coral_types);
            REEF(n).CORAL_juvenile_growth = CORAL.juvenile_growth_rate;
        end
        
        %%%% --------------------------------------------------------------
        %%%% Genetics of coral larvae (recruitment occurs at the end of the time step)
        %%%% --------------------------------------------------------------
        if season == 0 && META.doing_genetics == 1 % Only recruit in summer
            
            for s = 1:META.nb_coral_types
                
                if RESULT.coral_larval_supply(:,t,s)==0
                    continue
                else
                    
                    if META.nb_reefs > 1
                        larval_origins = floor(1000*coral_conn_matrix(:,n).*output_larvae(:,s)/sum(coral_conn_matrix(:,n).*output_larvae(:,s)));
                    else
                        larval_origins = 1000;
                    end
                    
                    all_reefs = 1:META.nb_reefs;
                    select_sources = all_reefs(larval_origins>0); %ID of source reefs
                    select_sources = select_sources(select_sources~=n); %(excluding reef n)
                    
                    if isempty(select_sources)==0
                        
                        c = 1;
                        
                        % Update the pool of genotypes produced by reef n with a random selection comming
                        % from the source reefs (they replace the existing ones
                        for nn = 1:length(select_sources) % for every source reef
                            
                            Nb_larvae = larval_origins(select_sources(nn)); % nb of larval genotypes to create for this source reef
                            
                            REEF(n).coral(s).QTL_pool_IN(c:(c+Nb_larvae-1),:,:) = ...
                                REEF(select_sources(nn)).coral(s).QTL_pool_OUT(randi(META.genetics_pop_size, Nb_larvae,1),:,:);
                            
                            c = Nb_larvae+1;
                        end
                    end
                end
            end
        end
        
        
        %%%% --------------------------------------------------------------
        %%%% POPULATION DYNAMICS
        %%%% --------------------------------------------------------------
        
        %%%%%%%%% EFFECTS OF SST ON CORAL GROWTH RATE %%%%%%%%%%%%%%%%%%%%%
        % NOTE: use the pre-specified scenario of growth rates as in Bozec & Mumby (2015)
        if META.doing_calcification == 1
            % Predicted SST for the current time step
            current_sst = REEF(n).SST(1,t);  % sd_sst = 1.47 % useful ??
            % Then calculate the relative reduction of calcification rate for each coral species
            Rel_change = exp( - 0.5 *( (current_sst - CORAL.SST_OPT) ./ CORAL.SD_relative_calci) .^2 );
            % Applies change to linear extension
            CORAL.growth_rate = init_coral_growth_rate .* CORAL.relative_calci .* Rel_change;
        end
        
        %%%%%%%%% EFFECTS OF RUBBLE ON CORAL JUVENILES
        if META.tracking_rubble==1
            % Elevate juvenile mortality with the amount of rubble present (RRAP feasibility study)
%             REEF(n).juv_whole_mortality_rate = CORAL.juv_whole_mortality_rate + ...
%                 (1-CORAL.juv_whole_mortality_rate) * RESULT.rubble_cover_pct2D(n,t)/100 ;        

            REEF(n).CORAL_recruit_survival = REEF(n).CORAL_recruit_survival*(1-RESULT.rubble_cover_pct2D(n,t)/100) ;
        end
        
        %%%%%%%%% THEN PROCESS CORAL POPULATION OF REEF N %%%%%%%%%%%%%%%%%
        % 11/07/2016: First estimate the balance amount of grazing needed (previously ALGALREMOVAL)
        algal_removal = (REEF(n).diadema).*ALGAL.diadema_props + (REEF(n).herbivory(1,t)).*ALGAL.herbivory_props(:,t);
        
        % Calculate difference between current SST and the Topt baseline of coral's phenotype for that reef
        SST_diff = REEF(n).predicted_SST(t)-REEF(n).Topt_baseline ; 
%         seed=rng;
        [metapop(n).coral, metapop(n).algal, metapop(n).genes, last_surface_area_grazed] = ...
            f_process_population (metapop(n).coral, metapop(n).algal, metapop(n).genes, season, algal_removal, META, REEF(n), CORAL, ALGAL, SST_diff);
% rng(seed);
        
        %%%% --------------------------------------------------------------
        %%%% CROWN-OF-THORN STARFISH PREDATION
        %%%% --------------------------------------------------------------
        % COTS population at t was previously estimated by the model
        % (either initialized for t=0 or processed through demographics
        if META.doing_COTS == 1
            
            % Force outbreaks to stop after several years, otherwise outbreaks would only stop when there is no coral left
            % while COTS populations could crash due to diseases (Pratchett 2001, Zann et al. 1990)
            % If META.COTS_outbreak_duration is a vector of two scalars, outbreak duration is chosen at random between those two bounds;
            % If META.COTS_outbreak_duration is a scalar then outbreak duration is fixed to that value
            % (see settings_COTS.m)
            
            % First tests if at the start of time step t the density of 18+
            % months COTS is above threshold for disease (ie, above outbreak threshold)
            if RESULT.COTS_adult_densities(n,t) >= META.COTS_density_threshold_for_disease
                    
                % If that the case, track recent history of that density
                if length(META.COTS_outbreak_duration)>1
                    COTS_outbreak_duration = randi(META.COTS_outbreak_duration); % in years
                else
                    COTS_outbreak_duration = META.COTS_outbreak_duration;
                end
                
                % Process die-back if outbreak has lasted the max duration
                if t > 2*COTS_outbreak_duration+1
                    
%                     COTS_density_history = sum(RESULT.COTS_all_densities(n,(t-COTS_outbreak_duration*2):t,3:end),3);
                    COTS_density_history = RESULT.COTS_adult_densities(n,(t-COTS_outbreak_duration*2):t);
                    
                    if min(COTS_density_history)>= META.COTS_density_threshold_for_disease
                        % if COTS have been above the outbreak density threshold over the tested period
                        % then population crashes due to disease - just leave CoTS recruits
                        RESULT.COTS_all_densities(n,t, META.COTS_dieback_classes)=REEF(n).COTS_background_density * ...
                            RESULT.COTS_all_densities(n,t,META.COTS_dieback_classes)/sum(RESULT.COTS_all_densities(n,t,META.COTS_dieback_classes)) ;
                        % RESULT.COTS_larval_supply(n,t) = 0 ; % prevent settlement if disease
                    end
                end
            end
                         
            % If forcing of CoTS density is available, then erase the predicted pop
%             if t > 1 && REEF_COTS.densities_M(n,t) > -1 % if empirical estimate available (otherwise do nothing)
%                 %(starts at t=2 to not override initial conditions)
            if t> 1 && isnan(REEF_COTS.densities_M(n,t))==0  % if empirical observation is available (otherwise do nothing)
                %(starts at t=2 to not override initial conditions)
                
                COTS_total_density = normrnd(REEF_COTS.densities_M(n,t),REEF_COTS.densities_SD(n,t));
                % COTS_total_density = poissrnd(25*REEF_COTS.densities_M(n,t))/25;
                COTS_total_density(COTS_total_density<0)=0;
                % COTS_densities = COTS_total_density * META.COTS_init_age_distri_OUTBREAK(season+1,:) ;
                COTS_densities = COTS_total_density * META.COTS_init_age_distri_OUTBREAK(2-season,:) ;
                
                RESULT.COTS_all_densities(n,t,:) = COTS_densities ;
                RESULT.COTS_total_perceived_density(n,t) = sum(COTS_densities(1,META.COTS_adult_min_age:end).*META.COTS_detectability(META.COTS_adult_min_age:end));
                RESULT.COTS_adult_densities(n,t) = sum(RESULT.COTS_all_densities(n,t,META.COTS_adult_min_age:end),3); % for CoTS control
                
                % Note: we populate at t because this is CoTS population before the processing of mortality
            end
                 
            
            % Add COTS recruits before processing CoTS pop dyn
            if season == 0 % Only recruit in summer
                RESULT.COTS_settler_densities(n,t) = META.COTS_BH_alpha*(RESULT.COTS_larval_supply(n,t)/(META.COTS_BH_beta + RESULT.COTS_larval_supply(n,t)));
            else
                RESULT.COTS_settler_densities(n,t) = 0;
            end
            
            
            % Then process COTS dynamics (mortality) and predation
            COTS_feeding_rates = META.COTS_feeding_rates(season+1,:); % current feeding rates (seasonal)
                           
            [metapop(n).coral, metapop(n).algal, NEW_COTS_densities ,total_coral_loss_COTS] = f_apply_COTS_predation(metapop(n).coral, metapop(n).algal,...
                RESULT.COTS_all_densities(n,t,:), RESULT.COTS_settler_densities(n,t), COTS_feeding_rates, CORAL.SI, REEF(n).COTS_background_density, META) ;
                       
            mature_COTS_density = sum(NEW_COTS_densities(1,META.COTS_fecundity~=0),2) ;
            
            % Fertilization success (0-1) from number of mature CoTS per hectare (Babcock et al. 2014)
            fertilization_success = 0.14 * (convert_area_ha*mature_COTS_density).^0.61 ;
            fertilization_success(fertilization_success>0.9)=0.9;  %fertilization cap (max obtained by Babcock)
            
            % Population estimates at the end of t (ie, t+1)
            RESULT.COTS_total_perceived_density(n,t+1) = sum(NEW_COTS_densities(1,META.COTS_adult_min_age:end).*META.COTS_detectability(META.COTS_adult_min_age:end));
            RESULT.COTS_all_densities(n,t+1,:) = NEW_COTS_densities ;
            RESULT.COTS_adult_densities(n,t+1) = sum(RESULT.COTS_all_densities(n,t+1,META.COTS_adult_min_age:end),3); % for CoTS control
            RESULT.COTS_fecundity(n,t+1) = sum(NEW_COTS_densities.*META.COTS_fecundity).*fertilization_success ;
            RECORD.coral_pct2D_lost_COTS(n,t,:) = 100*total_coral_loss_COTS/sum(REEF(n).substrate_SA_cm2);  %now per species
            
        end
        
        %%%% --------------------------------------------------------------
        %%%% HURRICANES
        %%%% --------------------------------------------------------------
        if META.doing_hurricanes == 1 && season == 0 % only in summer
            
            RECORD.hurricane_events(n,t) = REEF(n).hurricane_chronology(1,t);
            
            if RECORD.hurricane_events(n,t)>0
                
%                 seed = rng;

                % Estimate hurricane impact with the specified category (includes category 0 = no hurricane -> no effect
                [metapop(n).coral, metapop(n).algal, total_coral_loss_cyclones] = f_hurricane_effect...
                    (metapop(n).coral, metapop(n).algal, RECORD.hurricane_events(n,t), RECORD.hurricane_parm_k(n,t), META, CORAL) ;
                
                RECORD.coral_pct2D_lost_cyclones(n,t,:) = 100*total_coral_loss_cyclones/sum(REEF(n).substrate_SA_cm2); %now per species
% rng(seed)
            end          
        end
        
        %%%% --------------------------------------------------------------
        %%%% BLEACHING
        %%%% --------------------------------------------------------------
        if META.doing_bleaching == 1 && season == 0
            
            if REEF(n).predicted_DHWs(t)>= META.DHW_threshold && RECORD.hurricane_events(n,t) == 0 % only bleach if no hurricane
                            
                [metapop(n).coral, metapop(n).genes, metapop(n).algal, total_coral_loss_bleaching, total_mortality_bleaching] = ...
                    f_bleaching_new(metapop(n).coral, metapop(n).genes, metapop(n).algal, bleaching_mortalities(n,t),...
                    CORAL, META.doing_3D, META.nb_coral_types, META.doing_clades, META.doing_genetics, ...
                    META.bleaching_whole_offset, META.bleaching_partial_offset,REEF(n).Topt_baseline, META.Topt2index);
                
                
                RECORD.applied_DHWs(n,t) = REEF(n).predicted_DHWs(1,t);
                RECORD.applied_bleaching_mortality(n,t) = bleaching_mortalities(n,t); % record the bleaching mortality effectively applied
                RECORD.coral_pct2D_lost_bleaching(n,t,:) = 100*total_coral_loss_bleaching/sum(REEF(n).substrate_SA_cm2);
            end
        end

        %%%% --------------------------------------------------------------
        %%%% RUBBLE NATURAL STABILIZATION + FORMATION
        %%%% --------------------------------------------------------------
        % Track rubble formation and natural cementation
        if META.tracking_rubble == 1
            
            RESULT.rubble_cover_pct2D(n,t+1) = RESULT.rubble_cover_pct2D(n,t)*(1-META.rubble_decay_rate) ;
            
            if RECORD.hurricane_events(n,t)>0
                RESULT.rubble_cover_pct2D(n,t+1) = RESULT.rubble_cover_pct2D(n,t+1) + ...
                    META.convert_rubble * sum(RECORD.coral_pct2D_lost_cyclones(n,t,:),3);
            end
            
            if t>META.convert_rubble_lag && sum(RECORD.coral_pct2D_lost_bleaching(n,t-META.convert_rubble_lag,:),3) >0 % Assuming structural collapse after 3 years 
                RESULT.rubble_cover_pct2D(n,t+1) = RESULT.rubble_cover_pct2D(n,t+1) + ...
                    META.convert_rubble * sum(RECORD.coral_pct2D_lost_bleaching(n,t-META.convert_rubble_lag,:),3);
            end 
                                    
            if t>META.convert_rubble_lag && sum(RECORD.coral_pct2D_lost_COTS(n,t-META.convert_rubble_lag,:),3) >0 % Assuming structural collapse after 3 years 
                RESULT.rubble_cover_pct2D(n,t+1) = RESULT.rubble_cover_pct2D(n,t+1) + ...
                    META.convert_rubble * sum(RECORD.coral_pct2D_lost_COTS(n,t-META.convert_rubble_lag,:),3);
            end  
        end
        
        %%%% --------------------------------------------------------------
        %%%% CORAL RECRUITMENT
        %%%% --------------------------------------------------------------        
        % Recruitment of 6 month old corals, ie, the settlers that passed
        % the bottleneck of posettlement mortality through the whole summer time step       
        
        if season == 0 % Only recruit in summer
                      
            if META.recruitment_type == 1
                % Max density of settlers is density dependent and given by Beverton-Holt function of available coral larvae per unit reef area (400m2)
%                 max_density_settlers = REEF(n).CORAL_recruit_survival .* CORAL.BH_alpha.*(squeeze(RESULT.coral_larval_supply(n,t,:))'./(CORAL.BH_beta + squeeze(RESULT.coral_larval_supply(n,t,:))'));    
                max_density_settlers = REEF(n).CORAL_recruit_survival' .* CORAL.BH_alpha.*(squeeze(RESULT.coral_larval_supply(n,t,:))./(CORAL.BH_beta + squeeze(RESULT.coral_larval_supply(n,t,:))));               
            
            else
                max_density_settlers = REEF(n).CORAL_recruit_survival' .* CORAL.BH_alpha;
            end
            
            
            % Apply recruitment
            [metapop(n).coral, metapop(n).genes, total_settled, ID_colony_tracking(n,:), metapop(n).algal] = f_apply_recruitment(metapop(n).coral, metapop(n).algal,...
                metapop(n).genes, META, REEF(n), max_density_settlers, CORAL.clade_prop, ID_colony_tracking(n,:));

            RESULT.coral_settler_count(n,t+1,:) = total_settled ; % record the number of settlers before they are processed
            
        else
            % No recruitment in winter
            RESULT.coral_settler_count(n,t+1,:) = 0;
        end
        
        %%%% --------------------------------------------------------------
        %%%% RE-ARRANGE CORAL MATRICES FOR OPTIMIZATION
        %%%% --------------------------------------------------------------
        if arrange_step(t)==1 % OPTIMIZATION only occurs at the selected time steps)
            [metapop(n).coral, metapop(n).genes] = f_struct_arrange(metapop(n).coral, metapop(n).genes, META);
        end
        
        %%%% --------------------------------------------------------------
        %%%% REEF RESTORATION
        %%%% --------------------------------------------------------------
        if META.doing_restoration == 1
            
            %% Split the reef into restored/non-restored sites the first time it is restored (Nov 2021)
            if RECORD.outplanted_reefs(n,t)==1 || RECORD.enriched_reefs(n,t)==1
                
                if never_deployed(n,1)==1 % if deploy for the first time on this reef, need to define cells for deployment
                    
                    coral_pct2D_cells = zeros(nb_cells, META.nb_coral_types);
                    
                    for s=1:META.nb_coral_types
                        
                        living_planar_coral_cover_cm2_cells = sum(metapop(n).coral(s).cover_cm2,2) ;
                        coral_pct2D_cells(:,s) = 100*living_planar_coral_cover_cm2_cells./REEF(n).substrate_SA_cm2  ;
                    end
                    
                    total_coral_pct2D_cells = sum(coral_pct2D_cells(REEF(n).grazable_cell==1,:),2);
                    my_cells = list_cell(REEF(n).grazable_cell==1);
                    
                    [~,J] = sort(total_coral_pct2D_cells,'ascend');
                    my_sorted_cells = my_cells(J);
                    nb_cells_to_treat = META.coral_deployment.NumberCellsTreated(n);
                    
                    if nb_cells_to_treat > length(my_sorted_cells) % if there is not enough cells without sand to match the expected deployment
                        nb_cells_to_treat = length(my_sorted_cells); % then fix the nb of cells for deployment to the number of grazable cells
                    end
                    
                    REEF(n).restored_cells = my_sorted_cells(1:nb_cells_to_treat);
                    
                    never_deployed(n,1)=0; % so that restored cells are now fixed for the rest of the simulation
                end
            end
      
            %Need to re-arrange coral matrix in case if wasn't previously done (depends on t) 
            [metapop(n).coral, metapop(n).genes] = f_struct_arrange(metapop(n).coral, metapop(n).genes, META);

            %% CORAL OUTPLANTING
            if RECORD.outplanted_reefs(n,t)==1
                
                if META.outplant_density_variable==1
                    % Deployed density was calculated at t-1 depending on total coral cover but does not necessarily match the number of
                    % corals to deploy (depends on what was available on the last deployed reef). So we need to recalculate density:
                    Density_to_outplant = META.outplant_species_prop .* All_reef_nb_corals_to_outplant(n)/META.coral_deployment.DeploymentArea_km2(n)/1e6;
                else
                    Density_to_outplant = META.outplant_species_prop .* META.outplanted_density;
                end
                
                seed=rng;
                
                [metapop(n).coral, metapop(n).algal, metapop(n).genes, ID_colony_tracking(n,:), RECORD.total_outplanted(n,t,:)] = ...
                    f_coral_deployment(metapop(n).coral, metapop(n).algal, metapop(n).genes, REEF(n), ID_colony_tracking(n,:),...
                    META, Density_to_outplant, META.outplant_coral_diameter_mean, META.outplant_coral_diameter_sd, REEF(n).restored_cells) ;
                
                rng(seed);
            end
            
            %% LARVAL ENRICHMENT (March 2022)
            if RECORD.enriched_reefs(n,t)==1
                
                if META.enrichment_density_variable==1
                    % Deployed density was calculated at t-1 depending on total coral cover but does not necessarily match the number of
                    % corals to deploy (depends on what was available on the last deployed reef). So we need to recalculate density:
                    Density_to_enrich = META.enriched_species_prop .* All_reef_nb_corals_to_enrich(n)/META.coral_deployment.DeploymentArea_km2(n)/1e6;
                else
                    Density_to_enrich = META.enriched_species_prop .* META.enriched_density;
                end
                                
                seed=rng;
                
                [metapop(n).coral, metapop(n).algal, metapop(n).genes, ID_colony_tracking(n,:), RECORD.total_enriched(n,t,:)] = ...
                    f_coral_deployment(metapop(n).coral, metapop(n).algal, metapop(n).genes, REEF(n), ID_colony_tracking(n,:),...
                    META, Density_to_enrich, META.enriched_coral_diameter_mean, META.enriched_coral_diameter_sd, REEF(n).restored_cells) ;
                
                rng(seed);
            end
            
            %% RUBBLE STABILISATION
            if RECORD.stabilised_reefs(n,t)==1 % decision to undertake rubble stab was taken at t-1 (see below)
                
                RECORD.rubble_cover_pct2D_stabilised(n,t+1) = META.proportion_rubble_stabilised*RESULT.rubble_cover_pct2D(n,t+1);
                RESULT.rubble_cover_pct2D(n,t+1) = RESULT.rubble_cover_pct2D(n,t+1) - META.proportion_rubble_stabilised*RESULT.rubble_cover_pct2D(n,t+1) ;
            end
            
            %Record coral cover on the restored cells
            if never_deployed(n,1)==0 % only after a first deployment has already occurred
                
                for s=1:META.nb_coral_types
                    
                    living_planar_coral_cover_cm2 = sum(metapop(n).coral(s).cover_cm2(REEF(n).restored_cells,:),2) ;
                    RESULT.coral_pct2D_restored_sites(n,t+1,s) = 100*sum(living_planar_coral_cover_cm2)./sum(REEF(n).substrate_SA_cm2(REEF(n).restored_cells))  ;
                end
            end
        end
        
        %%%% --------------------------------------------------------------
        %% %% FINAL ESTIMATES FOR REEF n
        %%%% --------------------------------------------------------------
        %%%% Estimate coral fecundity of the reef
        %%%% --------------------------------------------------------------
        for s=1:META.nb_coral_types % this is to be used for the next time step, so t+1
            
            find_colonies =  find(metapop(n).coral(s).cover_cm2>0);
            
            if isempty(find_colonies)==0
                
                if META.doing_genetics == 1
                    
                    list_coral_ID = metapop(n).genes(s).list_coral_ID ;
                    list_new = metapop(n).coral(s).colony_ID(find_colonies) ;
                    
                    check = ismember(list_coral_ID,list_new) ;
                    metapop(n).genes(s).QTLs(check==0,:,:)=[];
                    metapop(n).genes(s).list_coral_ID(check==0,:)=[];
                    metapop(n).genes(s).phenotypes(check==0,:)=[];
                    
                    % Difference between phenotype (thermal optimum) and environment
                    D = REEF(n).predicted_SST(t) - REEF(n).Topt_baseline - metapop(n).genes(s).phenotypes ;
                    
                    % Estimate relative fitness with the breadth of thermal tolerance (SIGMA)
                    F_list = ones(size(D));
                    F_list(D>0) = normpdf(D(D>0),0,META.genetics.SIGMA_HOT(s))/normpdf(0,0,META.genetics.SIGMA_HOT(s));
                    F_list(D<0) = normpdf(D(D<0),0,META.genetics.SIGMA_COLD(s))/normpdf(0,0,META.genetics.SIGMA_COLD(s));
                     
                    RESULT.relative_fitness(n,t+1,s) = mean(F_list);
                    
                    All_Topts = round(10*(REEF(n).Topt_baseline + metapop(n).genes(s).phenotypes))/10 ;
                    
                    [count_Topt, ~] = hist(All_Topts,Topt_bins) ;
                    RESULT.Topt_distri(n,t+1,s,:) = count_Topt;
                    RESULT.Topt_mean(n,t+1,s) = mean(All_Topts);
                    
                    % Then create the pool of genotypes to be spread
                    
                    % Genotype of parents
                    QTLs = metapop(n).genes(s).QTLs ;
                    
                    select1 = randsample(1:length(list_new), META.genetics_pop_size, 'true', F_list')'; % pick up parent 1 for each larva
                    select2 = randsample(1:length(list_new), META.genetics_pop_size, 'true', F_list')'; % pick up parent 2 for each larva (note this does not avoid selfing)

                    % Look for hermaphrodites
                    test_self =  select1 - select2;
                    select2(test_self==0) = randi(length(list_new), length(test_self(test_self==0)),1);
                    % Note this won't completely exclude self-fertilization, though it should be much reduced'
                    % except when local population is made of only few individuals (in which case the number of larvae
                    % should be relatively low so that the impact of self-fertiization on the network shoud be negligible
                    
                    REEF(n).coral(s).QTL_pool_OUT(1:META.genetics_pop_size,:,1) = squeeze(QTLs(select1,:,1));
                    REEF(n).coral(s).QTL_pool_OUT(1:META.genetics_pop_size,:,2) = squeeze(QTLs(select2,:,1));
                    
                    % Record the mutation rate on each source reef
                    % Mutation rate is scaled to the estimated population size, which is nb colonies * habitat area
                    LocalMuteRate = REEF(n).scale_mutation_rate*META.genetics.MuteRate(s)*ID_colony_tracking(n,s)/META.genetics_pop_size;
                    rand_MuteIncidence = rand(META.genetics_pop_size, META.genetics.nb_loci, 2)<LocalMuteRate; % rand_MuteIncidence is logical!
                    mutations = zeros(size(rand_MuteIncidence));
                    
                    mutations(rand_MuteIncidence==1) = normrnd(0, META.genetics.MuteEffect(s),sum(sum(sum(rand_MuteIncidence))),1);
                    REEF(n).coral(s).QTL_pool_OUT = REEF(n).coral(s).QTL_pool_OUT + mutations ;
                    
                    % Update the pool of genotypes coming in = due to retention
                    % Note this pool will be updated at the next step with larvae coming from other reefs
                    REEF(n).coral(s).QTL_pool_IN = REEF(n).coral(s).QTL_pool_OUT ;
                    
                else
                    % (need F_list to run fecundity any way)
                    F_list = ones(length(metapop(n).coral(s).cover_cm2(find_colonies)),1); % everyone is fit
                    
                end
                
                % Then estimate the number of larvae produced by the reef (per 400m2 of reef)
                % = total number of larvae produced per species for each reef at every time step
                [RESULT.coral_total_fecundity(n,t+1,s)] = f_estimate_fecundity(metapop(n).coral(s).cover_cm2(find_colonies), F_list,...
                CORAL.fecund_min_size(s), CORAL.fecund_a(s), CORAL.fecund_b(s));
                
            end
        end
        
        %%%% --------------------------------------------------------------
        %%%% Estimate 3D areas and RUGOSITY
        %%%% --------------------------------------------------------------
        if META.doing_3D == 1
            do_erosion=1;
            [metapop(n).coral,volume_dead_eroded_cm3] = f_estimate_3D_colony(metapop(n).coral, do_erosion, REEF(n).fish_bioerosion, CORAL, META, last_surface_area_grazed);
            [metapop(n).coral, metapop(n).algal, REEF(n)]= f_estimate_3D_reef(metapop(n).coral, metapop(n).algal, META.cell_area_cm2, REEF(n),META);
            
            RESULT.volume_eroded(n,t+1) = sum(volume_dead_eroded_cm3) ;
            RESULT.substrate_cm2(n,t+1) = sum(REEF(n).substrate_SA_cm2);
        end
        
        %%%% --------------------------------------------------------------
        %%%% Report current reef status
        %%%% --------------------------------------------------------------
        
        % 1) Algal cover: total cover (%) of each algal type over the grid at time step t
        RESULT.algal_pct(n,t+1,:) = 100*full(sum([metapop(n).algal(:).cover_cm2]))./sum(REEF(n).substrate_SA_cm2)  ;
        
       
        % 2) Coral cover: size structure of coral populations at t
        for s=1:META.nb_coral_types
            
            living_planar_coral_cover_cm2 = metapop(n).coral(s).cover_cm2(metapop(n).coral(s).cover_cm2>0) ;
            RESULT.coral_pct2D(n,t+1,s) = 100*sum(living_planar_coral_cover_cm2)./sum(REEF(n).substrate_SA_cm2)  ;
            
            if META.doing_clades == 1
                living_coral_clades = metapop(n).coral(s).clade(metapop(n).coral(s).cover_cm2>0) ;
                RESULT.clade_prop(n,t+1,s) = nnz(living_coral_clades==1)/nnz(living_coral_clades);
            end
            
            % Estimaton of living and dead carbonate volumes (per species)
            if META.doing_3D == 1
                living_coral_volume_cm3 = metapop(n).coral(s).volume_cm3(metapop(n).coral(s).cover_cm2>0) ;
                dead_coral_volume_cm3 = metapop(n).coral(s).volume_cm3(metapop(n).coral(s).cover_cm2<0) ;
                RESULT.live_CaCO3(n,t+1,s)= sum(living_coral_volume_cm3);
                RESULT.dead_CaCO3(n,t+1,s)= sum(dead_coral_volume_cm3);
            else % if not doing 3D we remove the dead colonies to speed up the code
                metapop(n).coral(s).cover_cm2(metapop(n).coral(s).cover_cm2<0)=0;
            end
            
            % Coral size structure at time step t
            if META.doing_size_frequency == 1
                [count, ~] = f_count_sizefreq(living_planar_coral_cover_cm2, CORAL) ;
                RESULT.coral_juv_count(n,t+1,s,:)=count.juv;
                RESULT.coral_adol_count(n,t+1,s,:)=count.adol;
                RESULT.coral_adult_count(n,t+1,s,:)=count.adult;
            end
            
        end
        
        % 3) Rugosity at time step t
        if META.doing_3D == 1
            coral_sp_pct2D = squeeze(squeeze(RESULT.coral_pct2D(n,t+1,:)));
            [rugosity, SI_reef] = f_estimate_rugosity(REEF(n), META, CORAL, coral_sp_pct2D);
            RESULT.rugosity(n,t+1)=rugosity;
            RESULT.SI_reef(n,t+1)=SI_reef;
        end
        
    end % end of the loop over the n reefs
%     rng(seed)
    %%%% --------------------------------------------------------------
    %%%% REGISTER RESTORATION NEEDS FOR THE NEXT STEP
    %%%% --------------------------------------------------------------
    if META.doing_restoration == 1
        
        % 1) Coral outplanting -----------------------------------------  
        if META.nb_reefs_outplanted > 0
            
            if t < META.nb_time_steps && META.doing_coral_outplanting(t+1)==1
                
                if META.outplant_density_variable==1
                    
                    All_reef_states = squeeze(sum(RESULT.coral_pct2D(:,t+1,:),3));
                    All_reef_outplanting_densities = META.outplant_density_aquaculture.intercept + ...
                        META.outplant_density_aquaculture.slope * All_reef_states;
                    
                    Priority_reef_ID = META.reef_ID(META.priority_list_Outplant);
                    Priority_reef_states = All_reef_states(META.priority_list_Outplant);
                    select_reefs = find(Priority_reef_states < META.threshold_for_outplanting.max & Priority_reef_states>META.threshold_for_outplanting.min);
                    Priority_reef_ID = Priority_reef_ID(select_reefs);
                    % Exclude from priority list reefs with large amount of corals (no need to outplant on those)
                    %  Priority_reefs(All_reef_states(META.priority_list_Outplant)>META.threshold_for_outplanting.max)=[];
                    %  Priority_reefs(All_reef_states(META.priority_list_Outplant)<META.threshold_for_outplanting.min)=[];
                    % Note here could narrow down to META.nb_reefs_outplanted if this is a seleted option
                    
                    r = 0;
                    Total_available_outplants = META.total_nb_outplants;
                    All_reef_nb_corals_to_outplant = zeros(size(All_reef_states)); % total nb of outplants that should be deployed
                    
                    while Total_available_outplants > 0 && r+1 <= length(Priority_reef_ID) % this stops when no more priorty reefs
                        r = r + 1;
                        reef_index = find(META.reef_ID == Priority_reef_ID(r));
                        
                        All_reef_nb_corals_to_outplant(reef_index) = round(All_reef_outplanting_densities(reef_index)*Total_available_outplants...
                            *META.coral_deployment.DeploymentArea_km2(reef_index));
                        Total_available_outplants = Total_available_outplants - All_reef_nb_corals_to_outplant(reef_index);
                        RECORD.outplanted_reefs(reef_index,t+1)=1;
                    end
                    
                    if sum(All_reef_nb_corals_to_outplant) > META.total_nb_outplants
                        All_reef_nb_corals_to_outplant(reef_index) = All_reef_nb_corals_to_outplant(reef_index)+...
                            META.total_nb_outplants-sum(All_reef_nb_corals_to_outplant);
                    end
                    
                else % fixed density of outplants
                    
                    RECORD.outplanted_reefs(META.priority_list_Outplant,t+1)=1;

                end
            end
        end

       
        % 2) Larval enrichment -----------------------------------------  
        if META.nb_reefs_enriched > 0
            
            if t < META.nb_time_steps && META.doing_larval_enrichment(t+1)==1
                
                if META.enrichment_density_variable==1 % if 1, deployed density depends on how much coral cover is already there
                    
                    All_reef_states = squeeze(sum(RESULT.coral_pct2D(:,t+1,:),3));
                    
                    Priority_reef_ID = META.reef_ID(META.priority_list_LarvalEnrich);
                    Priority_reef_states = All_reef_states(META.priority_list_LarvalEnrich);
                    select_reefs = find(Priority_reef_states < META.threshold_for_larval_enrichment);
                    Priority_reef_ID = Priority_reef_ID(select_reefs);
                    
                    r = 0;
                    Total_available_larvae = META.total_nb_larvae;
                    All_reef_nb_corals_to_enrich = zeros(size(All_reef_states));

                    while Total_available_larvae > 0 && r+1 <= length(Priority_reef_ID) % this stops when no more priorty reefs
                        r = r + 1;
                        reef_index = find(META.reef_ID == Priority_reef_ID(r));
                        
                        % Note here that the enriched density is the same for all reefs
                        All_reef_nb_corals_to_enrich(reef_index) = round(META.enriched_density*Total_available_larvae...
                            *META.coral_deployment.DeploymentArea_km2(reef_index));
                        Total_available_larvae = Total_available_larvae - All_reef_nb_corals_to_enrich(reef_index);
                        RECORD.enriched_reefs(reef_index,t+1)=1;
                    end
                    
                    if sum(All_reef_nb_corals_to_enrich) > META.total_nb_larvae
                        All_reef_nb_corals_to_enrich(reef_index) = All_reef_nb_corals_to_enrich(reef_index)+...
                            META.total_nb_larvae-sum(All_reef_nb_corals_to_enrich);
                    end
                    
                else % fixed density of larvae
                    
                    RECORD.enriched_reefs(META.priority_list_LarvalEnrich,t+1)=1;

                end
            end
        end

        
        % 3) Fogging stabilisation -------------------------------------
        if META.nb_reefs_fogged > 0
            
            if t < META.nb_time_steps && iseven(t+1)==0 && META.doing_fogging(t+1)==1

                RECORD.fogged_reefs(META.priority_list_Fogging,t+1)=1;

                % fogging just reduces the bleaching mortality predicted by
                % DHW (see initialisation above)
                bleaching_mortalities(META.priority_list_Fogging,t+1) = ...
                    bleaching_mortalities(META.priority_list_Fogging,t+1)*META.bleaching_mortality_under_fogging;     
            end
        end

       
        % 3) Rubble stabilisation -------------------------------------
        if META.nb_reefs_stabilised > 0
            
            if t < META.nb_time_steps && META.doing_rubble_stabilisation(t+1) == 1
                
                All_reef_rubble = RESULT.rubble_cover_pct2D(:,t+1);
                
                Priority_reef_ID = META.priority_list_RubbleStab; % default priority list defined in settings_RESTORATION.m
                
                % Exclude from priority list reefs with small amount of rubble on reefs
                Priority_reef_ID(All_reef_rubble(META.priority_list_RubbleStab) < META.threshold_for_stabilisation)=[];
                
                if isempty(Priority_reef_ID)==0
                    
                    nb_restored_reefs = META.nb_reefs_stabilised ;
                    nb_restored_reefs(nb_restored_reefs>length(Priority_reef_ID)) = length(Priority_reef_ID);
                    
                    % Give a green light for for rubble stabilisation at the next step (if timing enables it) 
                    RECORD.stabilised_reefs(Priority_reef_ID(1:nb_restored_reefs),t+1)=1;
                                        
                end
            end
        end
        
    end
    
    %%%% --------------------------------------------------------------
    %%%% COTS control (Karlo & Caro's code 2020-21)
    %%%% --------------------------------------------------------------
    %%New code with COTS_control_CSIRO function
    if META.doing_COTS == 1 && META.doing_COTS_control == 1
        if t >= META.COTS_control_start
            
            if META.do_CSIRO_COTSctrl == 1
                [ RESULT ] = f_COTS_control_CSIRO( META, RESULT, t, REEF_POP, CONNECT);
            else
                [ RESULT ] = f_COTS_control_K( META, RESULT, t, REEF_POP, CONNECT);
            end
            
            if META.report_COTS_control == 0 && isfield(RESULT,{'COTS_records'}) == 1 % if not interested by the detailed results of culling
                RESULT=rmfield(RESULT,'COTS_records'); % then delete them
            end
        end
    end
    
    %%%% --------------------------------------------------------------
    %%%% Record individual colonies
    %%%% --------------------------------------------------------------
    if META.track_populations == 1       
        [environ_list,colony_list] = f_track_populations(environ_list,colony_list,metapop(n).algal,metapop(n).coral,t,META,REEF(n).grazable_cell);
    end
        
end % end of the time steps loop for a single simulation


if META.track_populations == 1
    f_generate_track_files(META, REEF, CORAL, ALGAL, RECORD, colony_list, environ_list)
end
