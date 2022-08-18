%__________________________________________________________________________
%
% REEFMOD-GBR - populate reef parameters and initial state 
%
% Yves-Marie Bozec, y.bozec@uq.edu.au, 09/2019
%__________________________________________________________________________


%% Populate REEF parameters
for n = 1:length(META.reef_ID)
        
    REEF(n).initial_coral_cover = init_coral_cover(n,:);
    REEF(n).initial_algal_cover = [ 0 ; 0.01 ; 0.01 ; 0 ] ;
    REEF(n).nongrazable_substratum = init_sand_cover(n,1) ;
    REEF(n).initial_rubble_pct = 100*init_rubble_cover(n,1) ;
    
    % Default values (allows for reef-specific modifications in runmodel)
    REEF(n).juv_whole_mortality_rate = CORAL.juv_whole_mortality_rate;
    REEF(n).adol_whole_mortality_rate = CORAL.adol_whole_mortality_rate;
    REEF(n).adult_whole_mortality_rate = CORAL.adult_whole_mortality_rate;
    
    % Default values for herbivory
    REEF(n).diadema = REEF(1).diadema ;
    REEF(n).herbivory = REEF(1).herbivory ;
    REEF(n).dictyota_declines_seasonally = REEF(1).dictyota_declines_seasonally ;
    
    % Store the bleaching scenario
    if META.doing_genetics == 0
        
        REEF(n).predicted_DHWs = squeeze(DHW(n,1:META.nb_time_steps))+ RESTORATION.doing_cooling*RESTORATION.cooling_factor*12;
        REEF(n).Topt_baseline = [];
        REEF(n).predicted_SST = zeros(1,META.nb_time_steps);
        
    else
        REEF(n).predicted_DHWs = squeeze(DHW(n,1:META.nb_time_steps,:))'+ RESTORATION.doing_cooling*RESTORATION.cooling_factor*12; % with all the Topt stacks
        REEF(n).Topt_baseline = MMM_CoRTAD(META.reef_ID(n),1) - META.Topt_offset ;
        REEF(n).predicted_SST = SST_GBR(n,:) + RESTORATION.doing_cooling*RESTORATION.cooling_factor/2;  % Decreases SST for only 3 months in summer

        for s=1:META.nb_coral_types
            
            % select genotypes from the bank and add random random variation (creates genetic diversity among species)
            REEF(n).coral(s).QTL_pool_IN = squeeze(QTL1000(META.reef_ID(n),1:META.genetics_pop_size,:,:)) + ...
                    normrnd(META.initial_push_QTL_mu(s), META.initial_push_QTL_sd(s),META.genetics_pop_size, META.genetics.nb_loci,2);
                
        end
        
    end
       
    % Store the scenario of cyclones
    REEF(n).hurricane_chronology = CYCLONE_CAT(n,1:META.nb_time_steps);
    
    % Assign to every reef the defaut background density of CoTS
    if META.doing_COTS == 1
        REEF(n).COTS_background_density = META.COTS_background_density;
    end
    
end
