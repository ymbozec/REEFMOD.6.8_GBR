%__________________________________________________________________________
%
% REEFMOD-GBR MAIN SCRIPT
%
% This is the clean version of 6.7 used for RRAP counterfactuals (Feb 2022) and interventions (Apr 2022).
% Hindcast has been extended to 2022, assuming no severe bleaching & cyclones in 2021 and 2022.
% Uses the latest available geomorphic maps (3D) for coral habitats (4 geomorphic classes).
% Enables the simulation of a specific region (with a proxy of larval input from excluded reefs).
% Allows tracking corals on specific reef sites for restoration interventions.
% The random generation of coral mortality due to bleaching and cyclones is now specific to run ID
%
% Yves-Marie Bozec, y.bozec@uq.edu.au, 08/2022
%__________________________________________________________________________
clear

addpath(genpath('.'));
SaveDir ='';

 
% NB_TIME_STEPS has to be an even number <= 26+158 (length of projected DHW time series is 79 years)
% We always run the hindcast (26 time steps) before future projections
% Initialisation = winter 2007

% NB_TIME_STEPS = 30; % HINDCAST: summer 2008 to winter 2022
NB_TIME_STEPS = 26+60; % HINDCAST+FORECAST summer 2008 - winter 2050
% NB_TIME_STEPS = 26+100; % HINDCAST+FORECAST summer 2008 - winter 2070
% NB_TIME_STEPS = 26+158; % HINDCAST+FORECAST summer 2008 - winter 2099

% OutputName = 'R0_FORECAST_GBR'; options = [1 1 1 1 0]; % see options below
OutputName = 'R0_BCA_Moore'; options = [1 1 1 1 1]; % Counterfactual of the Business Case Assessment (35x35 grid cells)
% Requires turning ON restoration, outplanted_density = 0 (ie, deployment of zero coral) and the year of first (ghost) 
% 'deployment' to keep track of intervention sites within reefs (ie, cells with deployment of zero coral) that are
% the exact same as in the intervention scenarios
% OutputName = 'R1_BCA_Moore'; options = [1 1 1 1 1]; % see options below

% select the Global Circulation Model for climate change projection 
GCM = 3; % 1=CCSM4, 2=CESM1_WACCM, 3=MIROC5, 4=GFDL_ESM2M, 5=HadGEM2, 6=GISS_E2_R
% Select the carbon emission pathway
RCP = 2; % 1=RCP2.6, 2=RCP4.5, 3=RCP6.0, 4=RCP8.5

%% --------------------------------------------------------------------------------
GCM_list=["CCSM4";"CESM1_WACCM";"MIROC5";"GFDL_ESM2M";"HadGEM2";"GISS_E2_R"];
RCP_list = ["26"; "45"; "60"; "85"];
OPTIONS.GCM = GCM_list(GCM);
OPTIONS.RCP = RCP_list(RCP);

% Stressor options: yes(1)/no(0)
OPTIONS.doing_cyclones = options(1);
OPTIONS.doing_bleaching = options(2) ; 
OPTIONS.doing_COTS = options(3);
OPTIONS.doing_WQ = options(4);

OPTIONS.doing_size_frequency = 1; % for tracking population size structure (incl. juveniles)
OPTIONS.doing_adaptation = 0 ;
OPTIONS.adaptation_parms = [ 1 1 1.5 ]; % sigma cold/sigma hot/esd
OPTIONS.doing_restoration = options(5) ;

% Below options are used to force simulations with specific starting conditions (eg, building LUT for the RRAP-RE)
% Keep them empty if of no use
OPTIONS.init_coral_cover = []; %0.01*ones(1,6); % as proportional cover (vector of 6 values)
OPTIONS.init_sand_cover = []; % as proportional cover 
OPTIONS.init_rubble_cover = []; % as proportional cover 
OPTIONS.ssc = []; %0.1; %in mg/L

%% CoTS control
OPTIONS.doing_COTS_control= 1; % Note control is set to start in 2019 (after 23 time steps)
OPTIONS.CoTS_control_scenarios = csvread('parsList2.csv', 0, 1);
% Caro: 5 boats, whole GBR, start control in 2019. New list adjusted for matching Coconet
% Option names in parsList2.csv. Second value is spatial strategy (14: whole GBR under CoTS control).
% Fourth value is number of boats (5)
OPTIONS.CoTS_control_scenarios(4)=5;

%% Outplanting
% Set the restoration effort: maximimum number of reefs where coral outplanting is undertaken at random at each time step
RESTORATION.nb_reefs_outplanted = Inf ;
% Set to Inf if unlimited OR if specific reefs are restored (listed in SETTINGS_RESTORATION
% If 0, outplanting cannot happen. For the counterfactual, set to Inf with outplanted_density = 0 for ghost deployment

RESTORATION.total_nb_outplants = Inf; % Max number of outplants available for all reefs at each time step.
% Set to Inf if outplant density is the driver

RESTORATION.outplanted_density = 0;  %only in the case of fixed density of outplants (ignores RESTORATION.total_nb_outplants)

RESTORATION.doing_coral_outplanting = zeros(1,NB_TIME_STEPS); % if zero, no coral deployment at time step for all reefs
RESTORATION.doing_coral_outplanting(1,38:2:46) = 1 ; % set to 1 to indicate outplanting: first in summer 2026 and last in summer 2030 (BCA)
% RESTORATION.doing_coral_outplanting(1,[38:2:46 78:2:86 118:2:126 158:2:166] ) = 1 ; % outplanting starts 2026-2030, 2046-2050, 2066-2070, 2086-2090

RESTORATION.thermal_tolerance_outplants = 0 ; % DegC above thermal tolerance of native corals (WITH GENETICS)
RESTORATION.DHW_tolerance_outplants = 4 ; % DHW tolerance relative to native corals (WITHOUT GENETICS)
RESTORATION.tradeoff_growth_outplant = 1; %0.6; % 40% reduction of growth rate relative to native corals

RESTORATION.outplant_species = [1 2 3 4 5 6]; % ID of outplanted coral types - this will create as many 'new' types after the 6 default types
RESTORATION.outplant_species_prop = [0.02 0.14 0.14 0 0.7 0]; % Taxonomic composition of outplanted corals: must sum to 1
RESTORATION.outplant_diameter_mean = [2.56 2.56 2.56 1.41 1.41 1.41] ; % mean diameter (in cm) of outplants of outplants of each deployed type
RESTORATION.outplant_diameter_sd = [0.26 0.26 0.26 0.14 0.14 0.14] ; % sd diameter (in cm) of outplants of outplants of each deployed type

%% Rubble stabilisation
% Set the restoration effort: number of reefs where rubble is stabilised at each time step
RESTORATION.nb_reefs_stabilised = 0 ; % (if 0 rubble stabilisation cannot happen) 
% Set the timing of intervention (if 1 intervention is deployed at step t, if 0 no intervention at t)
RESTORATION.doing_rubble_stabilisation = zeros(1,NB_TIME_STEPS);
% RESTORATION.doing_rubble_stabilisation(1,38:1:end) = 1 ; % set to 1 to do outplanting: first in summer 2026 and last in summer 2030

%% Larval enrichment
RESTORATION.nb_reefs_enriched = Inf ; % max number of reefs where larval enrichment is undertaken at each time step
% Set to Inf if unlimited OR if specific reefs are restored (listed in SETTINGS_RESTORATION)
% If 0, enrichment cannot happen. For the counterfactual, set to Inf with total_nb_larvae = 0 for ghost deployment

RESTORATION.total_nb_larvae = 1e6; % Max number of 'larvae' (ie, 1 yr old corals) available at each time step. Set to Inf if unlimited.

RESTORATION.doing_larval_enrichment = zeros(1,NB_TIME_STEPS);
% RESTORATION.doing_larval_enrichment(1,38:2:46) = 1 ; % set to 1 to do outplanting

RESTORATION.DHW_tolerance_larvae = 0;

%% SRM (cooling)
% Cloud brightening from RRAP feasibility study (needs revision)
RESTORATION.doing_cooling = 0 ;
RESTORATION.cooling_factor = 0 ; % otherwise [-0.3 ; -0.7 ; -1.3];

% Fogging for 2022 intervention simulations
RESTORATION.nb_reefs_fogged = Inf;
% Select reefs for the deployment of fogging: 695-Moore; 697-Elford; 698-Briggs ; 969-Milln; 970-Thetford
RESTORATION.fogged_reef_ID = [969 970]; % 1 fogging unit (4.5 km2) covering 4.6 km2 (equivalent 2D reef areas)
% RESTORATION.fogged_reef_ID = 695; % 2 fogging units (9 km2) covering 8.7 km2 (equivalent 2D reef areas)
% RESTORATION.fogged_reef_ID = [695 697 969 970]; % 5 fogging units (22.5 km2) covering 23 km2 (equivalent 2D reef areas)

RESTORATION.doing_fogging = zeros(1,NB_TIME_STEPS); %if zero don't do fogging at time step
% RESTORATION.doing_fogging(1,1:2:end) = 1 ; % set to 1 to do fogging - only in summer!!
RESTORATION.bleaching_mortality_under_fogging = 0.8 ; % fogging reduces bleaching mortality by 20%

%% Saving options
if NB_TIME_STEPS > 30   
    OPTIONS.OutputFileName = [SaveDir OutputName '_' char(OPTIONS.GCM) '_' char(OPTIONS.RCP) '.mat'];
else
    OPTIONS.OutputFileName = [SaveDir OutputName '.mat'];
end

%% --------------------------------------------------------------------------------
OUTPUTS = struct('REEF', [],'RESULT', [],'RECORD', []);
TEMP_META = struct('META', []);

% parfor run_id = 1:NB_SIMULATIONS
for run_id = 1:NB_SIMULATIONS
    
    run_id
    
    [meta, REEF, RESULT, RECORD] = f_multiple_reef(OPTIONS, RESTORATION, NB_TIME_STEPS, run_id);
      
    OUTPUTS(run_id).RESULT = RESULT ;
    OUTPUTS(run_id).RECORD = RECORD ;
    OUTPUTS(run_id).REEF = REEF ;
    TEMP_META(run_id).META = meta ;
    
end

META = TEMP_META(1).META; % Keep only one META because common to all simulations
clear TEMP_META ADAPT run_id meta RECORD REEF RESULT GCM GCM_list options SaveDir OutputName RCP RCP_list


%% --------------------------------------------------------------------------------
%% Memory allocation
% 1) coral outputs
coral_cover_per_taxa = zeros(NB_SIMULATIONS,META.nb_reefs,META.nb_time_steps+1,META.nb_coral_types,'single');
coral_larval_supply = coral_cover_per_taxa;
nb_coral_offspring = coral_cover_per_taxa;
nb_coral_recruit = zeros(NB_SIMULATIONS,META.nb_reefs,META.nb_time_steps+1,META.nb_coral_types,'uint16');
    
if OPTIONS.doing_size_frequency == 1   
    nb_coral_juv = zeros(NB_SIMULATIONS, META.nb_reefs, META.nb_time_steps+1, META.nb_coral_types, 2, 'uint16') ;
    nb_coral_adol = zeros(NB_SIMULATIONS, META.nb_reefs, META.nb_time_steps+1, META.nb_coral_types, 3, 'uint16') ;
    nb_coral_adult = zeros(NB_SIMULATIONS, META.nb_reefs, META.nb_time_steps+1, META.nb_coral_types, 6, 'uint16') ; 
end

% 2) Coral cover loss following stressors
coral_cover_lost_bleaching = zeros(NB_SIMULATIONS, META.nb_reefs, META.nb_time_steps, META.nb_coral_types, 'single');
coral_cover_lost_cyclones = coral_cover_lost_bleaching;
coral_cover_lost_COTS = coral_cover_lost_bleaching;

% 3) Stress records
record_applied_DHWs = zeros(NB_SIMULATIONS, META.nb_reefs, META.nb_time_steps,'single');
record_applied_bleaching_mortality  = zeros(NB_SIMULATIONS, META.nb_reefs, META.nb_time_steps,'single');
record_applied_cyclones = zeros(NB_SIMULATIONS, META.nb_reefs, META.nb_time_steps,'uint16');

% 4) Restoration records
if OPTIONS.doing_restoration==1
    
    % Total nb of outplants per reef per time step
    record_total_outplants_deployed = zeros(NB_SIMULATIONS, META.nb_reefs, META.nb_time_steps, length(META.outplant_species),'uint16');
    record_outplanted_reefs = zeros(NB_SIMULATIONS, META.nb_reefs, META.nb_time_steps+1,'uint16');
    
    record_total_larvae_deployed = zeros(NB_SIMULATIONS, META.nb_reefs, META.nb_time_steps, length(META.enriched_species),'uint16');
    record_enriched_reefs = zeros(NB_SIMULATIONS, META.nb_reefs, META.nb_time_steps+1,'uint16');  
    
    record_rubble_pct2D_stabilised = zeros(NB_SIMULATIONS, META.nb_reefs, META.nb_time_steps,'single');
    record_stabilised_reefs = zeros(NB_SIMULATIONS, META.nb_reefs, META.nb_time_steps+1,'uint16');
    
    coral_cover_per_taxa_restored_sites = zeros(NB_SIMULATIONS,META.nb_reefs,META.nb_time_steps+1,META.nb_coral_types,'single');
end

% 5) Other variables
rubble = zeros(NB_SIMULATIONS, META.nb_reefs, META.nb_time_steps+1,'single');
nongrazable = zeros(NB_SIMULATIONS, META.nb_reefs,'single');
macroTurf = zeros(NB_SIMULATIONS, META.nb_reefs, META.nb_time_steps+1,'single');
macroEncrustFleshy = zeros(NB_SIMULATIONS, META.nb_reefs, META.nb_time_steps+1,'single');
macroUprightFleshy = zeros(NB_SIMULATIONS, META.nb_reefs, META.nb_time_steps+1,'single');

% 6) CoTS outputs
if OPTIONS.doing_COTS == 1
    COTS_mantatow = zeros(NB_SIMULATIONS, META.nb_reefs, META.nb_time_steps+1, 'single');
    COTS_densities0 = zeros(NB_SIMULATIONS, META.nb_reefs, META.nb_time_steps+1, 16, 'single'); % 16 age classes
    COTS_settler_density = zeros(NB_SIMULATIONS, META.nb_reefs, META.nb_time_steps+1, 'uint16');
    COTS_larval_supply = COTS_mantatow;
    COTS_larval_output = COTS_mantatow;
end

if OPTIONS.doing_COTS_control == 1 && META.nb_time_steps > 23
    COTS_CONTROL_culled_reefs = zeros(NB_SIMULATIONS, META.nb_reefs, META.nb_time_steps); % gives which reefs were culled or not (1/0)
    COTS_CONTROL_remaining_dives = zeros(NB_SIMULATIONS, META.nb_time_steps); %remaining number of control dives available after intervention
    COTS_CONTROL_culled_density = zeros(NB_SIMULATIONS, META.nb_reefs, META.nb_time_steps); % density extracted by the intervention (adults)
end

%% Populate outputs
for simul = 1:NB_SIMULATIONS
    
    coral_cover_per_taxa(simul,:,:,:) = squeeze(cat(4,OUTPUTS(simul).RESULT.coral_pct2D));
    coral_larval_supply(simul,:,:,:) = squeeze(cat(4,OUTPUTS(simul).RESULT.coral_larval_supply)); % nb of incoming larvae per unit of reef area (400m2)
    nb_coral_recruit(simul,:,:,:) = squeeze(cat(4,OUTPUTS(simul).RESULT.coral_settler_count));
    nb_coral_offspring(simul,:,:,:) = squeeze(cat(4,OUTPUTS(simul).RESULT.coral_total_fecundity)); % nb of larvae produced per unit of reef area (400m2)
    
    if OPTIONS.doing_size_frequency == 1
        nb_coral_juv(simul,:,:,:,:)= squeeze(cat(4,OUTPUTS(simul).RESULT.coral_juv_count(:,:,:,:)));
        nb_coral_adol(simul,:,:,:,:) = squeeze(cat(4,OUTPUTS(simul).RESULT.coral_adol_count(:,:,:,:)));
        nb_coral_adult(simul,:,:,:,:) = squeeze(cat(4,OUTPUTS(simul).RESULT.coral_adult_count(:,:,:,:)));
    end
    
    coral_cover_lost_bleaching(simul,:,:,:) = squeeze(OUTPUTS(simul).RECORD.coral_pct2D_lost_bleaching);
    coral_cover_lost_cyclones(simul,:,:,:) = squeeze(OUTPUTS(simul).RECORD.coral_pct2D_lost_cyclones);
    coral_cover_lost_COTS(simul,:,:,:) = squeeze(OUTPUTS(simul).RECORD.coral_pct2D_lost_COTS);
    
    record_applied_DHWs(simul,:,:) = squeeze(OUTPUTS(simul).RECORD.applied_DHWs);
    record_applied_bleaching_mortality(simul,:,:) = squeeze(OUTPUTS(simul).RECORD.applied_bleaching_mortality);
    record_applied_cyclones(simul,:,:) = squeeze(OUTPUTS(simul).RECORD.hurricane_events);
    
    if  OPTIONS.doing_restoration==1
        record_total_outplants_deployed(simul,:,:,:) = squeeze(OUTPUTS(simul).RECORD.total_outplanted);
        record_outplanted_reefs(simul,:,:) = squeeze(OUTPUTS(simul).RECORD.outplanted_reefs);
        record_total_larvae_deployed(simul,:,:,:) = squeeze(OUTPUTS(simul).RECORD.total_enriched);
        record_enriched_reefs(simul,:,:) = squeeze(OUTPUTS(simul).RECORD.enriched_reefs);
        record_rubble_pct2D_stabilised(simul,:,:) = squeeze(OUTPUTS(simul).RECORD.rubble_cover_pct2D_stabilised(:,1:end));
        record_stabilised_reefs(simul,:,:) = squeeze(OUTPUTS(simul).RECORD.stabilised_reefs);
        
        coral_cover_per_taxa_restored_sites(simul,:,:,:) = squeeze(cat(4,OUTPUTS(simul).RESULT.coral_pct2D_restored_sites));
    end
    
    rubble(simul,:,:)= squeeze(OUTPUTS(simul).RESULT.rubble_cover_pct2D);
    nongrazable(simul,:) = squeeze(cat(4,OUTPUTS(simul).REEF.nongrazable_substratum));
    
    macroTurf(simul,:,:) = squeeze(cat(4,OUTPUTS(simul).RESULT.algal_pct(:,:,4)));
    macroEncrustFleshy(simul,:,:) = squeeze(cat(4,OUTPUTS(simul).RESULT.algal_pct(:,:,2)));
    macroUprightFleshy(simul,:,:) = squeeze(cat(4,OUTPUTS(simul).RESULT.algal_pct(:,:,3)));
      
    if OPTIONS.doing_COTS == 1
        % Assuming 0.6 CoTS per grid ~ 0.22 CoTS per tow
        % (0.22 per tow is equivalent to 1500 COTS per km2 (Moran & De'ath 92), so that 1 COTS per grid (400m2) is equivalent to 0.22*2500/1500
        COTS_mantatow(simul,:,:) = (0.22/0.6)*squeeze(cat(4,OUTPUTS(simul).RESULT.COTS_total_perceived_density));
        COTS_densities0(simul,:,:,:) = squeeze(OUTPUTS(simul).RESULT.COTS_all_densities); % Density for 400m2
        COTS_settler_density(simul,:,:) = squeeze(OUTPUTS(simul).RESULT.COTS_settler_densities); % Density for 400m2
        COTS_larval_supply(simul,:,:) = squeeze(OUTPUTS(simul).RESULT.COTS_larval_supply); % Density for 400m2
        COTS_larval_output(simul,:,:) = squeeze(OUTPUTS(simul).RESULT.COTS_larval_output); % Density for 400m2
    end
    
    if OPTIONS.doing_COTS_control == 1 && META.nb_time_steps > 23
        COTS_CONTROL_culled_reefs(simul,:,:) = squeeze(OUTPUTS(simul).RESULT.COTS_culled_reefs);
        COTS_CONTROL_remaining_dives(simul,:) = squeeze(OUTPUTS(simul).RESULT.COTS_control_remaining_dives);
        COTS_CONTROL_culled_density(simul,:,:) = squeeze(OUTPUTS(simul).RESULT.COTS_culled_density);
    end

end

% New (08/2021): only record COTS densities by yearly classes to reduce output size
% COTS_densities = COTS_densities0(:,:,:,1:2:end)+COTS_densities0(:,:,:,2:2:end);
% 09/2021: now just sum across all juveniles and across all adults
if OPTIONS.doing_COTS == 1
    COTS_juv_densities = sum(COTS_densities0(:,:,:,1:(META.COTS_adult_min_age-1)),4);
    COTS_adult_densities = sum(COTS_densities0(:,:,:,META.COTS_adult_min_age:end),4);
end

clear OUTPUTS simul s COTS_densities0 NB_TIME_STEPS

save (OPTIONS.OutputFileName)
