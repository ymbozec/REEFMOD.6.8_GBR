% REEFMOD-GBR restoration settings
%
% Yves-Marie Bozec, y.bozec@uq.edu.au, 09/2019 (RRAP investment case)
%
% Re-worked in 08/2021 for greater flexibility + extended options
%__________________________________________________________________________

%% 1) Options for coral deployment
% Restoration effort in terms of total number of outplants available at each time step
META.total_nb_outplants = RESTORATION.total_nb_outplants;
% Restoration effort in terms of number of reefs where coral outplanting is undertaken at each time step
META.nb_reefs_outplanted = RESTORATION.nb_reefs_outplanted ;
% Timing of effort
META.doing_coral_outplanting = RESTORATION.doing_coral_outplanting ; % timing of coral deployment
META.threshold_for_outplanting.min = 2.5 ; % Maximum percent coral cover (all corals) under which reef is selected for outplanting
META.threshold_for_outplanting.max = 15 ; % Maximum percent coral cover (all corals) under which reef is selected for outplanting

META.restore_random_density = 1 ;



% Old way for setting density of outplants (keep the code for now)
META.outplanted_coral_per_cell = zeros(1, META.nb_coral_types);

% Create demographic parameters of coral outplants as duplicate of their native counterparts
% This extend specific fields in the structure array CORAL
% eg, CORAL.growth_rate is initially 6x1 but becomes 8x1 if 2 types of outplants are used (the last 2 values are the added ones)
% The corresponding type is given in META.outplanted_species
if isempty(RESTORATION.outplanted_species)==0
 
    lengths = structfun(@(x) size(x,1), CORAL);
    I=find(lengths==META.nb_coral_types);
    fields = fieldnames(CORAL);
    select_fields = fields(I);
    
    for s = 1:length(META.nb_coral_types)
        
        for f=1:length(select_fields)
            
            V = extractfield(CORAL,cell2mat(select_fields(f))) ;
            newV = [V  V(META.outplanted_species(s))];
            CORAL = setfield(CORAL,cell2mat(select_fields(f)),newV') ;
            % Note that prop_settlers (proportion of each coral type in maximum settlement) summed to 1 initially
            % Now above 1 but probably not a big deal as an addition of a new population?
        end
    end
    
    META.nb_coral_types = 2*META.nb_coral_types; % duplicate all coral types
    META.outplanted_species = zeros(1,META.nb_coral_types);
    META.outplanted_coral_diameter = META.outplanted_species;
    META.outplant_species_prop = META.outplanted_species;
    
    META.outplanted_species(RESTORATION.outplanted_species+5)=1;
    META.outplanted_coral_diameter(RESTORATION.outplanted_species+5) = RESTORATION.outplant_diameters;
    META.outplant_species_prop(RESTORATION.outplanted_species+5) = RESTORATION.outplant_species_prop;
    
%     META.nb_coral_types = META.nb_coral_types + length(META.outplanted_species); 
%     META.outplanted_coral_per_cell(META.outplanted_species + 5) = RESTORATION.outplant_densities;
%     META.outplanted_coral_diameter(META.outplanted_species + 5) = RESTORATION.outplant_diameters;
%     META.outplant_species_prop(META.outplanted_species + 5) = RESTORATION.outplant_species_prop ;
end

% New way: density decreases linearly with total coral cover (from 10/m2 to 1/m2 between min and max thresholds)
META.outplanted_density_aquaculture.slope = -0.72 ; % slope of the deployed density/coral cover relationship
META.outplanted_density_aquaculture.intercept = 11.8 ;

%% REFINE HERE SPECIFIC OUTPLANT DEMOGRAPHICS
% CORAL.growth_rate(7)= 2;
% CORAL.sensitivity_bleaching(7) = CORAL.sensitivity_bleaching(2)/2;
% CORAL.sensitivity_bleaching(8) = CORAL.sensitivity_bleaching(3)/2;
CORAL.sensitivity_bleaching(META.outplanted_species+5) = CORAL.sensitivity_bleaching(META.outplanted_species)*exp(-0.35*RESTORATION.DHW_tolerance_outplants);
CORAL.growth_rate(META.outplanted_species+5)=CORAL.growth_rate(META.outplanted_species)*RESTORATION.tradeoff_growth_outplant;

%% DEPLOYMENT OPTIONS
MinDeploymentArea_km2 = 0.01;
MinDeploymentCells = 4 ;

GridSize = META.grid_x_count*META.grid_y_count;
NumberCellsTreated = uint16(GridSize*MinDeploymentArea_km2./META.area_habitat);
DeploymentArea_km2 = MinDeploymentArea_km2*ones(size(NumberCellsTreated));
Reef_ID = META.reef_ID;
Reef_area_km2 = META.area_habitat;

% Create table for deployment characteristics
TMP = table(Reef_ID, Reef_area_km2, DeploymentArea_km2, NumberCellsTreated);
I = find(TMP.Reef_area_km2 < MinDeploymentArea_km2); % find reefs than are smaller than the min deployment area -> will be excluded
J = find(TMP.NumberCellsTreated < MinDeploymentCells); % find reefs than are too big so that nb of deployed cells is too small

TMP.NumberCellsTreated(J) = MinDeploymentCells;
TMP.DeploymentArea_km2(J) = (MinDeploymentCells/GridSize)*TMP.Reef_area_km2(J);
TMP.DeploymentArea_km2(I) = 0;
TMP.NumberCellsTreated(I) = 0;

% Store the final table
META.outplanting_deployment = TMP;

%% 2) Options for rubble stabilisation
% Restoration effort: number of reefs where rubble is stabilised at each time step
META.nb_reefs_stabilised = RESTORATION.nb_reefs_stabilised ;
% Timing of effort
META.doing_rubble_stabilisation = RESTORATION.doing_rubble_stabilisation ; % timing of rubble stab
META.proportion_rubble_stabilised = 0.5 ; % proportion of rubble that will be stabilised on a reef
META.threshold_for_stabilisation = 2 ; % Minimum percent rubble cover above which reef is selected for rubble stabilisation

%% 3) Options for larval enrichment (TODO)
% Restoration effort: number of reefs where larval enrichment is undertaken at each time step
META.nb_reefs_enriched = RESTORATION.nb_reefs_enriched ;
% Timing of effort
META.doing_larval_enrichment = RESTORATION.doing_larval_enrichment ; % timing of rubble stab
META.threshold_for_larval_enrichment = 10 ; % Maximum percent coral cover (all corals) under which reef is selected for larval enrichment
% Set the proportion of the different coral types in the pool of larvae delivered
META.prop_enriched_larvae = [ 0.05  0.25  0.25  0.15  0.15  0.15]; % ~average props on midshelf/outer reefs (Jez) ;
% Set the total amount of larvae delivered
META.enriched_larvae = 1e6*META.prop_enriched_larvae ; %for a 400m2 reef, 

%% 4) Options for Solar Radiation Management (cooling)
META.doing_cooling = RESTORATION.doing_cooling ; % timing of coling
META.cooling_factor = RESTORATION.cooling_factor; % Consider the different levels (RRAP1) [-0.3 ; -0.7 ; -1.3];

%% CREATE PRIORITY LISTS FOR RESTORATION (1 for each technique)
% Currently five different criteria to prioritise reefs for restoration.

% 1) Regional
META.priority_option_Outplant.region = 0;%[1 2 3];
% Takes either 0(no regional priority) or a vector with 1(North) 2(Central) 3(South) in any priority order

% 2) Shelf position
META.priority_option_Outplant.shelf = 0;%[2 3 1];
% Takes either 0(no shelf priority) or a vector with 1(inshore) 2(midshelf) 3(outer) in any priority order

% 3) Likelihood of larval sink as the sum of all inbound link strengths
META.priority_option_Outplant.link_strength = 0 ; %1 (increasing, min external supply first) or %2 (decreasing, max first)
% Takes either 0(no priority), 1(decreasing, highest potential first) or 2(increasing, lowest potential first)

% 4) Likelihood of larval sink as the total number of inbound links
META.priority_option_Outplant.link_number = 0 ; %1 (increasing, min external supply first) or %2 (decreasing, max first)
% Takes either 0(no priority), 1(decreasing, highest potential first) or 2(increasing, lowest potential first)

% 5) Reef size (area)
META.priority_option_Outplant.reef_area = 0 ;
% Takes either 0(no priority), 
% Using GBRMPA reef outline as habitat area: 1(increasing, smallest area first) or 2(decreasing, largest area first)
% Using Gemorphic map as habitat area: 3(increasing, smallest area first) or 4(decreasing, largest area first)

% 5) Focused on a group of reef (from a focal reef, visit all reefs at increasing distance)
META.priority_option_Outplant.focal_reef = 695 ; %Moore Reef

%% Populate same options for Rubble Stabilisation...
META.priority_option_RubbleStab = META.priority_option_Outplant;
%... And refine if necessary:
% META.priority_option_RubbleStab.region = 0 ; % [1 2 3]
% META.priority_option_RubbleStab.shelf = 0 ; % [1 2 3]
% META.priority_option_RubbleStab.link_strength = 0 ;
% META.priority_option_RubbleStab.link_number = 0 ;
% META.priority_option_RubbleStab.reef_area = 0 ;

%% Populate same options for Larval enrichment...
META.priority_option_LarvalEnrich = META.priority_option_Outplant;
%... And refine if necessary:
% META.priority_option_LarvalEnrich.region = 0 ; % [1 2 3]
% META.priority_option_LarvalEnrich.shelf = 0 ; % [1 2 3]
% META.priority_option_LarvalEnrich.link_strength = 0 ;
% META.priority_option_LarvalEnrich.link_number = 0 ; 
% META.priority_option_LarvalEnrich.reef_area = 0 ; 

