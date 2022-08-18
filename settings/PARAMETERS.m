% Y.-M. Bozec, MSEL, created Nov 2011.
%
% GENERATES DEFAULT PARAMETERS
%
% Better to keep these values unchanged here and to refine specific parametrisation in MAIN
%_________________________________________________________________________________________
%_________________________________________________________________________________________
%
%       SET UP META-LEVEL PARAMETERS
%_________________________________________________________________________________________

META.nb_reefs = 1 ; % Number of populations (= reefs)

META.nb_time_steps = 60 ; % Number of time steps (= 6 mo)

META.nb_simul = 1 ; % Number of simulation runs to repeat

META.max_colonies = 40 ; % Maximum number of colonies of any given species in a cell

META.recruitment_type = 0 ; % set to 0 for fixed, 1 for variable. Refined with actual free space.

META.space_limited_settlement = 1 ; % set to 0 for no effect of available space on the probability of recruitment

META.grid_x_count = 20; % number of grid cells along the x-edge
META.grid_y_count = 20; % number of grid cells along the y-edge

META.cell_x_size = 100 ; % size of a single cell in centimeters along its x-edge
META.cell_y_size = 100 ; % size of a single cell in centimeters along its y-edge

META.track_populations = 0 ;% set to 1 for tracking the size of every single coral colony (for John Hedley reef viz)

META.doing_size_frequency = 0 ; % set to 1 to estimate the number of colonies in different size classes

META.doing_coral_competition = 0 ; % For OA, based on Nick's experiments in Moorea (requires update of f_coral_growth.m)

META.doing_parrotfish_fishing = 0 ; % set to 1 to allow graing impact changing as a response to fishing
% NOTE this requires running f_parrotfish_dynamics to estimate grazing over time

META.doing_COTS = 0 ; % set to 1 for simulating coral predation by the crown-of-thorn starfish (Pacific)

META.doing_coral_connectivity = 0 ; % set to 0/1 for using/not using a connectivity matrix

META.doing_COTS_connectivity = 0 ; % set 1 to simulate COTS larval connectivity (requires connectivity matrices)

META.isverbose = 0 ; % set to 0 to stop displaying setup otpions each time the model runs 

META.reef_ID = 1;

META.doing_water_quality = 0;

META.randomize_WQ_chronology = 0;

META.doing_clades = 0 ;

META.doing_genetics = 0;

%_________________________________________________________________________________________
%
%       SET UP REEF PARAMETERS
%_________________________________________________________________________________________

REEF.exposure = 0 ; % Only affects algal dynamics (growth rate)
% (replaces the previous 'is_leeward' that was affecting Dictyota)
% proportion of exposure between 0 (leeward glovers) and 1 (windward glovers).
% exposure used to slow down rate of growth of Dictyota to maximum of 43% of windward value (when exposure =0);

REEF.dictyota_declines_seasonally = 0 ; % set to 0/1 to switch off/on Dictyota dying back in winter

REEF.nongrazable_substratum = 0.1 ; % non-grazeable cover

REEF.herbivory = 0.33 ; % total amount of algae that can be grazed by fish per time step
%expressed as a proportion of the reef area - now acts as a constraint to max amount of seabed

REEF.diadema = 0 ; % total amount of algae that can be grazed by Diadema per time step
%expressed as a proportion of the reef area.

%_________________________________________________________________________________________
%
%       SET UP CORALS
%_________________________________________________________________________________________

%%%%% 1) SPECIES-SPECIFIC PARAMETERS (need as many values there are coral types)
% Code indicating coral types that are brooder (1) or spawner (0)
CORAL.is_brooder = [ 0 ; 0 ; 0 ; 0 ; 0 ;0 ];

% Initial coral covers (as proportions)
REEF.initial_coral_cover = [ 0.05 ; 0.05; 0.05 ; 0.05 ; 0.05 ; 0.05];  % 3% total cover

% Horizontal growth rate (radius extension) in cm per 6 month
CORAL.growth_rate = [ 2.2 ; 2.2 ; 1.5 ; 1.2 ; 0.4 ; 0.6]; % NEW (2018)
% YM: following Ken's suggestion, horizontal extension is the cosinus of 45
% degree angle between branch elongation and the reef substrate
% (cos(45)*4.1~2.2 for arborescent)
% YM: coral group definition (since 2018):
%     1. ARBORESCENT (STAGHORN) CORALS (e.g., Acropora muricata, Acropora nobilis, Acropora robusta);
%     2. PLATING CORALS (e.g., Acropora hyacinthus, Acropora cytherea);
%     3. CORYMBOSE/SMALL BRANCHING ACROPORIDS (e.g., Acropora millepora, Acropora humilis);
%     4. POCILLOPORIDS/NON-ACROPORID CORYMBOSE CORALS (e.g., Stylophora pistillata);
%     5. SMALL MASSIVE/SUBMASSIVE/ENCRUSTING (e.g., Lobophylliidae, favids, Goniastrea);
%     6. LARGE MASSIVE (Porites lutea, Porites lobata, Porites australiensis).

% FOR MOOREA (Nick):
% -> 1.75 for Pocillopora spp. (based on verrucosa)
% -> 2.5 for Acropora spp. (based on young A. hyacinthus and A. nasuta)

% NEW (2018):
CORAL.juvenile_growth_rate = 0.5; % radial extension in cm per 6 month from Doropoulos et al. (2015) - same for all species. 

% CURRENTLY NOT USED:
% surface index to convert planar into actual surface area of corals for COTS consumption:
CORAL.SI = [ 6.16 ; 5.43 ; 5.65 ; 5.12 ; 3.7/2 ; 3.7 ] ; % Holmes et al. 2008, Keesing 1990

%%%%%%% NEW (08/2015): Proportion of recruits from each coral group (sums to 1)
CORAL.prop_settlers = [ 0.05 ; 0.25 ; 0.25 ; 0.15 ; 0.15 ; 0.15]; % ~average props on midshelf/outer reefs (Jez)

% Maximum diameter of coral colonies to estimate max size (cm2)
CORAL.max_diameter = [120 ; 100 ; 50 ; 40; 60 ; 200 ]; % Further limited by cell size

% Define the initial size distribution of coral colonies
% 1 for a random selection of corals of different sizes - i.e. equal probability per size
% 2 for lognormal distribution (based on Meesters et al. 2001 in curacao)
CORAL.size_model = 2*[ 1 ; 1 ; 1 ; 1 ; 1 ; 1 ]; 

% Below are parameters for generating coral colonies at initial steps based
% on a lognormal distribution (using Meesters et al 2001 estimates).
% Note that generated sizes will be further limited by CORAL.max_diameter and META.cell_area_cm2
% (the size of a colony cannot exceed the size of a cell). 
% Note that CORAL.size_mean / CORAL.size_var can be used to parameterize other kind of distribution
% given a new model is specified in f_derivecoralcover

% Mean of colony size (cm2)
% CORAL.size_mean = [log(35) ; log(314) ; log(35) ; log(773)]; % BASED ON COZUMEL SIZE DISTRI (planar area) - but small!!
CORAL.size_mean = [log(700) ; log(700) ; log(500) ; log(700); log(700) ; log(700) ]; %For GBRF
% Meesters et al. (2001) - geometric means mu of surface area (cm2)
% Agaricia agaricites: 70 / Porites astreoides: 50 / Siderastrea siderea: 470
% Montastrea cavernosa: 770 / Mont. faveolata: 3100 / Mont. annularis: 1800 ; Porites porites: NA
% Mycetophyllia spp. NA / Meandrina meandrites (~)100
% -> size follows a lognormal distribution of mean log(mu)
% -> we use M. meandrites size parameters for coral #3 (normally mycetophyllia)
% WARNING: Meesters coral sizes are actual surface areas, not planimetric. Ask Erik for planimetric values.

% Variance of colony size
% CORAL.size_var = [0.92 ; 1.58 ; 0.92 ; 1.49]; % BASED ON COZUMEL SIZE DISTRI (planar area) - but small!!
% CORAL.size_var = [log(1.6^2) ; log(1.8^2) ; log(1.9^2) ; log(1.7^2) ; log(1.9^2) ; log(1.7^2)]; % BASED ON CURACAO - Meesters (surface area)
CORAL.size_var = [log(2^2) ; log(2^2) ; log(2^2) ; log(2^2) ; log(2^2) ; log(1.8^2)];  %For GBRF

% Meesters et al. (2001) - assumed to be geometric SD (cm2)
% Agaricia agaricites: 2.2 / Porites astreoides: 2.0 / Siderastrea siderea: 2.5
% Montastrea cavernosa: 1.7 / Mont. faveolata: 1.7 / Mont. annularis: 1.7 ; Porites porites: NA
% Mycetophyllia spp. NA / Meandrina meandrites 1.9
% -> size follows a lognormal distribution of variance log(SD^2)
% -> we use M. meandrites size parameters for coral #3 (normally mycetophyllia)
% WARNING: Meesters coral sizes are actual surface areas, not planimetric. Ask Erik for planimetric values.

% Rate at which coral recedes when in contact with Lobophora;
% 4 from Agaricia / Porites; % 0.75 from Meandrina;
% 0.75 from Mycetophyllia (but use value from Agaricia if assume few Mycetophyllia);
% 0.5 from Montastraea
CORAL.lobophora_reduce_rate = [ 0.75 ; 0.75 ; 0.75 ; 0.75 ; 0.75 ; 0.75]/7; % Ortiz et al 2014

% Rate at which coral recedes when in contact with Dictyota;
% = 0.25*6 months from Lirman 2001 (ranges from 0.25 cm2 per month to 0.43) (replaces DICT_OVER_SPAWNER?)
CORAL.dictyota_reduce_rate = [ 1.5 ; 1.5 ; 1.5 ; 0; 1.5 ; 1.5 ]; % for spawners only


%%%%% 2) THE FOLLOWING CORAL PARAMETERS ARE NOT SPECIES-SPECIFIC %%%%%%%%

% Percent of colony growth that is actively overgrow other corals
% set arbitrarily so that up to 20% of the colony's projected growth can overgrow smaller corals
CORAL.fractional_overgrowth = 0.2 ; 

% Size of full fecundity (OLD Caribbean stuff)
% threshold size for changing mortality as fecunditiy is a direct function of size)
% CORAL.adult_size = 110; % Caribbean: 250; % MOOREA: 50 ;

% 2018: Fecundity is now group-specific, with number of eggs predicted by colony size
CORAL.fecund_min_size = [123; 123 ; 134 ; 31 ; 38 ; 38 ];
% From Hall and Hughes 1996: size above which all colonies are gravid:
% A. hyacinthus: 123 cm2
% A. millepora: 134 cm2
% S. pistillata: 31 cm2
% Goniastrea retiformis: 38 cm2

% Parms for empirical relationship linking colony size to total egg volume (Hall & Hughes 1996)
% Egg volume = exp(a+b*log(size)), further divided by mean egg volume to get number of eggs
CORAL.fecund_a = [1.03 ; 1.03 ; 1.69 ; -1.20 ; 0.86 ; 0.86 ];
CORAL.fecund_b = [1.28 ; 1.28 ; 1.05 ; 2.27 ; 1.21 ; 1.21 ];

% Intercept for incidence relationship from Meesters
CORAL.partial_mortality_inci_int = 88.9;% Meesters et al. 1997 (see Mumby et al. 2013)

% Gradient for incidence relationship from Meesters
CORAL.partial_mortality_inci_gra = -11.2; % % Meesters et al. 1997 (see Mumby et al. 2013)

% Intercept for area relationship from Meesters
CORAL.partial_mortality_area_int = -2.9 ; % Mumby et al. 2013

% Gradient for area relationship from Meesters
CORAL.partial_mortality_area_gra = 1.59 ; % Mumby et al. 2013

%% Chronic whole-colony mortality (see size thresholds thereafter)
CORAL.adult_whole_mortality_rate = 0.01 ; % based on Porites astreoides (Bythell et al 1993)
CORAL.adol_whole_mortality_rate = 0.02 ; % based on 50-200cm2 Porites astreoides (Bythell et al 1993) % MOOREA: 0.04
CORAL.juv_whole_mortality_rate = 0.1 ; % Doropoulos et al (2016)

% New stuff (August 2013) -> species-specific sensitivity to natural mortality
CORAL.sensitivity_whole_natural = [0.1 ; 5 ; 4.5 ; 6 ; 1; 1] ; % to adjust CORAL.adult_whole_mortality_rate to Pacific corals (Ortiz et al 2014)
CORAL.sensitivity_partial_natural = [1 ; 1 ; 1 ; 1 ; 1; 1] ; % placeholder for when parmameters are available
CORAL.extent_partial_natural = [2 ; 0.4 ; 1.4 ; 1.8 ; 1 ; 1] ; % multiplicator of the extent of area lost (Ortiz et al 2014)

%% for GBR (2019)
CORAL.size_threshold_wcm = 60; % Size threshold for whole-colony mortality due to cyclones and bleaching (Edwards et al. 2011)

CORAL.size_threshold_pm = 250; % Size threshold for partial mortality (natural and due to cyclones) (Edwards et al. 2011)

CORAL.juv_max_size = 13 ; % at least 2 years old, escaped the post-settlement bottlenecks (Doropoulos et al. 2016)
% Used for cyclone mortality, natural mortality and growth

% Only for tracking size distribution: threshold diameters for transitions juv -> adol -> adults
CORAL.juv_max_diam = floor(sqrt(CORAL.juv_max_size/pi)*2); % maximum diameter of juveniles
CORAL.adol_max_diam = floor(sqrt(CORAL.size_threshold_pm/pi)*2); % maximum diameter of adolescents = lower bound of mature corals in the Caribbean 
CORAL.adult_max_diam = sqrt(META.cell_x_size*META.cell_y_size); % defined by cell size

% Define bin size for the size distributions for (1) juvenile, (2) adolescent and (3) adult colonies
CORAL.size_bins = [3 ; 20 ; 500]; % in cm2
% This gives the following bin edges 
% Juvenile (cm2):   0     3     6     9    12
% Adolescents (cm2):  12    32    52    72    92   112   132   152   172   192   212   232   240
% So to get juveniles <= 5cm take the first two size classes of adolescents (0 -10, 10-30cm2)
% Note this includes recruits (1cm2)
CORAL.diam_bins = [1 ; 1 ; 10]; % in cm
% This gives the following bin edges for juvenile bins (cm): 0  1  2  3  4 (4 size classes)
% Adolescents (cm):  4  5  6  7  8  9  10  11  12  13  14  15  16  17 (13 size classes)
% Adults (cm): 17 27 37 47 57 67 77 87 97 100 (9 size classes)

%% UPDATE March 2022 to optimise the design of colony size distribution for less memory usage
CORAL.juv_max_diam = 5; % every diameter <5cm is a juvenile
CORAL.diam_bins = [2 ; 4 ; 14]; % in cm
% This gives the following bin edges for juvenile bins (cm): 1  3  5 (2 size classes)
% Adolescents (cm):  5  9  13  17 (3 size classes)
% Adults (cm): 17 31 45 59 73 87 101 (6 size classes)

%__________________________________________________________________________________________
%
%       SET UP ALGAE - (FEB 2017) NEW PARAMETRISATION
%__________________________________________________________________________________________
% Initial algal covers in the following order after Bozec et al. 2019
% (1) EAM ; (2) DICTYOTA ; (3) LOBOPHORA ; (4) THICK TURF
REEF.initial_algal_cover = [ 0 ; 0.05 ; 0.05 ; 0 ]; % Turf/EAM are estimated during iniitalisation automatically

% Relative proportion of consumption of the different algal types by fish herbivory
ALGAL.herbivory_props = [ 0 ; 0.1 ; 0.05 ; 0.85 ]; % Palau

% Rules for switching among algal preys
ALGAL.feeding_prefs = [2 ; 3 ; 4 ; 1 ];

% Relative proportion of consumption of the different algal types by Diadema
ALGAL.diadema_props = [ 0.8 ; 0.2 ; 0 ; 0 ]; % TO RE-VISIT (Caribbean stuff)

% Proportional reduction in algal growth rate from contact with coral (actually used as 1-ALGAL.coral_reduce_macrogrowth)
ALGAL.coral_reduce_macrogrowth = 0.25 ;

% Local neighborhoud of algal that smothers adolescent corals
ALGAL.critical_algal_contact = 0.40 ;

% use 0.8 but could drop to 50% based on Foster et al
ALGAL.vcritical_algal_contact = 0.80 ;

% Effect of macroalgae on growth rate of coral recruits (based on Box)
ALGAL.macroalgal_coral_recruit_growth_rate = 0.3 ;

% Effect of macroalgae on growth rate of adol and adult corals (based on Lirman for Dicty vs Ag and Pa)
ALGAL.macroalgal_coral_growth_rate = 0.1 ;

% Intrinsic growth rate (logistic)
% Will be collated in INITIALISATION under 'ALGAL.growth_rate"
algal_growth_rate_summer = [ 0 ; 0.843 ; 0.598 ; 0 ] ; % month-1 !!! Palau (Lighthouse, high productivity site)
algal_growth_rate_winter = [ 0 ; 0.843 ; 0.598 ; 0 ] ; % month-1 !!! Palau (Lighthouse, high productivity site)

ALGAL.nb_step_algal_dynamics = 6 ; % run internally within a model step (coherent with growth rates)
META.convert_to_canopy = 1 ; % Multiplicator of DICT canopy - use 1 if Dict does not overtop Lob
ALGAL.settlement_lob_cm2 = 2.5; % Monthly settlement of Lobophora over 1m2 (Diaz-Pulido & McCook 2004)
ALGAL.margins_lob_overgrowth_cm2 = 45; % Monthly overgrowth over 1m2 from neighbouring cells with 60% Lob (de Ruyter van Steveninck & Breeman 1987)

%__________________________________________________________________________________________
%
%       SET UP HURRICANES
%__________________________________________________________________________________________

META.doing_hurricanes = 0 ; % set to 0/1 to switch off/on hurricanes

META.randomize_hurricane_chronology = 0 ; % set to 1 to randomize hurricane chronology for each simulation

META.random_hurricanes = 1 ; % if no prescribed scenario, just random hurricanes

REEF.hurr_strike_proba = 0.1; % Frequency at which hurricanes occur in summer, if no prescribed hurricane scenario (eg, Cozumel)

META.hurricane_effect_on_macroalgae = 0.9 ; % proportion of macroalgae removed by a hurricane

META.hurricane_effect_on_recruits = 0.8 ; % prob that recruits die if less than 10 cm diameter

C_adjust_GBR = 5 ; % calibration of cyclone mortality with AIMS LTMP transect data
CORAL.sensitivity_hurricane = C_adjust_GBR*[1.3 ; 1.3 ; 1.5 ; 1.1 ; 0.6 ; 0.5]; % relative proportions of damages from AIMS LTMP

% CORAL.extent_hurricane = [1 ; 1 ; 1 ; 1 ; 1 ; 1]; % Placeholder multiplicator affecting the extent of partial mortality due to hurricanes

META.hurri_whole_mort_dead = 0; % Multiplicator of CORAL.hurr_mortality_rate for the dead colonies
% 0 produces no effect of hurricane on dead colonies (they are not removed)
% 1 applies the same relative effect than specified for living colonies
% 2 doubles the relative effect for every species
% etc.

%__________________________________________________________________________________________
%
%       SET UP CORAL BLEACHING
%__________________________________________________________________________________________

META.doing_bleaching = 0 ; % set to 0/1 to switch off/on bleaching

CORAL.sensitivity_bleaching = [1.5 ; 1.6 ; 1.4 ; 1.7 ; 0.25 ; 0.25]; % GBR after calibration with Hughes et al (2018)
% A multiplicator of the bleaching mortality probabilities (whole and partial colonies)

CORAL.bleaching_depth = 1;  % coefficient to account for reduction in bleaching with depth due to light attentation (Baird et al. 2018)
% Use 1 to get bleaching mortality observed by Hughes et al. (2018) at ~2m -> used for the GBR hindcast (Bozec et al. 2022)
% Use 0.5 to extrapolate to ~7m depth according to Baird et al. (2018)

CORAL.bleaching_partial_extent = [0.05 ; 0.05 ; 0.05 ; 0.05 ; 0.4 ; 0.2]; % GBR from Baird and Marshall (2002) - minimal for branching
% Extent of lesions as a proportion of colony size

% Only if doing clades
CORAL.clade_prop = 0.95 ; % proportion of clade thermally sensitive over tolerant at initial step
CORAL.clade_reduced_growth = 0.5 ; % proportional reduction of coral growth rate when clade 2

CORAL.proba_switching = [0.5 ; 0.5 ; 0.5 ; 0.5 ; 0.5 ; 0.5]; % probability of a switching clade from sensitive to tolerant after a bleaching event
% NOTE THIS IS INDEPENDENT OF BLEACHING MORTALITY

CORAL.bleaching_tolerance_clade = 0.29 ; % Reduces mortality due to bleaching if clade 2
% NOTE this now replaces the reduction in mortality due to past bleaching experience

% Coeff for adjusting Hughes' initial bleaching mortality to
% long-term (6mo) mortality (based on calibration with post-bleaching cover loss)
META.bleaching_whole_offset = 6 ; % extrapolates Hughes' initial mortality to whole colony mortality for the entire event
META.bleaching_partial_offset = 6 ; % extrapolates Hughes' initial mortality to partial mortality for the entire event

%__________________________________________________________________________________________
%
%       SET UP 3D (FOR RUGOSITY - COZUMEL, not in the GBR model yet
%_________________________________________________________________________________________

META.doing_3D = 0 ; % set to 0/1 to switch off/on coral 'calcification'

REEF.initial_rugosity = 1.5 ; % Rugosity of the non-living substrate for initialisation
% Actual rugosity will be higher after adding living coral cover, so that several trials are
% be necessary to obtain the desired reef rugosity. 
% Reef dominated by Montastraea: R~1.2 for coral cover=5% (Alvarez-Filip et al. 2011),
% then add 0.2 units for every 5 additional percent of coral cover

CORAL.min_height = 3; % minimum height of dead coral colonies after erosion (when below, colony is removed)
CORAL.min_diameter = 2 ; % minimum diameter of dead coral colonies (when below, colony is ignored, 
% and exposes the underlying substrate)

%__________________________________________________________________________________________
%
%       SET UP RANDOMIZATION
%_________________________________________________________________________________________

META.randomize_inputs = 0 ; % set to 1 to randomize input parameters (see f_randomize_inputs.m)
% If ==1, the standard deviations below must be valued

% Standard deviation of initial coral cover for randomization
% REEF.initial_coral_cover_SD2 = [ 0 ; 0 ; 0 ; 0 ; 0 ; 0 ]; % null because no randomization by default -> needs parametrisation
% 
% REEF.nongrazable_substratum_SD2 = 0 ; % standard deviation of the non-grazable cover
% 
% CORAL.growth_rate_SD2 = [ 0 ; 0 ; 0 ; 0 ; 0 ; 0]; % Standard deviation of coral growth rate
% 
% REEF.herbivory_SD2 = 0 ; % Standard deviation of hebivory
% 
% REEF.fish_bioerosion_SD2 = 0; % for randomization; null by default

%__________________________________________________________________________________________
%
%       SET UP CALCIFICATION *
%_________________________________________________________________________________________
% * For the moment only estimates changes in coral calcification due to SST rises
% Parameters for carbonate budget may be implemented later

META.doing_calcification = 0 ; % set to 0/1 to switch off/on coral 'calcification'

%__________________________________________________________________________________________
%
%       SET UP BIOEROSION (COZUMEL)
%_________________________________________________________________________________________
% Bioeorsion on dead substratum (mixing past and recently dead)
% REEF.fish_bioerosion = 0.1989*9.33*365/(2*10000); % in cm3/cm2/6 month
% 0.1989 cm3/m2/hour obtained from  extractdataforgrazers(SWC_unfished_dat,SWC_unfished_txt)
% 9.33 is the number of daylight hours (see Pete_bioerosion2.m) and /10000 to convert in cm2
% -> gives 0.0339 cm3/cm2/6mo
% Mallela & Perry (2007): 20 g CaCO3/m2/yr ->  20/(1.5*10000) = 0.0013 cm3/cm2/year
 
% REEF.sponge_bioerosion = 0.507*1000/(10000*2*1.5);
% 0.507 kg/m2/year which is erosion rate on dead colonies (see Pete_bioerosion2.m)
% 1.5 is coral density (g/cm3) (average from Pete_bioerosion2.m). Use the species specific densities?
% 1000 converts kg to g
% /10000 to convert in cm2
%NOTE this produces bioerosion rate for fish much higher than for sponges(3 times higher)


%__________________________________________________________________________________________
%
%       SET UP CORAL DISEASES - not implemented yet (but see 2010 Caribbean model)
%__________________________________________________________________________________________

META.doing_disease = 0 ; % set to 0/1 to switch off/on diseases

% specific parametrisation has to be done in settings_DISEASES.m
