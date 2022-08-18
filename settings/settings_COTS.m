%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Y.-M. Bozec, MSEL, created March 2017.
%
% Parametrisation of COTS demographics and consumption
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

META.COTS_maximum_age = 16 ; % maximum age in number of time steps of 6 months

%%% COTS MORTALITY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('COTS_mortalities.mat') % Death fraction for each age (6mo age classes) from meta-analysis (YM)
% The first class gives mortality to apply to settlers in order to estimate 6mo old recruits
% Calculate assuming representative age is the median of each class (3mo, 9mo old etc.)
META.COTS_mortality = COTS_6mo_mortalities(1:META.COTS_maximum_age);
META.COTS_mortality(5:end) = META.COTS_mortality(5:end)/6; % AFTER CALIBRATION WITH LIZARD IS

META.COTS_outbreak_duration = [4,6] ; % min and max duration of an outbreak (random duration within this range)
% META.COTS_outbreak_duration = 2 ; % fixed duration to be used to reproduce Lizard outbreak (Pratchett)
% META.COTS_density_threshold_for_disease = 0.6 ;  % per 400m2 (~0.22 per tow)
META.COTS_density_threshold_for_disease = 2.7 ;  % per 400m2 (~1 per tow)
% A CoTS population is considered outbreaking when the density of 18 mo+ starfish reaches 0.3 individuals per 200 m2 
% (Keesing and Lucas 1992, Moran and De’Ath 1992). If adult COTS density is above this threshold for the specified duration
% then population crashes due to disease  

META.COTS_dieback_classes = 4:META.COTS_maximum_age; % CoTS classes to reset to the background density
% After the end of an outbreak (disease or starvation) the sum of all these
% classes eqauls the background density

META.COTS_background_density = 0.01 ; %total density (nb CoTS per grid) of CoTS, excluding recruits
% CoTS population is reset to this number if the cover of preferential coral preys is too low
META.COTS_background_density_INIT_BOX = 0.05; % Background density in the CoTS initiation box

%%% SIZE AND FECUNDITY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
COTS_ages = ((6:6:META.COTS_maximum_age*6)-3)/12 ; % in years (representative)

% Estimate diameter (in mm) of COTS ofr every age
META.COTS_diameters = (-0.6717*(COTS_ages.^2)+12.638*COTS_ages - 2.1623)*10 ; % in mm (Engelhardt et al. 1999)
% Refine settlers and juveniles
META.COTS_diameters(1)= 3.6; % from Wilmes et al. (2017)
META.COTS_diameters(2)= 16.8; % from Wilmes et al. (2017)

% Estimate female fecundity (millions eggs produced for every COT)
Weights = 6.29*(1e-5)*META.COTS_diameters.^2.929; %(Kettle & Lucas 1987)
Weights(Weights>5000)=5000 ; % cap weight (heaviest found by Kettle and Lucas 1987 was 4.5kg
META.COTS_fecundity = 558*(Weights.^1.439); % female gamete production million eggs (Kettle & Lucas 1987)
META.COTS_fecundity(1:4) = 0; % not mature before 2yr+ (200-250mm, Lucas 1984).

META.COTS_adult_min_age = 4; % Youngest adults to include in the calculation of perceived (mantatows) total adult density
% For COTS control essentially. Class #4 (=subadult1) has a median size of 179mm with detectability of 0.54

%%% COTS CONSUMPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
META.COTS_feeding_prefs = [0.12 ; 0.393 ; 0.258 ; 0.12 ; 0.084 ; 0.025]; % from Karlo

% Per capita consumption rate, in cm2 of corals consumed over 6 month
% COTS start feeding on corals when they are 1yr+
Summer1 = 30.4*6*(17.4*(META.COTS_diameters/10)-355) ; % Oct 87 Keesing and Lucas 1992
Summer2 = 30.4*6*(10.4*(META.COTS_diameters/10)-112) ; % Oct 88 Keesing and Lucas 1992
Summer3 = 30.4*6*(8.6*(META.COTS_diameters/10)-112) ; % Jan 88 Keesing and Lucas 1992
% Winter =  30.4*6*(7.7*(META.COTS_diameters/10)-214) ; % Jun 88 Keesing and Lucas 1992

COTS_feeding_rate_summer = (Summer1 + Summer2 + Summer3)/3 ;
COTS_feeding_rate_winter = COTS_feeding_rate_summer/2; % COTS eat ~twice less in winter (Keesing and Lucas 1992)

META.COTS_feeding_rates = [COTS_feeding_rate_summer ; COTS_feeding_rate_winter]; 
META.COTS_feeding_rates(META.COTS_feeding_rates<0)=0;

META.COTS_feeding_rates = 1.8*META.COTS_feeding_rates; % AFTER CALIBRATION WITH LIZARD IS


%% Parametrisation of COTS outbreaks and connectivity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ID of preferred corals to evaluate the coral threshold
META.COTS_pref_corals = [1:4] ; % as total of corals #1 to #3 (or #4)

%% COTS detectability (extracted from MacNeil et al (2016))
load('COTS_detectability.mat')
COTS_detectability = [COTS_detectability ; 0.98*ones(10,1)]; %just extending to extra classes (trimmed below)
META.COTS_detectability = COTS_detectability(1:META.COTS_maximum_age)';

META.doing_COTS_connectivity = 1 ;

META.COTS_immigration = 1e6*ones(1, META.nb_time_steps) ; % Forced larval input for a 400m2 area when connectivity is OFF
% only works if META.doing_COTS_connectivity = 0 (impact on COTS settlement depends on META.COTS_BH_beta)

% Parameter a of the B-H function (same for all reefs)
META.COTS_BH_alpha = 4*1e4 ; % which is 100 settlers per m2 (before mortality)
% Gives the following densities fro 200m-2:
% ~350 recruits (0-6 month old), ~71 juveniles (<15 cm), ~10 sub-adults (15–25 cm) and ~32 adults (>25 cm), 
META.COTS_BH_beta = 0.5*1e7 ; % fitted to Lizard Is (Pratchett 2005)

% Minimum survival rate of CoTS larvae before dispersal with and without
% Chlorophyll forcing. If no forcing with Chlorophyll, all reefs are set to this value
META.COTS_min_larval_survival = 0.01;

% Force self-seeding of COTS larvae
META.COTS_min_selfseed = 0.28 ; % same as for corals

% For the COTS_control module (must be 1 if META.COTS_min_selfseed>0):
META.force_COTS_selfseed = 1 ; % 

% Minimum proportional coral cover below which CoTS population crashes (end of the outbreak)
META.COTS_coral_threshold = 0.05 ; % as total cover of preferred corals (coral types defined in META.COTS_pref_corals)

% Placeholder for future implementation (needs to be synchronized with WQ layers)
META.connectivity_regime = 1 ; %determine which matrices to use for both coral and COTS

%% Population structure of COTS for initialisation (relative abundance)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Because recruitment is seasonal, at a given step a cohort is made of densities in 1 every 2 age class.
% With constant recruitment, 15 iterations (7.5 years) is enough to achieve equilibrium with 50 recruits. 
% We retain the age distri in winter, so there is no individuals in the recruit class.
% The resulting relative age distri is independent to the number of recruits.

iter = size(META.COTS_mortality,1); % number of iterations for desired stage of population build-up
% 1 iter = 6 month; starts with summer (recruitment) and finishes in summer
R1 = repmat([100 0], 1, iter/2) ; % first value for summer - here 50 settlers every summer for 5 years
COTS_pop_init = zeros(length(R1)+1,size(META.COTS_mortality,1)); 

for t = 2:(iter+1)
    temp = [ R1(t-1) COTS_pop_init(t-1,1:(end-1)) ] ;
    COTS_pop_init(t,:) = (1-META.COTS_mortality').*temp ;
end

% Keep the two last predictions of pop structure -> first row = pop in summer, second = winter
X = COTS_pop_init((end-1):end,:); 
rel_X_sub_and_adults = X(:,3:end)./sum(X(:,3:end),2); 

% Store for use during each season:
META.COTS_init_age_distri_OUTBREAK = zeros(2,size(COTS_pop_init,2));
% Correct for imperfect detectability from manta tow surveys
META.COTS_init_age_distri_OUTBREAK(:,3:end) = rel_X_sub_and_adults./META.COTS_detectability(3:end);
% Extrapolate the expected juveniles (class 2) in winter:
META.COTS_init_age_distri_OUTBREAK(2,2) = X(2,2)*sum(META.COTS_init_age_distri_OUTBREAK(2,3:end));
% Note: no recruits in winter; in summer, nb of recruits is estimated by the model from fecondation at
% the previous time steps so here left blank

% Store for use during each season:
Xref = zeros(2,size(COTS_pop_init,2));
META.COTS_ref_age_distri = zeros(2,size(COTS_pop_init,2));
% Scale to an observed density of 1
Xref(1,:) = X(1,:)/sum(X(1,META.COTS_detectability>0));
Xref(2,:) = X(2,:)/sum(X(2,META.COTS_detectability>0));
% Correct for imperfect detectability
META.COTS_ref_age_distri(1,:) = Xref(1,:)/sum(Xref(1,3:end).*META.COTS_detectability(3:end));
META.COTS_ref_age_distri(2,:) = Xref(2,:)/sum(Xref(2,3:end).*META.COTS_detectability(3:end));
% This gives the density of representative CoTS pop for an observed total density of 1

clear iter reef R1 R2 COTS_init_age_distri_INCIPIENT COTS_init_age_distri_OUTBREAK t n X rel_X_sub_and_adults temp COTS_pop_init
clear COTS_ages Weights COTS_feeding_rate_summer COTS_feeding_rate_winter COTS_detectability COTS_6mo_mortalities
clear Summer1 Summer2 Summer3 Winter Xref
