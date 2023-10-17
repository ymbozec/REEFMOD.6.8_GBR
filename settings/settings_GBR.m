%__________________________________________________________________________
%
% REEFMOD-GBR combined Hindcast (2008-2020.5) and Forecast (2021-...)
%
% Yves-Marie Bozec, y.bozec@uq.edu.au, 12/2019
%__________________________________________________________________________


%% INITIAL REEF STATES WITH AIMS LTMP TRANSECT AND MANTA TOW DATA
load('AIMS_LTMP_Transect_means.mat') 
% Complete (1992-2017) LTMP transect dataset at 5-6m averaged per reef (2-3 sites per reef).
% Total cover (TOT) + cover of each coral group (G1....G6).
% Reefs assigned to Karlo's polygons based on shortest distance between survey and centroid coordinates
% (see 'RE-ARRANGE_LTMP_TRANSECTS.m'). There are discrepancies in reef names and shelf position 
% between AIMS database and Karlo's list of 3806. using Karlo's designations here for consistency
% 'YEAR_Reefmod' is modelled year (surveys assigned to each season from their actual dates) 

%% 1) Build the average community compo for the corresponding years
survey_year = LTMP_Transects_means.YEAR_Reefmod;
select_LTMP_compo = find(survey_year>=2004 & survey_year<2009); % select transect data from 2004 to 2008 to initialise community compo

TOT_COVER = table2array(LTMP_Transects_means(select_LTMP_compo,9)); % 341 transects surveyed between 2004-08 (123 reefs)
Z = table2array(LTMP_Transects_means(select_LTMP_compo,10:15))./TOT_COVER(:,ones(6,1)); % relative proportion of each group (each row sums to 1)
LTMP_Transects_compos = [LTMP_Transects_means(select_LTMP_compo,[16 17 3]) array2table(Z)]; % Extract AIMS sector (3), not Karlo's (18) as they differ
LTMP_Transects_compos.AIMS_SHELF(LTMP_Transects_compos.SHELF=='I')=1;
LTMP_Transects_compos.AIMS_SHELF(LTMP_Transects_compos.SHELF=='M')=2;
LTMP_Transects_compos.AIMS_SHELF(LTMP_Transects_compos.SHELF=='O')=3;

omean = @(x) mean(x,'omitnan');
% Calculate average community composition (relative props) for each sector(1-11)*shelf positon(1-3)
AIMS_AVERAGECOMPO = table2array(varfun(omean, LTMP_Transects_compos,'GroupingVariables',{'SectorCode','AIMS_SHELF'},...
    'InputVariables',{'Z1','Z2','Z3','Z4','Z5','Z6'})); % 20 combinations of sector*shelf (out of 33)

% Add the missing sector/shelf combinations
% For missing compo take the one of the preceding sector (north) for the same shelf position
SECTOR = AIMS_AVERAGECOMPO(:,1);
SHELF = AIMS_AVERAGECOMPO(:,2);

Missings = [1 1 NaN AIMS_AVERAGECOMPO(SECTOR==3&SHELF==1,4:9) ;...
            1 2 NaN AIMS_AVERAGECOMPO(SECTOR==3&SHELF==2,4:9) ;...
            1 3 NaN AIMS_AVERAGECOMPO(SECTOR==3&SHELF==3,4:9) ;...
            2 1 NaN AIMS_AVERAGECOMPO(SECTOR==3&SHELF==1,4:9) ;...
            2 2 NaN AIMS_AVERAGECOMPO(SECTOR==3&SHELF==2,4:9) ;...
            2 3 NaN AIMS_AVERAGECOMPO(SECTOR==3&SHELF==3,4:9) ;...
            7 1 NaN AIMS_AVERAGECOMPO(SECTOR==6&SHELF==1,4:9) ;...
            7 2 NaN AIMS_AVERAGECOMPO(SECTOR==6&SHELF==2,4:9) ;...
            7 3 NaN AIMS_AVERAGECOMPO(SECTOR==6&SHELF==3,4:9) ;...
            9 1 NaN AIMS_AVERAGECOMPO(SECTOR==8&SHELF==1,4:9) ;...
            9 3 NaN AIMS_AVERAGECOMPO(SECTOR==8&SHELF==3,4:9) ;...
            10 1 NaN AIMS_AVERAGECOMPO(SECTOR==11&SHELF==1,4:9) ;...
            11 2 NaN AIMS_AVERAGECOMPO(SECTOR==10&SHELF==2,4:9) ];
            
AIMS_AVERAGECOMPO = array2table([AIMS_AVERAGECOMPO ; Missings]) ; % Now 33 -> all possible sector*shelf combinations populated

% Assign to each of the 3806 reefs the average community compo representative of shelf and sector
relative_cover = table2array(outerjoin(GBR_REEFS,AIMS_AVERAGECOMPO, 'LeftKeys', [10 11],'RightKeys', [2 1], 'LeftVariables',[1 10 11], 'RightVariables',[4:9]));

[~,I]=sort(relative_cover(:,1));
relative_cover = relative_cover(I,:); % sort by increasing reef ID

% Then swap for the surveyed reefs (n=123) which compo was actually recorded (averaged over the selected period (2006-2008)
AIMS_AVERAGECOMPO_REEF = table2array(varfun(omean, LTMP_Transects_compos,'GroupingVariables',{'KarloID'},...
    'InputVariables',{'Z1','Z2','Z3','Z4','Z5','Z6'}));

relative_cover(AIMS_AVERAGECOMPO_REEF(:,1),4:end) = AIMS_AVERAGECOMPO_REEF(:,3:end); % including 123 reefs with actual community composition recorded

%% 2) Assign to each of the 3806 reefs the average total cover representative of shelf and sector between 2006-2008
StartDate = 2006 ; % earliest date used to initialize coral cover
EndDate = 2009 ; % exclusive

% Calculate average total cover of reefs surveyed with transects %%%%%%%%%%%%%%%%%%%
select_TR = find(LTMP_Transects_means.YEAR_Reefmod >= StartDate & LTMP_Transects_means.YEAR_Reefmod < EndDate);
AIMS_TR = LTMP_Transects_means(select_TR,[16 17 3 9 5]); %5: Years, 10: cover tot, 9: YEAR_reefmod
AIMS_TR.Properties.VariableNames = {'Reef_ID','SECTOR','SHELF','CCOVER','YEAR'};

AIMS_TR.AIMS_SHELF(AIMS_TR.SHELF=='I')=1;
AIMS_TR.AIMS_SHELF(AIMS_TR.SHELF=='M')=2;
AIMS_TR.AIMS_SHELF(AIMS_TR.SHELF=='O')=3;
AIMS_TR = AIMS_TR(:,[1 2 6 4 5]); 
NB_REEFS_TR_INIT = length(unique(AIMS_TR.Reef_ID)); %123 reefs surveyed using transects between 2006-2008

% Select manta tow surveys %%%%%%%%%%%%%%%
load('AIMS_MantaTow_1985_2020.mat') % with SHELF as defined by AIMS
select_MT = find(MantaTow_1985_2020.YearReefMod >= StartDate & MantaTow_1985_2020.YearReefMod < EndDate);
AIMS_MT = MantaTow_1985_2020(select_MT,[16 3 10 8]); % select AIMS shelf classification
AIMS_MT.Properties.VariableNames = {'Reef_ID','SHELF','MEAN_LIVE_CORAL','YEAR'};

AIMS_MT.AIMS_SHELF(char(AIMS_MT.SHELF)=='I')=1;
AIMS_MT.AIMS_SHELF(char(AIMS_MT.SHELF)=='M')=2;
AIMS_MT.AIMS_SHELF(char(AIMS_MT.SHELF)=='O')=3;
NB_REEFS_MT_INIT = length(unique(AIMS_MT.Reef_ID)); %140 reefs surveyed using transects between 2006-2008
NB_REEFS_TOT_INIT = length(unique([unique(AIMS_TR.Reef_ID) ; unique(AIMS_MT.Reef_ID)])); % 186 reefs in total used for intitialisation of coral cover

% convert Manta tows into transect equivalent %%%%%%%%%%%%%
load('LTMP_Transect2Tow_NEW.mat') % (linear model allowing to predict transect coral cover from manta-tow coral cover)
AIMS_MT.MEAN_LIVE_CORAL_EQ = predict(LTMP_Transect2Tow_Model,AIMS_MT.MEAN_LIVE_CORAL);
AIMS_MT.SECTOR = GBR_REEFS.AIMS_sector(AIMS_MT.Reef_ID);
AIMS_MT = AIMS_MT(:,[1 7 5 6 4]);
AIMS_MT.Properties.VariableNames = {'Reef_ID','SECTOR','AIMS_SHELF','CCOVER','YEAR'};

% Bind transect and tow data
AIMS_ALLtmp = [AIMS_MT ; AIMS_TR];
AIMS_ALL = varfun(omean, AIMS_ALLtmp,'GroupingVariables',...
    {'Reef_ID';'SECTOR';'AIMS_SHELF'}, 'InputVariables','CCOVER'); %average across sites/method

AIMS_AVERAGESTATE = table2array(varfun(omean, AIMS_ALL,'GroupingVariables',...
    {'SECTOR','AIMS_SHELF'}, 'InputVariables',{'Fun_CCOVER'})); % 22 combinations of sector*shelf (out of 33)
% Exclude sectors represented by less than 3 reefs
AIMS_AVERAGESTATE = AIMS_AVERAGESTATE(AIMS_AVERAGESTATE(:,3)>=3,:); % 21 combinations of sector*shelf (out of 33)

% Add the missing sector/shelf combinations
SECTOR = AIMS_AVERAGESTATE(:,1);
SHELF = AIMS_AVERAGESTATE(:,2);

Missings = [1 1 NaN AIMS_AVERAGESTATE(SECTOR==3&SHELF==1,4) ;...
    1 2 NaN AIMS_AVERAGESTATE(SECTOR==3&SHELF==2,4) ;...
    1 3 NaN AIMS_AVERAGESTATE(SECTOR==3&SHELF==3,4) ;...
    2 1 NaN AIMS_AVERAGESTATE(SECTOR==3&SHELF==1,4) ;...
    2 2 NaN AIMS_AVERAGESTATE(SECTOR==3&SHELF==2,4) ;...
    2 3 NaN AIMS_AVERAGESTATE(SECTOR==3&SHELF==3,4) ;...
    7 1 NaN AIMS_AVERAGESTATE(SECTOR==6&SHELF==1,4) ;...
    7 3 NaN AIMS_AVERAGESTATE(SECTOR==6&SHELF==3,4) ;...
    9 1 NaN AIMS_AVERAGESTATE(SECTOR==8&SHELF==1,4) ;...
    9 3 NaN AIMS_AVERAGESTATE(SECTOR==8&SHELF==3,4) ;...
    10 1 NaN AIMS_AVERAGESTATE(SECTOR==8&SHELF==1,4) ;...
    11 2 NaN AIMS_AVERAGESTATE(SECTOR==9&SHELF==2,4) ] ;
        
AIMS_AVERAGESTATE = array2table([AIMS_AVERAGESTATE ; Missings],'VariableNames', {'Sector';'Shelf';'N';'TOT'} ) ;

init_coral = table2array(outerjoin(GBR_REEFS,AIMS_AVERAGESTATE, 'LeftKeys', [10 11],'RightKeys', [2 1], 'LeftVariables',[1 10 11], 'RightVariables',[4]));
init_coral = init_coral(I,:); % sort by increasing reef ID (just in case)

clear Z TOT_COVER survey_year select_nan omean SECTOR SHELF Missings AIMS_MTtmp

% Then for each surveyed reef swap the assigned average by the actual value of total cover
% This gives x reefs initialised with the total cover observed, and y reefs initialised with a regional average
init_coral(AIMS_ALL.Reef_ID,4) = AIMS_ALL.Fun_CCOVER;

%% 3) Now populate each reef with a random cover (normally distributed around specified mean)
X = init_coral(META.reef_ID,4);
Y = relative_cover(META.reef_ID,4:end);
varSDcoral = 0.2*ones(1,META.nb_coral_types);

init_coral_cover = X(:,ones(size(Y,2),1)).* Y ;

for s=1:6
    init_coral_cover(:,s) = normrnd(init_coral_cover(:,s), varSDcoral(1,s)*init_coral_cover(:,s))/100;
end

init_coral_cover(init_coral_cover<0.005)=0.005; % minimum 0.5% cover for each group
TCmax = 0.8; %maximum initial total cover is 80%
select_TCmax = find(sum(init_coral_cover,2)>TCmax);
adjust_TCmax = sum(init_coral_cover(select_TCmax,:),2);
init_coral_cover(select_TCmax,:)=TCmax*init_coral_cover(select_TCmax,:)./adjust_TCmax(:,ones(size(init_coral_cover,2),1));

%% Now populate random rubble and sand covers
init_rubble = 0.1; % 10% on average 
init_sand = 0.3; % 30% sand on average

varSDother = 0.2;

init_rubble_cover = normrnd(init_rubble*100, varSDother*init_rubble*100, length(META.reef_ID), 1)/100 ;
init_rubble_cover(init_rubble_cover<0.01) = 0.01;

% init_sand_cover = normrnd(init_sand*100, varSDother*init_sand*100, length(META.reef_ID), 1)/100 ;
% init_sand_cover(init_sand_cover<0.05) = 0.05;

% Update Nov 2022: now sand determined by Roelfsema et al (2021)'s geomorphic and benthic maps
init_sand_cover = normrnd(GBR_REEFS.UNGRAZABLE(META.reef_ID)*100, varSDother*GBR_REEFS.UNGRAZABLE(META.reef_ID)*100)/100 ;
init_sand_cover(init_sand_cover<0.05) = 0.05;


% Adjust if too high (coral + sand > 95%). We leave 5% free space for a safe intialisation
CHECK = sum(init_coral_cover,2) + init_sand_cover;
init_sand_cover(CHECK>0.95) = 0.95 - sum(init_coral_cover(CHECK>0.95,:),2);

clear init_coral init_sand init_rubble relative_cover X Y varSD adjust_TCmax select_TCmax TCmax varSDother CHECK
clear AIMS_AVERAGECOMPO AIMS_AVERAGECOMPO_REEF AIMS_AVERAGESTATE AIMS_MANTA_TOW AIMS_TRANSECT LTMP_Transects_compos LTMP_Transects_means LTMPcoral
clear select_LTMP varSDcoral

%% COTS densities - INPUT DENSITIES MUST BE PER GRID (400m2), NOT PER TOW
REEF_COTS.densities_M = nan(META.nb_reefs,META.nb_time_steps+1);
REEF_COTS.densities_SD = nan(META.nb_reefs,META.nb_time_steps+1);

% First, initialize CoTS populations with CoCoNet hindcast densities of adult CoTS (mean and SD per grid) provided by Scott
load('GBR_past_COTS_CoCoNet.mat')
% Matrices of mean density and SD already converted into nb CoTS per grid for the period 1985-2017.
% Values are mean based on 50 replicate runs, interpolated from 3638 reefs - see 'REEF_COTS_INTERPOLATION_1985_2017.m'
% So 33 years * 2 seasons = 66 columns
CoCoNet_years = 1985:0.5:2017.5;

start_year = find(CoCoNet_years==2007.5) ; %first column is 1985, so 2007 is 1+22
REEF_COTS.densities_M(:,1) = PAST_COTS_M_COCONET(META.reef_ID,start_year); % data are already converted into CoTS densities per grid
REEF_COTS.densities_SD(:,1) = PAST_COTS_SD_COCONET(META.reef_ID,start_year); % data are already converted into CoTS densities per grid

% Then we force CoTS density to past observations where available
load('GBR_PAST_COTS_1992_2020.mat') % Mean coTS per tow between 1992 and 2000 (AIMS LTMP and FMP)
% This is the average number of CoTS per manta tow (or RHIS survey, assuming equivalence with manta tow)
% (Missing values are now identified as NaN, not -1 as previously)

GBR_PAST_COTS_PER_GRID = (GBR_PAST_COTS_NEW(META.reef_ID,:)/0.22)*(1500/2500); % Conversion: 1500 COTS per km2 ~0.22 per manta tow (Moran and De'ath 1986)
% which is further divided by 2500 to get number of COTS per 400m2 (reef grid size)

COTS_years = 1992:0.5:2020;
start_step = find(COTS_years==2007.5);
end_step = size(GBR_PAST_COTS_PER_GRID,2);
length_history = length(start_step:end_step);

RECENT_COTS = max(GBR_PAST_COTS_PER_GRID(:,(start_step-3):start_step),[],2); % maximum density observed from 2006 to 2007.5
tmp2 = find(isnan(RECENT_COTS)==0);

% Replace CoCoNet predictions with recent observations (2006-2007.5)
REEF_COTS.densities_M(find(isnan(RECENT_COTS)==0),1)= RECENT_COTS(find(isnan(RECENT_COTS)==0));
REEF_COTS.densities_SD(find(isnan(RECENT_COTS)==0),1)= 0.28 * RECENT_COTS(find(isnan(RECENT_COTS)==0)); 
% (with 0.28 = slope of linear model between SD and mean CoCoNet predictions)

% Populate from 2008 onwards with observations
REEF_COTS.densities_M(:,2:length_history) = GBR_PAST_COTS_PER_GRID(:,(start_step+1):end_step) ;
REEF_COTS.densities_SD(:,2:length_history) = 0.28 * GBR_PAST_COTS_PER_GRID(:,(start_step+1):end_step) ;


clear GBR_PAST_COTS about_GBR_COTS PAST_COTS_M_COCONET PAST_COTS_SD_COCONET CoCoNet_years end_step
clear GBR_PAST_COTS_NEW GBR_PAST_COTS_NEW_ID GBR_PAST_COTS_PER_GRID about_GBR_PAST_COTS MantaTow_1985_2020

%% SCENARIO OF THERMAL STRESS
if META.doing_genetics == 0 % NO ADAPTATION, NO THERMAL OPTIMUM
    
    DHW = zeros(length(META.reef_ID),META.nb_time_steps);
    
    %% PAST THERMAL STRESS REGIME 2008-2020
    % Using Coral Reef Watch 5km product: max DHW every year from 1985 to 2020 from closest 5x5 km pixel    
%     load('GBR_past_DHW_CRW_5km_1985_2020.mat') %
%     GBR_PAST_DHW = GBR_PAST_DHW(:,24:36); %select 2008 to 2020   
    load('GBR_past_DHW_CRW_5km_1985_2022.mat') % updated NOAA-CRW with 2021 (no bleaching) and 2022 (bleaching)
    GBR_PAST_DHW = GBR_PAST_DHW(:,24:38); %select 2008 to 2022
  
    end_hindcast = size(GBR_PAST_DHW,2);
    DHW(:,1:2:end_hindcast*2) = GBR_PAST_DHW(META.reef_ID,:);
    clear GBR_PAST_DHW
    
    %% FUTURE THERMAL STRESS (CMIP5)
%     if META.nb_time_steps>30 % now starts in 2023
%         
%         % Load the selected forecast scenario of DHW
%         load(['GBR_maxDHW_' char(OPTIONS.GCM) '_rcp' char(OPTIONS.RCP) '_2021_2099.mat'])
%         
%         if (META.nb_time_steps-30)/2 <= size(max_annual_DHW,2)
%             DHW(:,31:2:end) = max_annual_DHW(META.reef_ID,2:(META.nb_time_steps-28)/2);
%         else
%             error('!!REEFMOD ERROR!! Simulated timeframe inconsistent with the DHW time-series. NB_TIME_STEPS has to be <= 26+158 (length of projected DHM time series is 79 years)')
%         end
%     end
    
    %% FUTURE THERMAL STRESS (CMIP6)
    if META.nb_time_steps > 30 
        % Forecast starts in 2023 (step 31)
        % The DHW matrices start in 2000 (column 6) so year 2023 is column 29       
        
        % Load the selected forecast scenario of DHW
        load([char(OPTIONS.GCM) '_' char(OPTIONS.SSP) '_annual_DHW_max.mat'])
        
        if META.nb_time_steps-30 > 2*size(max_annual_DHW(:,29:end),2)
            
            error('!!REEFMOD ERROR!! Simulated timeframe inconsistent with the DHW forecast. NB_TIME_STEPS has to be <= 30+156')

        else
            
            % Select the specified timeframe in the available forecast
            DHW_FORECAST = max_annual_DHW(META.reef_ID,29:(29+(META.nb_time_steps-32)/2));
            % Let's shuffle available years within each decade
            DHW_FORECAST_shuffled = nan(size(DHW_FORECAST));
            start = 1; % first year of the selected forecast
            remain = size(DHW_FORECAST,2); % number of years still available
            
            while remain > 10
                
                sample = randperm(10); % sample at random within the decade
                DHW_FORECAST_shuffled(:,start:(start-1+length(sample)))= DHW_FORECAST(:,start-1+sample); % assign the shuffled years
                start = start + length(sample);
                remain = remain - length(sample);
            end
            
            % Last years available to be shuffled as well (Less than a
            % decade available)
            sample = randperm(remain);
            DHW_FORECAST_shuffled(:,start:(start-1+length(sample)))= DHW_FORECAST(:,start-1+sample);
            
            % Finally assign to the DHW matrix (only in summers)
            DHW(:,31:2:end) = DHW_FORECAST_shuffled;

        end
    end
    
% For genetic simulations (CMIP5)   
%     % First select the climate model
%     climatemodel = 'CCSM4'
%     dhw_method = 'Hotspot_0'
%     
%     % Load the DHM scenarios for all Topt values
%     % Then select the appropriate matrices depending on whether corals have a thermal optimum
%     META.Topt_offset = 3;
%     Topt_scenario_YM = ['GBR_ALL_DHW_' climatemodel '_' dhw_method '_Topt_offset_' num2str(META.Topt_offset) '_'];
%     % with Topt from 10 to 45 by 0.1°C -> 351 stacks of DHM for 1312 and 95 years
%     
%     switch warming_scenario
%         
%         case '2p6' ; load([Topt_scenario_YM 'RCP26.mat']) % Load DHM matrices for different thermal optima
%         case '4p5' ; load([Topt_scenario_YM 'RCP45.mat'])
%         case '6p0' ; load([Topt_scenario_YM 'RCP60.mat'])
%         case '8p5' ; load([Topt_scenario_YM 'RCP85.mat'])
%             
%     end
%     
%     % Only consider DHM for Topt between 15 and 40°C to speed up the code
%     resize_Topt_list = find(T_opt==15):1:find(T_opt==40);
%     GBR_DHM_Max = DHM_Max(META.reef_ID,:,resize_Topt_list) ; % now only 251 stacks
%     T_opt = T_opt(resize_Topt_list) ;
%     
%     % Now load the MMM
%     load('GBR_MMM_CoRTAD.mat') % Gives the baseline for creating Topt values (Topt + offset = MMM)
%     
%     % We now select the MMM from the Topt stacks assuming it corresponds to
%     % MMM ~ Topt-Topt_offset (eg, 3°C) -> this way ensures comparison between
%     % scenario with and without genetics (ie, same DHM used thermal stress)
%     MMM_equivalent = MMM_CoRTAD(META.reef_ID,1) - META.Topt_offset ;
%     Topt2index = 10*max(T_opt)-length(T_opt) ; % then the index of the stack of a given T_opt value is 10*T_opt - Topt2index
%     
%     select_Topt = 10*round(10*MMM_equivalent)/10 - Topt2index;
%     end_forecast = (META.nb_time_steps/2)-end_hindcast;
%     
%     for n=1:META.nb_reefs
%         DHW(META.reef_ID(n),(1+end_hindcast*2):2:META.nb_time_steps) = ...
%             squeeze(((365/12)/7)*GBR_DHM_Max(META.reef_ID(n),1:end_forecast,select_Topt(n)));
%     end

else
    error('Cannot run adaptation with this version of Reefmod')
end


%% PAST STORM REGIME 2008 to 2020
% Category cyclones (Saffir-Simpson scale) assigned to the spatial prediction of >4m wave height from Puotinen et al (2016).
% Original set of predictions is matrice of 0/1 compiled by Rob for the entire GBR (3806 reefs). 
% Then, occurence of damaging wave is blended with category cyclone estimated from maximum sustained wind and distance to cyclone track
% (with 10% increase to meet US standards) measured along each real cyclone track (data extracted from BoM Database of past cyclone tracks).
load('GBR_cyclones_2008-2020_NEW.mat')

CYCLONE_CAT = zeros(length(META.reef_ID),2*size(GBR_PAST_CYCLONES_NEW,2)) ;
CYCLONE_CAT(:,1:2:end) = GBR_PAST_CYCLONES_NEW(META.reef_ID,:);
clear GBR_PAST_CYCLONES_NEW

%% Update Fed 2022: PATCH FOR 2021 (assuming no severe cyclone)
CYCLONE_CAT = [CYCLONE_CAT zeros(length(META.reef_ID), 2)];
    
%% FUTURE STORM REGIME (from 2022 onwards)
if META.nb_time_steps>28
    
    load('Reef_Cyclone_TimeSeries_Count_Cat.mat')
    
    FUTURE_CYCLONE_CAT =  zeros(length(META.reef_ID),2*size(Cyc_cat,2)) ;
    FUTURE_CYCLONE_CAT(:,1:2:end) = Cyc_cat(GBR_REEFS.Nick_ID(META.reef_ID),:,simul);
    
    CYCLONE_CAT = [CYCLONE_CAT FUTURE_CYCLONE_CAT];
    
end

clear GBR_DHM_Max DHM_Max FUTURE_COTS_M FUTURE_COTS_SD CURRENT_CORAL_COVER CURRENT_COTS_DENSITIES
clear CURRENT_RUBBLE_COVER CURRENT_NONGRAZABLE ALL_CYCLONE_CAT Cyc_cat DHWtmp Model_reef_indices ModelReef_ID ModelReef_ID_label
clear SST_GBRtmp Years_unique yr nb_hindcast_runs rel T_opt
clear nb_time_steps Reef_SST_Mean_Yr total_coral_cover Yr Topt_scenario_YM varSD resize_Topt_list
clear FUTURE_CYCLONE_CAT AIMS_ALL AIMS_ALLtmp AIMS_MT AIMS_TR max_annual_DHW LTMP_Transect2Tow_Model
