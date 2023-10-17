
function[META] = f_settings_COTS_control(inputs, META)
%% Parameterize COTS control strategies 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %timestep in 6-month intervals when control should start, e.g. to allow for burning in of the model
    META.COTS_control_start=inputs(16);
    %choose control strategy
    META.COTS_control_strat=inputs(1);
    %choose method of selecting reefs to manage
    META.COTS_reefs2cull_strat=inputs(2);
    %cull surrounding reefs
    META.COTS_cull_surrounding=inputs(3);
    %specify the available effort in number of boats and days per boat; implement surveys later
    META.COTS_cull_boats=inputs(4);
    META.COTS_survey_boats=inputs(5);
    META.COTS_cull_days=inputs(6);%100 per AMPTO
    %number of discrete voyages per boat per 6 months
    META.COTS_cull_voyages=inputs(7);%13 per AMPTO
    META.COTS_survey_voyages=inputs(8);

    [META] = f_give_boatsProperties(META, inputs) ;

    %specify whether full state of the system is known; implement partial knowledge later
    META.COTS_state_known=1;
    %picking the highest priority reef
    META.top_reef_picks=inputs(10);
    META.top_reef_picks=inputs(10);
    
  
    %CS=load('COTS_sites.mat');
    CS=load('COTS_sites_new.mat');%%based on geomorphic habitat maps
    META.COTS_control_sites=CS.COTS_sites;
    
    
Newlist=[META.reef_ID META.COTS_control_sites(META.reef_ID,2)]; %%Add index to track other metrics %% YM correction for sizeable code
META.cntrl_sites=rmmissing(Newlist); %%Remove missing reefs
META.cntrl_reefID=[META.cntrl_sites(:,1)];%% get the reefID for reefs within MPA

    %META.COTS_pref_coral_groups=1:4;%which coral gorups are taken into account for ET
    META.do_CSIRO_COTSctrl=inputs(13);%whether we do CSIRO-like algorithms for COTS control
    META.COTS_global_trigger=inputs(14);%check global trigger for CSIRO COTS control
    META.COTS_reef_trigger=inputs(15);%check reef-level trigger for CSIRO COTS control
    META.COTS_postcontrol_proportions=[ 1 1 (1-META.COTS_detectability(3:end))./(sum(1-META.COTS_detectability(3:end)))];%this is the popualtion structure at ET after control; base don detectability, i.e. survivng population=1-detectability
end
    %ecological threshold above which a reef must be to be treated
    %META.COTS_ecological_threshold=0.22;
 
