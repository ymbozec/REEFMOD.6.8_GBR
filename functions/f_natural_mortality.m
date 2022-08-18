%__________________________________________________________________________
%
% Natural mortality (partial and whole-colony)
%
% Yves-Marie Bozec, y.bozec@uq.edu.au
% Modified in August 2013 after changes in parameterisation as shown in Mumby et al. 2013 appendix
% Modified July 2018 to integrate genetic adaptation to warming
% Modified back to previous version in Aug 2018 as adaptation now affects growth
%__________________________________________________________________________
function [coral_cm2]=f_natural_mortality(coral_cm2, CORAL, REEF, species_ID, nb_coral_types)

%%%% This is new stuff (Aug 2013) for implementing species-specific natural mortalities 
id0 = zeros(size(coral_cm2)) ;
sensitivity_whole_natural = id0;
sensitivity_partial_natural = id0;
extent_partial_natural = id0 ;

id1 = id0+1;
id1(coral_cm2 <= 0) = 0 ; % Remove the already dead ones

col_start = 1;
col_stop = 0;

for s = 1:nb_coral_types
    
    col_stop = col_stop + species_ID(s) ; 
    sensitivity_whole_natural(:,col_start:col_stop)= CORAL.sensitivity_whole_natural(s);
    sensitivity_partial_natural(:,col_start:col_stop)= CORAL.sensitivity_partial_natural(s);
    extent_partial_natural(:,col_start:col_stop)= CORAL.extent_partial_natural(s);
    col_start = col_start + species_ID(s) ;
    
end

%----------------------------------
% WHOLE COLONY MORTALITY
%----------------------------------
id_juv = id1;
id_adult = id1;

id_juv(coral_cm2 >= CORAL.juv_max_size) = 0 ; % remove adol and adults
% id_adult(coral_cm2 < CORAL.adult_size) = 0 ; % remove juv and adol
id_adult(coral_cm2 < CORAL.size_threshold_pm) = 0 ;
% We now use the threshold size used by Edwards et al 2011 (250cm2)
id_adol = id1 - id_juv - id_adult ; % deduce adol ie between CORAL.juv_max_size and CORAL.size_threshold_pm

% rand_mort = rand(size(coral_cm2));
rand_mort = sprand(id1);

rand_mort_juv = rand_mort.*id_juv ;
rand_mort_adol = rand_mort.*id_adol ;
rand_mort_adult = rand_mort.*id_adult ;

% switch to 0 coral stages escaping mortality by chance
id_juv(rand_mort_juv > REEF.juv_whole_mortality_rate) = 0 ; % not species specific at those sizes
id_adol(rand_mort_adol > REEF.adol_whole_mortality_rate*sensitivity_whole_natural) = 0 ;
id_adult(rand_mort_adult > REEF.adult_whole_mortality_rate*sensitivity_whole_natural) = 0 ;

% dead corals are those remaining (id=1)
id_dead = id_adol + id_adult + id_juv ;

%% TEMP - try optimization

% First way  - does not work! need to account for species sensitivity!!
% mort_juv = id_juv;
% mort_juv(mort_juv==1)=binornd(1,REEF.juv_whole_mortality_rate);
% mort_adol = id_adol;
% mort_adol(mort_adol==1)=binornd(1,REEF.adol_whole_mortality_rate*sensitivity_whole_natural);
% mort_adult = id_adult;
% mort_adult(mort_adult==1)=binornd(1,REEF.adult_whole_mortality_rate*sensitivity_whole_natural);
% 
% id_dead = id_adol.*mort_adol + id_adult.*mort_adult + id_juv.*mort_juv ;
% 

% Second way - does not work! need to account for species sensitivity!!
% mort_juv = binornd(1,REEF.juv_whole_mortality_rate,size(id1));
% mort_adol = binornd(1,REEF.adol_whole_mortality_rate,size(id1));
% mort_adult = binornd(1,REEF.adult_whole_mortality_rate,size(id1));

% surviv_juv = 1-binornd(1,REEF.juv_whole_mortality_rate,size(id1));
% surviv_adol = 1-binornd(1,REEF.adol_whole_mortality_rate,size(id1));
% surviv_adult = 1-binornd(1,REEF.adult_whole_mortality_rate,size(id1));

% id_dead = id_adol.*mort_adol + id_adult.*mort_adult + id_juv.*mort_juv ;


% Update coral cover (turn into negative size the dead colonies)
% coral_cm2(id_dead==1) = - coral_cm2(id_dead==1) ;
coral_cm2(id_dead==1) = 0 ; % IF NOT KEEPING TRACK OF THE DEADS
 
%----------------------------------
% PARTIAL COLONY MORTALITY
%----------------------------------
% rand_mort2 = rand(size(coral_cm2));
rand_mort2 = sprand(id1);

% allocate space
area_lost = id0 ;
% log_calc = id0;
% proba = id0;
% 
% log_calc(coral_cm2>0) = log(coral_cm2(coral_cm2>0));
% proba(coral_cm2>0) = 1-(CORAL.partial_mortality_inci_int + ...
%     CORAL.partial_mortality_inci_gra * log_calc(coral_cm2>0))/100 ;
% 
% % Adjust the mortality for each species
% proba = proba.*sensitivity_partial_natural ;
% 
% lost_factor(rand_mort2 < proba) = CORAL.partial_mortality_area_int + ...
%     CORAL.partial_mortality_area_gra * log_calc(rand_mort2 < proba) ;
% 
% lost_factor(lost_factor<0)=0;  % because the formula above generates negatives
% 
% % Since Mumby et al. (2014) we change the calculation of area lost for every affected:
% area_lost(rand_mort2 < proba) = (exp(lost_factor(rand_mort2 < proba))-1)/100 ;

%% Faster code:
log_calc = log(coral_cm2.*id1+0.01);
proba = id1.*(1-(CORAL.partial_mortality_inci_int + ...
    CORAL.partial_mortality_inci_gra * log_calc)/100) ;

% Adjust the mortality for each species
proba = proba.*sensitivity_partial_natural ;

area_lost(rand_mort2 < proba) = (exp(CORAL.partial_mortality_area_int)*...
    (coral_cm2(rand_mort2 < proba).^CORAL.partial_mortality_area_gra)-1)/100 ;

area_lost(area_lost<0)=0;

%%
area_lost = round(extent_partial_natural.*area_lost);

% cannot loose more than the colony area
area_lost(area_lost > coral_cm2 & coral_cm2>0) = coral_cm2(area_lost > coral_cm2 & coral_cm2>0) ;

%%%% Update the grid
coral_cm2 = coral_cm2 - area_lost;  % remove lost tissues
% algal_cm2(:,1) = algal_cm2(:,1) + sum(area_lost,2) ; % update turf with coral loss
