%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Y.-M. Bozec, MSEL, created March 2017 (merging of previous f_COTS_predation and f_COTS_dynamics).
%
% Last update (YM): 02/2021 (track coral loss per species)
%
% This first simulates COTS population dynamics (growth, mortality) over one time step (6 month),
% then process to comsumption of corals following population size. COTS might disappear if coral
% cover drops below a given threshold (end of COTS outbreak) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [coral,algal,new_density_COTS,total_coral_loss] = f_apply_COTS_predation(coral, algal, ...
    density_COTS, density_COTS_settlers, COTS_feeding_rates, SI, COTS_background_density, META)

%%% 1) PROCESS COTS DYNAMICS (GROWTH AND MORTALITY)
temp1_density_COTS = zeros(size(COTS_feeding_rates,2),1) ;
temp1_density_COTS(2:end) = squeeze(density_COTS(1:(end-1))) ;% Increment the age of all COTS before mortality
% (note this eradicates the oldest COTS)
temp1_density_COTS(1) = density_COTS_settlers ; % COTS settlers before mortality 

new_density_COTS = (1-META.COTS_mortality').*temp1_density_COTS' ;


%%% 2) ESTIMATE TOTAL MORTALITY OF EACH CORAL SPECIES
% This accounts for feeding preferences and availability of each coral species
nb_coral_species = size(coral,2) ;
coral_species_planar_area_cm2 = zeros(1,nb_coral_species);

% Note that coral surface is planar, while other models have considered the
% actual surface area (Scandoll 1999)
for s=1:nb_coral_species
    id = zeros(size(coral(s).cover_cm2));
    id(coral(s).cover_cm2>0) = 1 ;
    coral_species_planar_area_cm2(1,s) = sum(sum(coral(s).cover_cm2.*id)) ;
%     coral_species_planar_area_cm2(1,s) = sum(sum(coral(s).cover_cm2(coral(s).cover_cm2>0))) ;
end

% coral_species_surface_area_cm2 = coral_species_planar_area_cm2.*SI' ;

coral_species_prop_area = coral_species_planar_area_cm2/sum(coral_species_planar_area_cm2) ;
% coral_species_prop_area = coral_species_surface_area_cm2/sum(coral_species_surface_area_cm2) ;

COTS_consumption = new_density_COTS.*COTS_feeding_rates ;
weighted_feeding_prefs = (META.COTS_feeding_prefs.*coral_species_prop_area')/sum(META.COTS_feeding_prefs.*coral_species_prop_area');

eaten_coral_cm2 = floor(sum(COTS_consumption(ones(1,nb_coral_species),:).*...
    weighted_feeding_prefs(:,ones(1,size(COTS_consumption,2))),2)) ;

coral_loss = zeros(1, nb_coral_species);

%%% 3) PROCESS CONSUMPTION OF CORAL
% Use coral species-specific consumption rates, or is this covered by feeding prefs?
for s = 1:nb_coral_species % eaten_coral_cm2 has one value for each coral species
    
    if eaten_coral_cm2(s)==0
        
        continue
    else
        
        cover_cm2 = coral(s).cover_cm2 ; % temporary storage to speed up code
%         surface_cm2 = cover_cm2*SI(s); % convert into surface area
        
        I = find(cover_cm2>0) ; % spot the living colonies
        which_eaten = zeros(size(I)) ;
        r = uint32(randperm(length(I))) ;
        
        % Calculate the differential from total consumption and the cumulative sum of all colony
        % areas, picked-up in a random manner
%         cumsum_coral = eaten_coral_cm2(s)*ones(length(I),1) - cumsum(surface_cm2(I(r)), 1) ;
        cumsum_coral = eaten_coral_cm2(s)*ones(length(I),1) - cumsum(cover_cm2(I(r)), 1) ;

        % Positive values in cumsum_coral indicate those living coral colonies that will
        % be removed until the the total area of coral is consumed
        
        last_eaten = find(cumsum_coral<0,1,'first'); % the last eaten colony is the first to produce a negative differential
        % (the size of this colony overtakes the required total of areas consummed by COTs).
        % Here we assume this colony is entirely eaten, which makes the total effectively consumed a little bit higher than
        % expected (conservative assumption). As a result, there is no partial mortality due to COTs (they
        % finish the job before moving to the next colony
        
        cumsum_coral(last_eaten)= 1;
%         cumsum_coral(last_eaten)= 0; % Always leave at least one colony after CoTS predation
        
        which_eaten(I(r),1) = sign(cumsum_coral);
        which_eaten(which_eaten<0,1)=0; % excludes negatives (coral cm2 in excess from consumption)

        temp = zeros(size(cover_cm2));
        temp(which_eaten==1) = cover_cm2(which_eaten==1) ;
        algal(1).cover_cm2 = algal(1).cover_cm2 + sum(temp,2) ;
         
%         coral(s).cover_cm2(which_eaten==1) = - cover_cm2(which_eaten==1) ;
        coral_loss(1,s) = sum(coral(s).cover_cm2(which_eaten==1));
        coral(s).cover_cm2(which_eaten==1) = 0 ; % IF NOT KEEPING TRACK OF THE DEADS
        
%         coral(s).cover_cm2 = coral(s).cover_cm2 - temp ; % IF NOT KEEPING TRACK OF THE DEADS
%         coral_species_planar_area_cm2(1,s) = sum(sum(coral(s).cover_cm2(coral(s).cover_cm2>0)));
        coral_species_planar_area_cm2(1,s) = sum(sum(coral(s).cover_cm2(I)));

    end
end

total_coral_loss = coral_loss;
 

% DIE-BACK IF LESS THAN x% CORAL COVER AFTER CONSUMPTION
if sum(new_density_COTS(META.COTS_dieback_classes))>0 && sum(coral_species_planar_area_cm2(1,META.COTS_pref_corals))/META.total_area_cm2 < META.COTS_coral_threshold    

    new_density_COTS(META.COTS_dieback_classes) = COTS_background_density * new_density_COTS(META.COTS_dieback_classes)/sum(new_density_COTS(META.COTS_dieback_classes)) ;
   
end
