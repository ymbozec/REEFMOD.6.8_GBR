%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Y.-M. Bozec, MSEL, created Nov 2011.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% function [fecundity_adol, fecundity_adult] = f_estimate_fecundity (coral_cm2, CORAL)
% 
% % Calculate fecundity output for a reef in terms of total number of larvae produced over 
% % the grid for a coral species (will be transformed to larval input for other reefs using
% % the connectivity matrix).
% id1 = zeros(size(coral_cm2)) ;
% id2 = id1 ;
% 
% % Calculate fecundity for adults
% id1(coral_cm2 >= CORAL.adult_size) = 1 ;
% fecundity_adult = sum(sum(2*pi*((sqrt(coral_cm2.*id1/pi)).^2)*216)) ;
% 
% % calculate fecundity for adol
% id2(coral_cm2 >= CORAL.adol_size) = 1 ;
% id2 = id2 - id1;
% fecundity_adol = sum(sum(0.25*2*pi*((sqrt(coral_cm2.*id2/pi)).^2)*216)) ;

function [fecundity] = f_estimate_fecundity (coral_cm2, F_list, fecund_min_size, a, b)

% Select colony sizes with 100% gravid
I=find(coral_cm2>= fecund_min_size);
adult_coral_sizes = coral_cm2(I);
adult_relfitness = F_list(I);

% Allometric relationship based on Hall and Hughes 1996:
% Egg volume = exp(a+b*log(size))

all_egg_volumes = exp(a + b*log(adult_relfitness.*adult_coral_sizes)) ; % mm3 of eggs produced by each colony

fecundity = floor(sum(sum(all_egg_volumes))/0.1) ; %0.1 mm3 is the average volume of an egg

