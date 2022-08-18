% -------------------------------------------------------------------------
% Y.-M. Bozec, MSEL, created Aug 2015.
% For tracking population size structure
% -------------------------------------------------------------------------

function [count,class] = f_count_sizefreq(cover_cm2, CORAL)
% CALCULATES THE NUMBER OF CORALS OF A GIVEN SPECIES IN A SIZE CLASS

colony_sizes = cover_cm2(:) ;

juveniles = colony_sizes(colony_sizes < CORAL.juv_max_size) ;
adolescents = colony_sizes(colony_sizes < CORAL.size_threshold_pm & colony_sizes >= CORAL.juv_max_size) ;
adults = colony_sizes(colony_sizes >= CORAL.size_threshold_pm) ;

%% IF SIZE IS DEFINED AS CORAL SURFACE (in cm2)
% [count.juv, class.juv] = hist(juveniles,0:CORAL.size_bins(1):CORAL.juv_max_size);
% [count.adol, class.adol] = hist(adolescents,0:CORAL.size_bins(2):CORAL.size_threshold_pm);
% [count.adult, class.adult] = hist(adults,0:CORAL.size_bins(3):max(CORAL.adult_max_size));

%% IF SIZE IS DEFINED AS CORAL DIAMETER (in cm)
D_juveniles = sqrt(juveniles/pi)*2;
D_adolescents = sqrt(adolescents/pi)*2;
D_adults = sqrt(adults/pi)*2;

edges_juv = 0:CORAL.diam_bins(1):CORAL.juv_max_diam; % note this includes 6 mo old recruits (1cm2 corals)
edges_adol = CORAL.juv_max_diam:CORAL.diam_bins(2):CORAL.adol_max_diam;
edges_adult = CORAL.adol_max_diam:CORAL.diam_bins(3):(CORAL.adult_max_diam+CORAL.diam_bins(3));

[count.juv, class.juv] = histcounts(D_juveniles,edges_juv);
[count.adol, class.adol] = histcounts(D_adolescents,edges_adol);
[count.adult, class.adult] = histcounts(D_adults,edges_adult);