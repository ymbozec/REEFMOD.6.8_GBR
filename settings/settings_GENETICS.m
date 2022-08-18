% Y.-M. Bozec, MSEL, created June 2018.
% Last modified: 07/06/2018
%
% SETTINGS FOR GENETICALLY DRIVEN FITNESS
%_________________________________________________________________________________________

%% THERMAL TOLERANCE FOLLOWING MATZ ET AL (2018)
META.genetics.nb_loci = 10 ; % number of QTLs (loci) - will create 2 alleles for each locus
META.genetics_pop_size = 1000 ; % number of representative corals in a population (gene bank)

% Make sure there are as many values as there are coral types targeted by genetics
META.genetics.group = [ 1 1 1 1 1 1 ]; % ID of coral types
META.genetics.QTL_mu = [ 0 0 0 0 0 0 ] ; % mean of QTL effect size (normal distribution)
META.genetics.QTL_sd = [ 0.2 0.2 0.2 0.2 0.2 0.2 ]; % sd of QTL effect size (normal distribution)

META.genetics.esd = 1+[ 0 0 0 0 0 0 ]; % SD of environmental effect on fitness (mean=0), on top of genetic (0 means perfect heritability)
% Misha used two settings: esd=0 and esd=2 (low heritability). Use 0 to
% simulate the most efficient selection (perfect heritability, narrow tolerance)

META.genetics.SIGMA_HOT = 1+[ 0 0 0 0 0 0 ]; % 
META.genetics.SIGMA_COLD = 1+[ 0 0 0 0 0 0 ]; % for the cold side (when temp<Topt)
% Breadth of thermal tolerance (SD of the bell curve defining stabilising selection) - (named "pl" in Matz' model code)
% Misha used pl = 0.5, 1 and 2, corresponding to 86%, 40% and 13% fitness drop when temperature mismatched the phenotypic optimum by 1 degree C
% Here we introduce asymetry to implement a steeper slope on the hot side (use same value for both for a Gaussian curve)

META.genetics.MuteEffect = [ 0.2 0.2 0.2 0.2 0.2 0.2 ]; % SD of mutation effects at QTLs (mean=0)
META.genetics.MuteRate = 1*[ 1e-6 1e-6 1e-6 1e-6 1e-6 1e-6 ] ; % Mutation rate per locus per gamete


%% ADDITIONS
META.initial_push_QTL_mu = 0.0+META.genetics.QTL_mu; % Mean of QTL effect size to add to every QTL at initialisation
% This allows introducing variability across reefs + little push of thermal
% tolerance to account for past thermal changes (including 2017 bleaching
% which is likely to have selected some thermally resistant)
META.initial_push_QTL_sd = 0.1*META.genetics.QTL_sd ; % SD of QTL added (make it narrower than for initial TQLs)

% Proportion of natural mortality due to chronic thermal stress
% META.genetics.thermal_mort = 0.5 ; % eg, if 0.8 then adaptation affects 80% of natural mortality

% Power of the effect of thermal tolerance on bleaching resistance
% META.genetics.bleaching_resistance = 5 ; % eg, if 2 then resistance is relative fitness elevated at the power of 2
% To change sensitivity to bleaching need to aad a parm for steepness of the sigmoid function and/or position
% of the inflexion point (currently at 0.5 rel sensitivity
