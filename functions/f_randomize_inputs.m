%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Y.-M. Bozec, MSEL, created Aug 2012.
% Last modified: 28/08/2012
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function randomizes input parameters for Reefmod simulations.
% It takes the X from PARAMETERS_DEFAULT or MAIN (if updated)
% and returns a random number Y from a normal distribution defined by the mean assumed to be X 
% and a standard deviation SD which is defined hereafter
% Note that random values are generated for all simulations here rather
% than doing the generation each time a simulation is launched to reduce
% computation time. Values are stored in RANDOM_INPUTS which is a structure
% array
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function RANDOM_INPUTS = f_randomize_inputs(META, REEF, CORAL)

%%%%% Randomize initial coral cover %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
X_init_cover = REEF.initial_coral_cover(:,ones(1,META.nb_simul)) ;
SD2_init_cover = REEF.initial_coral_cover_SD2(:,ones(1,META.nb_simul)) ;

Y_init_cover = transpose(normrnd(X_init_cover,SD2_init_cover)); %creates a matrix of initial cover
% that will be read row after row (1 row for 1 simulation)

check_negatives = 0*Y_init_cover ;
check_negatives(Y_init_cover < 0) = 1 ;

if sum(sum(check_negatives)) >= 1
    
    Y_init_cover(check_negatives==1) = 0.005 ; % turns to 0.5% cover the zeros and the negatives
    
    disp(' ')
    disp('****** [REEFMOD WARNING] *******************************');
    disp('The randomization process has generated negative values')
    disp('of initial cover for at least one species.');
    disp('You may consider reducing the range of possible values')
    disp('for some coral species (ie, lowering standard deviation)');
    disp(' ')
    disp('Reefmod is now running with 0 instead of negatives')
    disp('********************************************************');
    disp(' ')
    
end

RANDOM_INPUTS.initial_cover = Y_init_cover ;

%%%%% Randomize coral growth rate %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
X_growth_rate = CORAL.growth_rate(:,ones(1,META.nb_simul)) ;
SD2_growth_rate = CORAL.growth_rate_SD2(:,ones(1,META.nb_simul)) ;

Y_growth_rate = transpose(normrnd(X_growth_rate,SD2_growth_rate));

check_zeros = 0*Y_growth_rate ;
check_zeros(Y_growth_rate <= 0) = 1 ;

if sum(sum(check_zeros)) >= 1
    
    disp(' ')
    disp('****** [REEFMOD ERROR] *************************************');
    disp('The randomization process has generated negative/zero values')
    disp('of coral growth rate for at least one species.');
    disp('Please consider reducing the range of possible values')
    disp('for some coral species (ie, lowering standard deviation)');
    disp(' ')
    disp('Reefmod cannot run with null growth rates')
    disp('************************************************************');
    disp(' ')
    
    error('>>>>> REEFMOD ERROR >>>>>')
    
end

RANDOM_INPUTS.growth_rate = Y_growth_rate ;


%%%%% Randomize ungrazable substratum %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
X_nongrazable = REEF.nongrazable_substratum(:,ones(1,META.nb_simul))  ;
SD2_nongrazable = REEF.nongrazable_substratum_SD2(:,ones(1,META.nb_simul))  ;

Y_nongrazable = transpose(normrnd(X_nongrazable,SD2_nongrazable));
Y_nongrazable = round(Y_nongrazable*100)/100; % because the proportion of
% nongrazable cannot have more than 2 decimals
Y_nongrazable(Y_nongrazable < 0) = 0 ; % cannot be negative

RANDOM_INPUTS.nongrazable_substratum = Y_nongrazable ;


%%%%% Randomize Herbivory %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
X_herbivory = REEF.herbivory(:,ones(1,META.nb_simul))  ;
SD2_herbivory = REEF.herbivory_SD2(:,ones(1,META.nb_simul)) ;

Y_herbivory = transpose(normrnd(X_herbivory,SD2_herbivory));
% Y_nongrazable = round(Y_nongrazable*100)/100; % because the proportion of
% nongrazable cannot have more than 2 decimals
Y_herbivory(Y_herbivory < 0) = 0 ; % cannot be negative

RANDOM_INPUTS.herbivory = Y_herbivory ;

%%%%% Randomize initial rugosity of the substrate %%%%%%%%%%%%%%%%%%%%%%%%%
X_rugosity = REEF.initial_rugosity(:,ones(1,META.nb_simul))  ;
SD2_rugosity = REEF.initial_rugosity_SD2(:,ones(1,META.nb_simul))  ;

Y_rugosity = transpose(normrnd(X_rugosity,SD2_rugosity));

Y_rugosity(Y_rugosity < 1) = 1 ; % cannot be below 1

RANDOM_INPUTS.initial_rugosity = Y_rugosity ;

% total_cover_init = sum(RANDOM_INPUTS.initial_cover,2) ;
% test = total_cover_init./RANDOM_INPUTS.initial_rugosity;
% RANDOM_INPUTS.initial_rugosity(test>0.25)= total_cover_init(test>0.25)/0.25;


%%%%% Randomize Fish bioerosion %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
X_fish_bioerosion = REEF.fish_bioerosion(:,ones(1,META.nb_simul))  ;
SD2_fish_bioerosion = REEF.fish_bioerosion_SD2(:,ones(1,META.nb_simul)) ;

Y_fish_bioerosion = transpose(normrnd(X_fish_bioerosion,SD2_fish_bioerosion));
% Y_nongrazable = round(Y_nongrazable*100)/100; % because the proportion of
% nongrazable cannot have more than 2 decimals
Y_fish_bioerosion(Y_fish_bioerosion < 0) = 0 ; % cannot be negative

RANDOM_INPUTS.fish_bioerosion = Y_fish_bioerosion ; 
