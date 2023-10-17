%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Y.-M. Bozec, MSEL, created Nov 2011.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This function outputs a vector of colony sizes (cm2) for a given coral species 
% c indicates which coral species is in use

function colony_sizes = f_derivecoralcover(c,REEFn,CORAL)

% Calculate the desired total amount of coral cover (cm2) on grid
initial_cover_cm2 = REEFn.initial_coral_cover(c)*sum(REEFn.substrate_SA_cm2);
% Note that sand is included in the substrate area

maximum_size_cm2 = CORAL.adult_max_size(c); % note maximum size of corals is bounded by cell size in INITIALISATION.m
colony_sizes = [];

switch CORAL.size_model(c)
    
    case 1 % Uniform size distribution model (equal probability for every size)
        % NOTE previous ReefMod (<3.1) did play with age instead of size (equal probability of age).
        % Because size is a power function of age and because age had discrete values, this produced more colonies in small than large size classes 
        %(although equal probability to observe the discrete size values). Now we play with sizes, not ages
        s = randi(maximum_size_cm2,1,2000); % Preg-generation of random sizes (save time)
        
    case 2 % Lognormal distribution of colony sizes
        s = lognrnd(CORAL.size_mean(c),CORAL.size_var(c),1,2000) ; % Preg-generation of random sizes (save time)
        
end

s(s > maximum_size_cm2) = []; % bound by max size of species c
s(s < 2) = 2;% bound by arbitrary min size -> avoids zeros, negatives, and a myriad of setllers (recruitment happens right after)

s = repmat(s,1,10); % duplicate s
sum_colonies = cumsum(s) ;
I = find(sum_colonies < (initial_cover_cm2 - maximum_size_cm2));
colony_sizes = s(1:(max(I)+1));

% The last colony matches the desired total cover
colony_sizes=round([colony_sizes, (initial_cover_cm2 - sum(colony_sizes))]) ;
colony_sizes(colony_sizes > maximum_size_cm2) = maximum_size_cm2; % bound by max size of species c
