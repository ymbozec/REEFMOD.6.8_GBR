%_________________________________________________________________________________________
%
%       WATER QUALITY SETTINGS
%_________________________________________________________________________________________

load('GBR_REEF_POP.mat') % Baseline scenarios (8 years of demographic layers)

% Only select the reefs to be simulated
for yr = 1:size(GBR_REEF_POP,2) % for every year
    
    REEF_POP(yr).CORAL_recruit_survival = GBR_REEF_POP(yr).CORAL_recruit_survival(META.reef_ID,1); % no need of winter values
    REEF_POP(yr).CORAL_juvenile_growth = GBR_REEF_POP(yr).CORAL_juvenile_growth(META.reef_ID,:);
    REEF_POP(yr).COTS_larvae_survival = GBR_REEF_POP(yr).COTS_larvae_survival(META.reef_ID,1); % no need of winter values
    
    if META.doing_COTS == 1
        % Need a lower bound for survival to allow for the slightest dispersal
        REEF_POP(yr).COTS_larvae_survival(REEF_POP(yr).COTS_larvae_survival<META.COTS_min_larval_survival)=META.COTS_min_larval_survival;
    end
end

% Coral spawning dates only available for 2011-2016
% -> Spawning season 2010-11 and 2017-18 are missing fill with one of the eReefs layers (2011-16) picked up at random
REEF_POP(1).CORAL_larvae_production = GBR_REEF_POP(randi(6)).CORAL_larvae_production(META.reef_ID,1); % no need of winter values
REEF_POP(2).CORAL_larvae_production = GBR_REEF_POP(1).CORAL_larvae_production(META.reef_ID,1); % no need of winter values
REEF_POP(3).CORAL_larvae_production = GBR_REEF_POP(2).CORAL_larvae_production(META.reef_ID,1); % no need of winter values
REEF_POP(4).CORAL_larvae_production = GBR_REEF_POP(3).CORAL_larvae_production(META.reef_ID,1); % no need of winter values
REEF_POP(5).CORAL_larvae_production = GBR_REEF_POP(4).CORAL_larvae_production(META.reef_ID,1); % no need of winter values
REEF_POP(6).CORAL_larvae_production = GBR_REEF_POP(5).CORAL_larvae_production(META.reef_ID,1); % no need of winter values
REEF_POP(7).CORAL_larvae_production = GBR_REEF_POP(6).CORAL_larvae_production(META.reef_ID,1); % no need of winter values
REEF_POP(8).CORAL_larvae_production = GBR_REEF_POP(randi(6)).CORAL_larvae_production(META.reef_ID,1); % no need of winter values