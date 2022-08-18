function[META] = f_give_boatsProperties(META, input)

numb = META.COTS_cull_boats ;

for i = 1:numb %here you could remove the loop, if numb is a scalar integer
    %identifyTheBoat
    META.boatProperties.boatID = 1:numb ;
    %alotted boat days per 6 months
    META.boatProperties.boatDays = repmat(input(6),1,numb);%randi([175 250],1,numb);
    %alotted voyages per 6 months
    META.boatProperties.voyages = repmat(input(7),1,numb);%randi([10 12],1,numb);
    %number of divers onboard 
    META.boatProperties.divers = repmat(8,1,numb); %randi([8 12],1,numb); %AMPTO = 10
    %maximum days the vessel can be at sea at one time
    META.boatProperties.maxDays_atSea = repmat(13,1,numb);%randi([20 25],1,numb);
    %the location on the coast that the vessel is stationed
    % 1 = Port Douglas, 2 = Cairns, 3 = Townsville, 4 = Mackay
    %homeports=[2 2 1 1 3 3 4 4 2 2 1 1 3 3 4 4]; % extend this vector if your wish to have more than 8 vessels
    homeports=[2 2 1 1 3 3 4 4 2 2 1 1 3 3 4 4 2 2 1 1 3 3 4 4 2 2 1 1 3 3 4 4 2 2 1 1 3 3 4 4]; % extended for modelling up to 40 vessels
    META.boatProperties.homePort = homeports(1:numb);%randi([1 4],1,numb);
    %Is this vessel a cull boat (1) or a survey vessel (0) currently not active
    META.boatProperties.mission = repmat(1,1,numb);%randi([0 1],1,numb);

    META.boatProperties.effortQuota = zeros(1, length(META.boatProperties.boatID));
    
    %does the boat visit reefs in specific fixed order of their raniking, or does it choose randomly from the top X/reefs2cull reefs where to go first
    %1 if strictly fixed order, 0 if sampling
    META.boatProperties.fixedOrder=repmat(input(11),1,numb);
    
    META = f_cull_effort_K(META);

end