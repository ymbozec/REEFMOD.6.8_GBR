% Generate priority list of reefs for restoration
% (extension of 'PRIORITY_LIST.m' created in 09/2019 for RRAP investment case)
%
% Yves-Marie Bozec, y.bozec@uq.edu.au, 08/2021 
%
%__________________________________________________________________________

function priority_list = f_generate_priority_list_NEW(priority_option, MY_REEFS, CONNECT)
% inputs are set of options specified in settings_RESTORATION.m
% returns the list of reef ID sorted from highest to lowest priority
% NOTE THE ORDER OF APPLICATION OF EACH CRITERION IS DETERMINANT (ie,
% EXAMPLE: priority North/South/Central and smallest reefs first will list first all Northern reefs
% of increasing size, then all Southern reefs of increasing size, etc.
REEF_ID = MY_REEFS.Reef_ID;
ORDERS = nan(length(REEF_ID), 6);

%% First criterion is regional
% 0: GBR-wide (no regional priority)
% otherwise vector of indices 1, 2 and 3 in the desired priority

if priority_option.region == 0 ; ORDERS(:,1) = 1;
    
else
    
    for r = 1:length(priority_option.region)
        
        switch priority_option.region(r)
            
            case 1; select_region = find((MY_REEFS.AIMS_sector==1|MY_REEFS.AIMS_sector==2|MY_REEFS.AIMS_sector==3));
            case 2; select_region = find(MY_REEFS.AIMS_sector==4|MY_REEFS.AIMS_sector==5|MY_REEFS.AIMS_sector==6|MY_REEFS.AIMS_sector==7|MY_REEFS.AIMS_sector==8);
            case 3; select_region = find(MY_REEFS.AIMS_sector==9|MY_REEFS.AIMS_sector==10|MY_REEFS.AIMS_sector==11);               
        end
        
        ORDERS(select_region,1) = r;       
    end
end

%% Second criterion is shelf position
% 0: cross-shelf (no shelf priority)
% otherwise vector of indices 1, 2 and 3 in the desired priority
if priority_option.shelf == 0 ; ORDERS(:,2) = 1;
    
else
    
    for s = 1:length(priority_option.shelf)
        
        select_shelf = find(MY_REEFS.Shelf_position == priority_option.shelf(s));       
        ORDERS(select_shelf,2) = s;     
    end
end


%% Third criterion is likelihood of larval sink as the sum of all inbound link strengths (excluding self supply)
% 0: no priority
% 1: highest potential of larval sink
% 2: lowest potential
if priority_option.link_strength == 0 ; ORDERS(:,3) = 1;
    
else
    
    % Calculate mean connectivity strength
    X_connect = CONNECT(1).ACROPORA;
    
    for y = 2:size(CONNECT,2)
        
        X_connect = X_connect + CONNECT(y).ACROPORA;    
    end
    
    self = 1 + diag(-1 * ones(1,length(REEF_ID))); % connectivity matrix of ones with a diagonal of zeros
    Mean_connect = (X_connect.*self)./size(CONNECT,2); % mean connectivity matrix with diagonal of zeros
    Lik1 = sum(Mean_connect,2); % sum of all inbound link stengths for every reef
    
    switch priority_option.link_strength
        
        case 1; [~,J] = sort(Lik1,'descend'); % highest potential first
        case 2; [~,J] = sort(Lik1,'ascend'); % lowest potential first           
    end
    
%     ORDERS(J,3) = REEF_ID; 
        ORDERS(:,3) = J; 
end

%% Fourth criterion is likelihood of larval sink as the number of all inbound links
% 0: no priority
% 1: highest potential of larval sink
% 2: lowest potential
if priority_option.link_number == 0 ; ORDERS(:,4) = 1;
    
else
    
    % Calculate mean connectivity strength
    X_connect = CONNECT(1).ACROPORA;
    
    for y = 2:size(CONNECT,2)
        
        X_connect = X_connect + CONNECT(y).ACROPORA;    
    end
    
    Mean_connect = X_connect/size(CONNECT,2); % mean connectivity matrix (no need to remove self supply here)
    ID_connect = zeros(size(Mean_connect));
    ID_connect(Mean_connect~=0)=1;
    Lik2 = sum(Mean_connect,2); % sum of all inbound link stengths for every reef
    
    switch priority_option.link_number
        
        case 1; [~,K] = sort(Lik2,'descend'); % highest potential first
        case 2; [~,K] = sort(Lik2,'ascend'); % lowest potential first           
    end
    
%     ORDERS(K,4) = REEF_ID; 
        ORDERS(:,4) = K;

end

%% Fifth criterion is reef area
% 0: no area-based priority
% 1: largest areas first
% 2: smallest reef areas first
    
seed=rng;
if priority_option.reef_area == 0

%     ORDERS(:,5) = randperm(length(REEF_ID)); % just assign random ranks
    ORDERS(:,5) = 1; % just assign random ranks

else
    
    switch priority_option.reef_area
        
        case 1; [~,I] = sort(MY_REEFS.Reference_Area_km2,'descend'); % largest area first
        case 2; [~,I] = sort(MY_REEFS.Reference_Area_km2,'ascend'); % smallest area first
        case 3; [~,I] = sort(MY_REEFS.GeomCH_3D_Area_km2,'descend'); % largest area first
        case 4; [~,I] = sort(MY_REEFS.GeomCH_3D_Area_km2,'ascend'); % smallest area first
    end
    
%     ORDERS(I,5) = REEF_ID;
    ORDERS(:,5) = I;

end

rng(seed)

%% Sixth criterion is distance to focal reef
% Requires the ID of the focal reef
if priority_option.focal_reef == 0

    ORDERS(:,6) = 1;

else
        
    F = find(MY_REEFS.Reef_ID == priority_option.focal_reef);
    DIST = distance(MY_REEFS.LAT(F),MY_REEFS.LON(F),MY_REEFS.LAT,MY_REEFS.LON);
[~,L] = sort(DIST,'ascend'); % largest area first

%     ORDERS(I,6) = REEF_ID;
    ORDERS(:,6) = L;

end

rng(seed)


% [A,priority_list] = sortrows(ORDERS,'ascend');
priority_list = ORDERS(:,6); %shortcut here to only use priority #6 (focal reef)

TEST = MY_REEFS(priority_list,:);

