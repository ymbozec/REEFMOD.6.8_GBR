% Y.-M. Bozec, MSEL, created Jan 2012.
%
% INITIALISATION OF PARAMETERS.
% 
% This performs automatic generation of parameters, 
% such as loading files, basics calculation.
%_________________________________________________________________________________________

%_________________________________________________________________________________________
%
%       MEMORY OPTIMIZATION
%_________________________________________________________________________________________
CORAL.is_brooder = uint8(CORAL.is_brooder) ;
CORAL.size_model = uint8(CORAL.size_model) ;

%_________________________________________________________________________________________
%
%       DISPLAYING THE SELECTED OPTIONS
%_________________________________________________________________________________________

if META.isverbose == 1
    
    disp('------------------------------------------')
    
    if META.nb_reefs == 1
        disp('---> SINGLE REEF SIMULATIONS');
        
        %     if REEF.doing_sedimentation==0
        %         disp('---> OFF sedimentation');
        %     else
        %         disp('---> ON  sedimentation');
        %
        %         if REEF.is_sedimentation_seasonal==0
        %             disp('---> sedimentation IS NOT seasonal');
        %         else
        %             disp('---> sedimentation IS seasonal');
        %         end
        %     end
        
    else
        disp('---> MULTIPLE REEF SIMULATIONS');
        disp('---> ON  seasonality')
        disp('---> OFF sedimentation');
    end
    
    
    if META.track_populations == 1
        
        disp('---> TRACKING COLONY SIZES');
        
        if META.nb_simul>1
            
            META.nb_simul = 1 ; % forces to only one simulation
            disp(' ########################################################')
            disp(' ## PROGRAM FORCED TO 1 SIMULATION BECAUSE OF TRACKING ##')
            disp(' ########################################################')
            
        end
    end
end
%_________________________________________________________________________________________
%
%       AUTOMATIC GENERATIONS
%_________________________________________________________________________________________

% Maximum age of coral colonies in time steps of 6 months so 30 years but constrained by
% the area of the cell - this is now OLD STUFF
% CORAL.max_age = floor(((META.cell_area_cm2/pi)^(1/2))./CORAL.growth_rate);
% Replaced by maximum size (necessary for initialing populations (f_derivecoralcover)
CORAL.adult_max_size = floor(pi*(CORAL.max_diameter/2).^2) ;
% diameter can't be bigger than a x% of a cell
rel = 1 ;
CORAL.adult_max_size(CORAL.adult_max_size>rel*pi*(META.cell_x_size/2)^2)=floor(rel*pi*(META.cell_x_size/2)^2);
% Number of coral species in the model
META.nb_coral_types = length(CORAL.growth_rate);

% Number of algal species in the model
META.nb_algal_types = length(ALGAL.feeding_prefs);

% Planar area of a cell in cm2
META.cell_area_cm2 = META.cell_x_size * META.cell_y_size ;

% Planar area of the reef (grid) in cm2
META.total_area_cm2 = META.grid_x_count * META.grid_y_count * META.cell_area_cm2 ; 

% Generate a scenario of constant herbivory over time if no parrotfish fishing
if META.doing_parrotfish_fishing == 0
    REEF.herbivory = REEF.herbivory(1,ones(META.nb_time_steps,1));
    ALGAL.herbivory_props = ALGAL.herbivory_props(:,ones(META.nb_time_steps,1));
end

ALGAL.growth_rate=[algal_growth_rate_summer algal_growth_rate_winter];
clear algal_growth_rate_summer algal_growth_rate_winter;


% Create the neighbouring "environment" for every cell (TORUS)
% This identifies the cells which are around every cell, based on the type defined in the function
% f_create_environment (currently 4-cell neighboring for the cellular automaton - but could be
% changed in the future (8-cell neighboring?)
META.environ = f_create_environ (META.grid_x_count, META.grid_x_count) ;

% clear all the temporary variables
clear partial_mortality_inci_int_sed partial_mortality_inci_gra_sed
clear partial_mortality_area_int_sed partial_mortality_area_gra_sed
clear recruit_prob_sed growth_rate_sed
%_________________________________________________________________________________________
%
%       SETTING UP BLEACHING, HURRICANES, CONNECTIVITY, DISEASE
%_________________________________________________________________________________________
  
% if META.doing_bleaching == 1   
%     settings_BLEACHING;
% end
% 
% if META.doing_hurricanes == 1   
%     settings_HURRICANES;
% end

if META.doing_3D==1
    settings_3D;
end

if META.doing_coral_competition==1
    settings_CORAL_COMPETITION;
end

if META.doing_COTS==1
    settings_COTS;   
end


if META.doing_coral_connectivity==1 ||  META.doing_COTS_connectivity==1      
    settings_CONNECTIVITY;
else
    CONNECT = [];
end

if META.doing_water_quality==1    
    settings_WATER_QUALITY;
else
    REEF_POP = [];
end


%_________________________________________________________________________________________
%
%       DISPLAY OPTIONS
%_________________________________________________________________________________________

if META.isverbose == 1
    
    if META.doing_coral_connectivity==0
        disp('---> OFF coral connectivity');
    else
        if META.nb_reefs == 1
            error('REEFMOD error: Define the number of reefs for connectivity (in META.nb_reefs)')
        end
        
        for s=1:META.nb_coral_types
            if size(META.ACROCNS(1).M,1) ~= META.nb_reefs
                error('REEFMOD error: connectivity matrices wrong size (must equal the number of reefs!)');
            end
        end
        disp('---> ON  Coral connectivity');
    end
    
    
    if META.doing_bleaching == 0
        disp('---> OFF  bleaching');
    else
        disp('---> ON  bleaching');
    end
    
    
    if META.doing_hurricanes==0
        disp('---> OFF hurricanes');
    else
        disp('---> ON  hurricanes');
    end
    
    
    if META.doing_3D==0
        disp('---> OFF 3D reef');
    else
        disp('---> ON  3D reef');
    end
    
    
    if META.doing_calcification==0
        disp('---> OFF calcification');
    else
        disp('---> ON  calcification');
    end
    
    if META.doing_coral_competition==0
        disp('---> OFF coral/coral competition');
    else
        disp('---> ON  coral/coral competition');
    end
    
    
    if META.doing_COTS_connectivity == 1
        disp('---> ON COTS connecitivty');
    else
        disp('---> OFF COTS connecitivty');
    end  
    
    if META.doing_water_quality == 1
        disp('---> ON Water quality');
    else
        disp('---> OFF Water quality');
    end  
    
    % Display the size of the cell
    str_x_cell = num2str(META.cell_x_size);
    str_y_cell = num2str(META.cell_y_size);
    str_tot_cell = num2str(META.cell_area_cm2/10000);
    disp( ['---> CELL SIZE: '  str_x_cell 'cm * '  str_y_cell 'cm' ' = ' str_tot_cell ' m2'] );
    
    % Display the size of the grid
    str_x_grid = num2str(META.cell_x_size * META.grid_x_count/100);
    str_y_grid = num2str(META.cell_y_size * META.grid_y_count/100);
    str_tot_grid = num2str(META.total_area_cm2/10000);
    disp( ['---> GRID SIZE: '  str_x_grid 'm * '  str_y_grid 'm' ' = ' str_tot_grid ' m2'] );
    
    clear str_tot_cell str_tot_grid str_x_cell str_x_grid str_y_cell str_y_grid
    
end
