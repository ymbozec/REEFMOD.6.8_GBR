% Compute grazing and macroalgal growth following Bozec et al. (2019)
%
% Yves-Marie Bozec, y.bozec@uq.edu.au, 11/2016
%__________________________________________________________________________

function [algal_cm2,last_surface_area_grazed] = f_algal_dynamics(algal_removal,ALGAL,total_coral_cm2, algal_cm2, ...
    coral_reduce_macrogrowth, growth_rate, cell_area_cm2, convert_to_canopy, season, dictyota_declines_seasonally, grazable_cell,substrate_SA_cm2,environ)

% Just need the growth rates of Dict and Lob
r_DICT = growth_rate(2,season+1) ;  % to be converted in cm2?
r_LOB = growth_rate(3,season+1) ;  % to be converted in cm2?
 
% List of grazable cells (exclude sand for cell loops)
list_cell = 1:(size(algal_cm2,1)) ;
list_cell = list_cell(grazable_cell==1);

substrate_env_cm2= substrate_SA_cm2(environ(:,2)) + substrate_SA_cm2(environ(:,3)) + ...
    substrate_SA_cm2(environ(:,4)) + substrate_SA_cm2(environ(:,5)) ;

% Simulate algal dynamics
for t=1:ALGAL.nb_step_algal_dynamics
    
    actual_algal_consumpt_pct = f_algal_removal(algal_cm2, algal_removal, ALGAL.feeding_prefs, sum(cell_area_cm2)) ;
    % actual_algal_consumpt_pct gives the total amount (in %) of each alga to be removed based on their availability
    % Note that REEF.substrate_SA_cm2 is the surface area of the substrate underneath live corals
    % -> we use now the actual surface area of the reef instead of the planar area (META.total_area_cm2)

    % Converts in total amount of cm2 algae that can be consumed
    max_fish_consump = actual_algal_consumpt_pct*sum(cell_area_cm2);
    
    %%%%%%%%%%%%%%%%%% Set up grazing for every cell %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    r = uint16(randperm(length(list_cell)));  % randomize the order of cell to visit
    cumulated_algae = cumsum(algal_cm2(list_cell(r),:), 1);
    cumulated_algae(end,1)=cumulated_algae(end,1)-1; % fix to leftover during full grazing
    cumsum_algal = max_fish_consump(ones(size(list_cell(r),2),1),:) - cumulated_algae ;    
    eaten_algal = 0*algal_cm2;
    eaten_algal(list_cell(r),:) = sign(cumsum_algal);
    eaten_algal(algal_cm2==0)=0;
    eaten_algal(eaten_algal<0)=0; % just remove the negatives% GRAZABLE=REEF.grazable_SA_cm2(1:3,1)
    % NOTE in the last selected cell the alga is completely removed, not just the required cover.   
    
    %%%%%%%%%%%%%%%%%% Remove LOBOPHORA due to grazing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % amount of dict sitting above lob already. Note max(x, 0) gives a 0 when x is negative
    existing_colocation_cm2 = max(sum(algal_cm2,2) + total_coral_cm2 - cell_area_cm2,0); % max(x, 0) gives a 0 when x is negative
    % where lobophora is eaten so is the dictyota that is with it
    algal_cm2(eaten_algal(:,3)==1,2) = max(algal_cm2(eaten_algal(:,3)==1,2) - existing_colocation_cm2(eaten_algal(:,3)==1),0) ;
    % all lobophora eaten so it turns to EAM
    algal_cm2(eaten_algal(:,3)==1,1) = algal_cm2(eaten_algal(:,3)==1,1) + algal_cm2(eaten_algal(:,3)==1,3) ;
    % no more lobophora
    algal_cm2(eaten_algal(:,3)==1,3)= 0 ;
    
    %%%%%%%%%%%%%%%%%% Now process DICTYOTA grazing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Update the amount of dictyota which co-located with lobophora already
    existing_colocation_cm2 = max(sum(algal_cm2,2) + total_coral_cm2 - cell_area_cm2,0); % max(x, 0) gives a 0 when x is negative
    
    % check existing_colocation_cm2 with lob_over_dict_cm2 t see how big the difference is
    if (dictyota_declines_seasonally == 1 && season == 1)  %season = 1|0 -> winter|summer
        % if winter, all dictyota killed so it turns to turf
        algal_cm2(:,1) = algal_cm2(:,1) + algal_cm2(:,2) - existing_colocation_cm2 ;
        % no more dictyota, any that was with lobophora is just lobophora now% GRAZABLE=REEF.grazable_SA_cm2(1:3,1)
        algal_cm2(:,2) = 0 ;
        
    else % otherwise all will be eaten "including dislodging any lob underneath"
        % all dictyota killed so it turns to turf
        algal_cm2(eaten_algal(:,2)==1,1) = algal_cm2(eaten_algal(:,2)==1,1) + algal_cm2(eaten_algal(:,2)==1,2) ;
        % also reduce any co-located lobophora
        algal_cm2(eaten_algal(:,2)==1,3) = max(algal_cm2(eaten_algal(:,2)==1,3) - existing_colocation_cm2(eaten_algal(:,2)==1),0) ;
        % no more dictyota
        algal_cm2(eaten_algal(:,2)==1,2) = 0 ;
    end
    
    % Then eat on THICK TURF
    algal_cm2(eaten_algal(:,4)==1,1) = algal_cm2(eaten_algal(:,4)==1,1) + algal_cm2(eaten_algal(:,4)==1,4) ;
    algal_cm2(eaten_algal(:,4)==1,4) = 0 ;
    
    %%%%%%%%%%%%%%%%%% Now proceed to algal growth %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Determine environment for lob growth (grow faster in dense vegetation)
    % Lob covers in the 4 cells surrounding every cell
    lob_env_cm2 = algal_cm2(environ(:,2),3) + algal_cm2(environ(:,3),3) + ...
        algal_cm2(environ(:,4),3) + algal_cm2(environ(:,5),3) ;
        
    lob_env_prop = lob_env_cm2./substrate_env_cm2;% Proportion (prop) of lob cover over 4 cells
    
    % Overgrowth from the margins (surrounding cells) follwoing de Ruyter van Steveninck and Breeman (1987)
    % Filtered later so that only ungrazed cells are overgrown
    lob_overgrowth_cm2 =  ALGAL.margins_lob_overgrowth_cm2*lob_env_prop/0.6 ;
    
    % Adjust for available space
    space_avail_cm2 = cell_area_cm2 - total_coral_cm2 - algal_cm2(:,3) - algal_cm2(:,2);
    lob_overgrowth_cm2(lob_overgrowth_cm2 > space_avail_cm2) = space_avail_cm2(lob_overgrowth_cm2 > space_avail_cm2);
    
    init_LOB_cm2 = algal_cm2(:,3);
    init_DICT_cm2 = algal_cm2(:,2)/convert_to_canopy;
    init_TURF_cm2 = algal_cm2(:,4);
    
    % Lob and Dict can grow only in cells that were not grazed, for any kind of algae
    id_nongrazed = 0*algal_cm2(:,1) ;
    id_nongrazed(sum(eaten_algal(:,[1 4]),2)==0) = 1 ; % id of cells that were not grazed
    id_nongrazed = id_nongrazed.*grazable_cell;  % exclude non-grazable cells
    
    % Set up thick turf: expands quickly in the absence of grazing
    init_TURF_cm2(id_nongrazed==1) = cell_area_cm2(id_nongrazed==1) - total_coral_cm2(id_nongrazed==1)...
        - init_LOB_cm2(id_nongrazed==1) - init_DICT_cm2(id_nongrazed==1) ;
    
    % Set up algal settlement within TURF (use the same settlement rate for UMA; scale settlement down to available turf
    seeding = 0*init_TURF_cm2;
    seeding(id_nongrazed==1) = 2*ALGAL.settlement_lob_cm2 * init_TURF_cm2(id_nongrazed==1)./cell_area_cm2(id_nongrazed==1) ;
    
    % Macroalgae can grow only in ungrazed cells
    id_growth=find(id_nongrazed==1 & init_TURF_cm2>0); % can grow only if ungrazed AND some turf is available 
    
    if isempty(id_growth) == 0 % if there are cells where algae can grow
        
        % Seeding with recruits
        dict_cm2 = init_DICT_cm2(id_growth) + seeding(id_growth)/2;
        lob_cm2 = init_LOB_cm2(id_growth) + seeding(id_growth)/2;
        turf_cm2 = init_TURF_cm2(id_growth) - seeding(id_growth) ;
        
        coral_effect = coral_reduce_macrogrowth(id_growth);
        
        % Lob grows first
        K_LOB =  turf_cm2 + lob_cm2 ;
        DELTA_LOB = r_LOB * lob_cm2.*(1 - (lob_cm2./K_LOB)).*(1 - coral_effect);
        lob_cm2 = floor(lob_cm2 + DELTA_LOB + lob_overgrowth_cm2(id_growth)) ;
        % (add lob_overgrowth after growth)
        
        % Re-calculate turf to update Dict carrying capacity
        turf_cm2 = ceil(turf_cm2 - DELTA_LOB - lob_overgrowth_cm2(id_growth)) ;
        
        % re-adjust where lob_overgrowth exceeds available TURF (ie, where TURF_cm2<0)
        lob_cm2(turf_cm2<0) = lob_cm2(turf_cm2<0) + turf_cm2(turf_cm2<0);
        turf_cm2(turf_cm2<0) = 0;
        
        % Dict grows after Lob
        K_DICT = turf_cm2 + dict_cm2 ;
        DELTA_DICT = r_DICT * dict_cm2.*(1 - (dict_cm2./K_DICT)).*(1 - coral_effect);
        dict_cm2 = floor(dict_cm2 + DELTA_DICT) ;
        
        % Update algal covers
        algal_cm2(id_nongrazed==1,1)=0 ; % Turn EAM into 0 in the non-grazed cells
        algal_cm2(id_growth,3) = lob_cm2;
        algal_cm2(id_growth,2) = dict_cm2;
        algal_cm2(id_nongrazed==1,4) = ceil(cell_area_cm2(id_nongrazed==1) - total_coral_cm2(id_nongrazed==1)...
            - algal_cm2(id_nongrazed==1,2) - algal_cm2(id_nongrazed==1,3)) ;
             
        % Finally estimate true cover of DICT considering it can overtop LOB (but easier to process this here)
        algal_cm2(id_growth,2) = floor(algal_cm2(id_growth,2)*convert_to_canopy);
        % A csq of this is that total benthic cover can override the area of a cell
        % As the canopy of Dict overtops corals, turf, EAM and LOB, need to re-adjust at the end 
        % what is the visible cover of each (not their true colonization)
    end
end

last_surface_area_grazed = sum(algal_cm2.*eaten_algal,2);

