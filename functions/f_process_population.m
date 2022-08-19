% -------------------------------------------------------------------------
% Y.-M. Bozec, MSEL, created Nov 2011.
% Last modified: 21/12/2018
%
% PROCESS CORAL AND ALGAL SPATIAL DEMOGRAPHICS
% July 2018: genetic adaptation (thermal) integrated into natural mortality
% -------------------------------------------------------------------------

function [coral, algal, genes, last_surface_area_grazed] = f_process_population (coral, algal, genes, season, algal_removal, META, REEF, CORAL, ALGAL, SST_diff)

% This function processes all neighbourhood type interactions and succession for every cell
% in the grid once and once only. All the cells are processed simultaneously (vectorisation).
% Also keeps records of mortality events and calculates the fecundity total for the population
% (based on the cells after they have transformed).

% % Dimensions of the grid
% m = META.grid_x_count ;
% n = META.grid_y_count ;
% 
% % List of grazable cells (exclude sand for cell loops)
% list_cell = 1:(m*n) ;
% list_cell = list_cell(REEF.grazable_cell==1) ;

% extracting corals and algae records for internal use
algal_cm2 = full([algal.cover_cm2]);  % NOT a sparse matrix otherwise slow down the code
[coral_cm2, surface_cm2, volume_cm3, clade, colony_ID,species_ID] = f_struct_deploy(coral);
clear coral algal

% Create an id matrix (1 for every colony)
id1 = zeros(size(coral_cm2));
id1(coral_cm2 > 0) = 1; 

% total coral cover in every cell
total_coral_cm2 = sum(coral_cm2.*id1,2) ;

%_________________________________________________________________________________________
%
% SET UP THE SURROUNDING ENVIRONMENT AT t-1
%_________________________________________________________________________________________

% Coral and algal environments are set before processing - they are fully determined by covers at t-1
% NOTE: what about dictlob_cm2 in the 4 surrounding cells?
environ = META.environ; % changing name for tractability

% Total coral cover over 5 cells (cm2)
coral_env_cm2 = total_coral_cm2(environ(:,1)) + ...
    total_coral_cm2(environ(:,2)) + total_coral_cm2(environ(:,3)) + ...
    total_coral_cm2(environ(:,4)) + total_coral_cm2(environ(:,5)) ;

coral_env_cm2(REEF.grazable_cell==0)=0;

% Proportion of coral cover over 5 cells
% From now we work out the effect of actual surface area
substrate_env_cm2=REEF.substrate_SA_cm2(environ(:,1)) + ...
    REEF.substrate_SA_cm2(environ(:,2)) + REEF.substrate_SA_cm2(environ(:,3)) + ...
    REEF.substrate_SA_cm2(environ(:,4)) + REEF.substrate_SA_cm2(environ(:,5)) ;
coral_env_prop = coral_env_cm2./substrate_env_cm2 ;

% Neighborhood coral cover limit macroalgal growth 
coral_reduce_macrogrowth = coral_env_prop*ALGAL.coral_reduce_macrogrowth ;

% Algal covers over 5 cells (cm2)
% algal_env_cm2 = algal_cm2(environ(:,1),:) + ...
%     algal_cm2(environ(:,2),:) + algal_cm2(environ(:,3),:) + ...
%     algal_cm2(environ(:,4),:) + algal_cm2(environ(:,5),:) ;
% algal_env_cm2(REEF.grazable_cell==0)=0;

% Proportion (prop) of lob cover over 5 cells
% lob_env_prop = algal_env_cm2(:,3)./substrate_env_cm2;

%_________________________________________________________________________________________
%
% PROCESS ALGAE
%_________________________________________________________________________________________
% First determine which cells are fully grazable, i.e. outside Liagora canopies
grazable_cell = REEF.grazable_cell;

[algal_cm2,last_surface_area_grazed] = f_algal_dynamics(algal_removal, ALGAL, total_coral_cm2, algal_cm2, ...
    coral_reduce_macrogrowth, ALGAL.growth_rate, REEF.substrate_SA_cm2, META.convert_to_canopy ,season, REEF.dictyota_declines_seasonally,grazable_cell,REEF.substrate_SA_cm2,environ);

% Re-estimate the new colocation of DICT and LOB (positive) % TO RE-VISIT!
new_colocation_cm2 = max(sum(algal_cm2,2) + total_coral_cm2 - REEF.substrate_SA_cm2,0);

% remove used objects 
% clear coral_env_prop coral_reduce_macrogrowth  algal_env_cm2 coral_env_cm2
%_________________________________________________________________________________________
%
% PROCESS CORALS
%_________________________________________________________________________________________

%%%%%%% 1) PARROTFISH PREDATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% id_pred = id1; % id1 is id matrix of living colonies (see above)
% % assigns a random proba of mortality to every colony
% rand_mort = rand(size(id_pred));
% % switch to 0 colonies escaping predation due to large size
% id_pred(coral_cm2 > CORAL.threshold_predation_size) = 0 ;
% % % switch to 0 colonies escaping predation by chance
% id_pred(rand_mort > CORAL.parrotfish_predation) = 0 ;
% 
% % update the grid (id_pred == 1 identifies the colonies eaten)
% coral_cm2(id_pred==1) = 0 ; % Eaten colonies are not kept - flat substratum (complete removal of small colonies)

%%%%%%% 2) MACROALGAL OVERGROWTH OF CORALS %%%%%%%%%%%%%%%%%%
% Update the proportions of macroalgae over the 5 cells after grazing and algal growth
% local macroalgal abundance for competitive interactions

% New algal covers over 5 cells (cm2)
algal_env_cm2 = algal_cm2(environ(:,1),:) + ...
    algal_cm2(environ(:,2),:) + algal_cm2(environ(:,3),:) + ...
    algal_cm2(environ(:,4),:) + algal_cm2(environ(:,5),:) ;

algal_env_cm2(REEF.grazable_cell==0)=0; %just to be sure there is no algae on sand

% Calculate proportion of macroalgal cover over 5 cells
lob_env_prop = algal_env_cm2(:,3)./substrate_env_cm2; % Proportion (prop) of lob cover over 5 cells
dict_env_prop = algal_env_cm2(:,2)./substrate_env_cm2; % Proportion (prop) of dict cover over 5 cells
macroalgal_env_prop = lob_env_prop + dict_env_prop ; % Proportion (prop) of all macroalgae

% Lobophora overgrowth
if max(lob_env_prop)>0.01 % process only if at least 1 cell has >1% of lob in its 5 cells environment
    [coral_cm2, algal_cm2] = f_macroalgae_overgrowth_corals(coral_cm2, algal_cm2, species_ID, 3, lob_env_prop, CORAL.lobophora_reduce_rate, META.nb_coral_types); % /7 is unclear (mult in the reference code)
end

% Dictyota overgrowth
if max(dict_env_prop)>0.01 % process only if at least 1 cell has >1% of dict in its 5 cells environment
    [coral_cm2, algal_cm2] = f_macroalgae_overgrowth_corals(coral_cm2, algal_cm2, species_ID, 2, dict_env_prop, CORAL.dictyota_reduce_rate, META.nb_coral_types);
end
% Note that overgrowth leads to shrinking so that if a colony dies it becomes 0, not negative

    
%%%%%%% 3) CORAL GROWTH %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NOTE that the reference code compares potential coral growth with the total cell area,
% while available space should account for macroalgal cover
% Here available space excludes lob and dict (but not thick turf)
total_coral_cm2 = sum(coral_cm2.*id1,2) ; % update total coral cover in every cell
avail_cell_areas_cm2 = REEF.substrate_SA_cm2 - sum(algal_cm2(:,2:3),2) + new_colocation_cm2 - total_coral_cm2; % Exclude thick turf -> can be overgrown by corals
avail_cell_areas_cm2(REEF.grazable_cell==0)=0; % cannot grow on sand

% SST below is actually the difference between current SST and baseline SST for the reef
[coral_cm2, genes] = f_coral_growth(coral_cm2, colony_ID, genes, species_ID, clade, macroalgal_env_prop, avail_cell_areas_cm2, ...
    META.nb_coral_types, CORAL, ALGAL, META.doing_clades, META.doing_genetics, META.genetics, SST_diff, REEF.CORAL_juvenile_growth);

%%%%%%% 4) NATURAL MORTALITY   %%%%%%%%%%%%%%%%%%%%%%%%%
% Includes now (27/08/13) partial and whole colony mortalities 
[coral_cm2] = f_natural_mortality(coral_cm2, CORAL, REEF, species_ID, META.nb_coral_types);
% coral_cm2_aft2 = coral_cm2;
%_________________________________________________________________________________________
%
%       FINALLY ADJUST THE TOTALS FOR CORAL GROWTH AND MORTALITY
%_________________________________________________________________________________________

% %%%%%%% Calculate overshoot for adjustment %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
id1(coral_cm2<0)=0;  % update the identity matrix of living colonies
cell_area_cm2 = REEF.substrate_SA_cm2 ;

cell_area_cm2(REEF.grazable_cell==0)= 0;

overshoot = sum(algal_cm2,2) + sum(coral_cm2.*id1,2) - new_colocation_cm2 - cell_area_cm2 ;

% NEGATIVE OVERSHOOT means sum of all covers < total space 
% (essentially due to coral mortality that needs to be filled with EAM)
algal_cm2(overshoot<0,1) = algal_cm2(overshoot<0,1) - overshoot(overshoot<0) ;

% POSITIVE OVERSHOOT means sum of all covers > total space
% Reduce EAM first (if available)
algal_cm2(overshoot>0 & algal_cm2(:,1)>0,1) = algal_cm2(overshoot>0 & algal_cm2(:,1)>0,1)...
    - overshoot(overshoot>0 & algal_cm2(:,1)>0) ;
algal_cm2(overshoot>0 & algal_cm2(:,1)>0,1) = algal_cm2(overshoot>0 & algal_cm2(:,1)>0,1)...
    - overshoot(overshoot>0 & algal_cm2(:,1)>0) ;
algal_cm2(algal_cm2(:,1)<0,1)=0; % turn negatives into 0!
% Update overshoot
overshoot = sum(algal_cm2,2) + sum(coral_cm2.*id1,2) - new_colocation_cm2 - cell_area_cm2 ;
% Then reduce TURF if we still have positive overshoot
algal_cm2(overshoot>0 & algal_cm2(:,4)>0,4) = algal_cm2(overshoot>0 & algal_cm2(:,4)>0,4)...
    - overshoot(overshoot>0 & algal_cm2(:,4)>0) ;
algal_cm2(algal_cm2(:,4)<0,4)=0; % turn negatives into 0!
% Update overshoot
overshoot = sum(algal_cm2,2) + sum(coral_cm2.*id1,2) - new_colocation_cm2 - cell_area_cm2 ;

 
% % IF OVERSHOOT IS NEGATIVE -> fill the empty space with EAM
% % IF OVERSHOOT IS POSITIVE -> reduce EAM covers
% algal_cm2(:,1) = algal_cm2(:,1) - overshoot ;
% % This may have generated negative turf size if the sum of all other covers
% % is still too large for the cell. Then force turf to be null...
% algal_cm2(algal_cm2(:,1)<0,1) = 0 ;
% overshoot = sum(algal_cm2,2) + sum(coral_cm2.*id1,2) - new_colocation_cm2 - cell_area_cm2 ;
% 
% if sum(overshoot)>0
%     % FOR CELLS WITH POSITIVE OVERSHOOT -> first reduce TURF covers
%     algal_cm2(:,4) = algal_cm2(:,4) - overshoot ;
%     algal_cm2(algal_cm2(:,4)<0,4) = 0 ;
%     
%     overshoot = sum(algal_cm2,2) + sum(coral_cm2.*id1,2) - new_colocation_cm2 - cell_area_cm2 ;
%     
%     if sum(overshoot)>0
%         
%         list_cell = 1:(size(algal_cm2,1)) ;
% % list_cell = list_cell(REEF.grazable_cell==1) ;
%         % Then runs through hierarchical sequence adjusting macroalgae for cell size
%         list_cell2 = list_cell(overshoot(REEF.grazable_cell==1) > 0) ;
%         
%         for cell_id=list_cell2
%             [algal_cm2(cell_id,:), new_colocation_cm2(cell_id,1)] = ...
%                 f_adjust_total_cover (algal_cm2(cell_id,:), new_colocation_cm2(cell_id,1), overshoot(cell_id,1)) ;
%         end
%     end
%     
% end


% Then store the new algal cover
algal(4).cover_cm2=[];
for a=1:META.nb_algal_types
    algal(a).cover_cm2(:,1) = sparse(floor((algal_cm2(:,a)))) ;
%     
%     if isnan(sum(algal_cm2(:,a)))==1
%         stop
%     end
end


%%%%%%% Record the new coral cover %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[coral] = f_struct_rebuild (coral_cm2, surface_cm2, volume_cm3, colony_ID, clade, species_ID, META.nb_coral_types, META.doing_clades, META.doing_3D);
