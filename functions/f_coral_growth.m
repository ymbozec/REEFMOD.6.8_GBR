%__________________________________________________________________________
%
% Process coral colony growth
%
% Yves-Marie Bozec, y.bozec@uq.edu.au, created Nov 2011
%
% Modified (08/2018) to integrate genetic adaptation to warming
% _________________________________________________________________________

function [coral_cm2, genes] = f_coral_growth(coral_cm2, colony_ID, genes, species_ID, clade, ...
    macroalgal_env_prop, avail_cell_areas_cm2, nb_coral_types, CORAL, ALGAL, doing_clades, doing_genetics, genetics, SST_diff, CORAL_juvenile_growth)

id0 = zeros(size(coral_cm2)) ;

growth_rates = id0;
max_size = id0;
% over_cell = id0;

id1 = id0+1;
id1(coral_cm2<=0)=0 ; % Remove the already dead ones (also exclude zeros)

J = ones(1,size(coral_cm2,2)) ; 

macroalgal_env_prop = id1.*macroalgal_env_prop(:,J); % macroalgal environment of every colony (same env if same cell)

algal_contact = id1;
algal_contact(macroalgal_env_prop >= ALGAL.critical_algal_contact)=2;
algal_contact(macroalgal_env_prop >= ALGAL.vcritical_algal_contact)=3;

deltagrowth = id1 ;% preallocating space for the proportional reduction in coral growth due to macroalgae
% deltagrowth(coral_cm2 > 0 & coral_cm2 < CORAL.adol_size & algal_contact==2) = ALGAL.macroalgal_coral_recruit_growth_rate ;
deltagrowth(coral_cm2 > 0 & coral_cm2 < CORAL.juv_max_size & algal_contact==2) = ALGAL.macroalgal_coral_recruit_growth_rate ;
% deltagrowth(coral_cm2 > 0 & coral_cm2 < CORAL.adol_size & algal_contact==3) = 0 ;
deltagrowth(coral_cm2 > 0 & coral_cm2 < CORAL.juv_max_size & algal_contact==3) = 0 ;
% deltagrowth(coral_cm2 >= CORAL.adol_size & algal_contact==2) = ALGAL.macroalgal_coral_growth_rate ;
deltagrowth(coral_cm2 >= CORAL.juv_max_size & algal_contact==2) = ALGAL.macroalgal_coral_growth_rate ;


col_start = 1;
col_stop = 0;

% genefit = zeros(1,nb_coral_types);

if doing_genetics == 1
      
    for s = 1:nb_coral_types
        
        if species_ID(s)>0
            
            col_stop = col_stop + species_ID(s) ;
            
            id1_tmp = id1(:,col_start:col_stop);
            id0_tmp = id0(:,col_start:col_stop);
            colony_ID_tmp = colony_ID(:,col_start:col_stop);
            
            if genetics.group(s)==1
                
                % First, update the list of QTLs to only keep the survivors from previous mortality
                list_old = genes(s).list_coral_ID;
                list_new = colony_ID_tmp(id1_tmp==1);

                check = ismember(list_old,list_new) ;

                genes(s).QTLs(check==0,:,:)=[];
                genes(s).list_coral_ID(check==0,:)=[];
                genes(s).phenotypes(check==0,:)=[];
                               
                % Difference between phenotype (thermal optimum) and environment
                D = SST_diff - genes(s).phenotypes ;
                
                % Estimate relative fitness with the breadth of thermal tolerance (SIGMA)
                F_list = ones(size(D));
                F_list(D>0) = normpdf(D(D>0),0,genetics.SIGMA_HOT(s))/normpdf(0,0,genetics.SIGMA_HOT(s));
                F_list(D<0) = normpdf(D(D<0),0,genetics.SIGMA_COLD(s))/normpdf(0,0,genetics.SIGMA_COLD(s));
                                
                % Assign to the corresponding colony in the matrix
                F = id0_tmp ;
                F(id1_tmp==1) = F_list;

            else
                % If not doing genetics for that species then assume perfect fitness    
                F = id1_tmp; 
            end
            
            growth_rates(:,col_start:col_stop) = CORAL.growth_rate(s).*F;

            max_size(:,col_start:col_stop)= CORAL.adult_max_size(s);           
            col_start = col_start + species_ID(s) ;
            
        end
    end
    
else
    
    % no thermal adaptation
    for s = 1:nb_coral_types
        
        col_stop = col_stop + species_ID(s) ;
        growth_rates(:,col_start:col_stop)= CORAL.growth_rate(s);
        max_size(:,col_start:col_stop)= CORAL.adult_max_size(s);
        col_start = col_start + species_ID(s) ;
        
    end
end

% refine growth rate for juveniles (1cm per 6 month linear extension)
% growth_rates(coral_cm2 > 0 & coral_cm2 < CORAL.adol_size) = CORAL_juvenile_growth;
growth_rates(coral_cm2 > 0 & coral_cm2 < CORAL.juv_max_size) = CORAL_juvenile_growth;


if doing_clades ==1
    
    % adjust growth rate following coral clade
    id_clade2 = (clade.*id1) - id1 ; % eliminates id of clade 1 ("1")
    id_clade1 = id1 - id_clade2 ; % eliminates id of clade 2 ("2")
    growth_rates = (growth_rates .* id_clade1) + (growth_rates .* id_clade2)*CORAL.clade_reduced_growth ;
    
end

% increases overall size by one growth increment -> this gives a first potential growth for a hypothetical free space
newsize = pi*(growth_rates + sqrt(coral_cm2 .* id1/pi)).^2 ; % note final sizes will be rounded at the end

% just check if new sizes do not exceed maximum size of each species (cannot grow bigger than the max size)
check_size = max_size - newsize ;
newsize(check_size<0) = max_size(check_size<0) ; % if too big just assign the maximum size

% Now, reduce the potential growth due to algae. Note deltagrowth takes values 1 (growth is not reduced), 0.1 or 0.3
growth_cm2 = (newsize - (coral_cm2 .* id1)).*deltagrowth ;

%% %%%%TRY
% Corals fill the available space proportionately to their growth. If no available space, corals don't grow 
total_growth = sum(growth_cm2,2); % sum of all potential growth per cell
total_growth_mat = total_growth(:,J);
avail_space_mat = avail_cell_areas_cm2(:,J);

% identify cells where corals will overload the cell area
over_cell = zeros(size(total_growth));
over_cell(avail_cell_areas_cm2 - total_growth < 0) = 1 ;

id_over = over_cell(:,J) .* id1 ; % gives 1 to cells where total_growth overtakes available space
id_ok = id1 - id_over ; % identify with 1 the cells where total_growth can take the available space

% for the OK cells, all corals just grow freely
coral_cm2 = coral_cm2 + growth_cm2.*id_ok ;

% For the cell overgrown, growth is reduced proportionnally to the available space
coral_cm2(id_over==1) = coral_cm2(id_over==1) + growth_cm2(id_over==1).*avail_space_mat(id_over==1)./total_growth_mat(id_over==1);

% Finally round all the sizes after growth
% Flooring let a couple of cm2 free of corals so future recruits may settle again but won't ever grow
coral_cm2 = floor(coral_cm2);

%% %%%%%%

% identify cells where corals will overload the cell area
% over_cell = zeros(size(avail_cell_areas_cm2));
% total_growth = sum(growth_cm2,2); % sum of all potential growth per cell
% over_cell(avail_cell_areas_cm2 - total_growth < 0) = 1 ;
% 
% id_over = over_cell(:,J) .* id1 ; % gives 1 for cells where total_growth overtake available space
% id_ok = id1 - id_over ; % identify with 1 the cells where total_growth can take the available space
% 
% % Corals fill the available space proportionately to their growth. If no available space, corals don't grow 
% total_growth = total_growth(:,J).*id_over; % tranform into a matrix with same dimensions as coral colonies
% avail_cell_areas_cm2 = avail_cell_areas_cm2(:,J).*id_over ; %turn it into a matrix form
% temp_coral_cm2 = coral_cm2 + avail_cell_areas_cm2.*growth_cm2./total_growth ;
% coral_cm2(id_over==1) = temp_coral_cm2(id_over==1);
% 
% % for the OK cells, all corals just grow freely
% coral_cm2 = coral_cm2 + growth_cm2.*id_ok ;
% 
% % Finally round all the sizes after growth
% % Flooring let a couple of cm2 free of corals so future recruits may settle again but won't ever grow
% coral_cm2 = floor(coral_cm2);
