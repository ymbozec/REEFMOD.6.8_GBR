function [sorted_indices, criteria, global_trigger, RESULT] = f_makeReefList(META, RESULT, t, connmetrics, remaining_COTS, thisboatorder, current_COTS_ET, COTS_site_densities)
%original code writen by Karlo Hock

%criteria used

criteria=[];
others=find(isnan(META.COTS_control_sites(:,2))==1); %% This is the index of reefs outside the MP
idx=META.cntrl_reefID; %% Reefs within MPA

switch META.COTS_reefs2cull_strat
    
    case 0%pick a reef with lots of COTS, use current COTS densities on reefs
        adult_cots_dens=sum(remaining_COTS(:,META.COTS_adult_min_age:end),2);
        [ ~, sorted_indices ] = sort( adult_cots_dens(:), 'descend' );
        criteria(:,1)= adult_cots_dens(:);
    case 1%estimated COTS larval export connectvity strategy
        %uses future matrix, knowledge of future connectivity perfect
        %pick reef based on short-range connectivity, select only on connectivity
        %see how much COTS larvae are about the be exported; this should
        %already include the size of the COTS population and the potnetial
        %to export; which means that it does not count coral on home reef
        %more, and could ignore large COTS popualtions on sinks,
        
        %Looking into the future!!
        
        
        new_density_COTS=zeros(META.nb_reefs,16);
        for rf=1:META.nb_reefs
            temp1_density_COTS = zeros(size(META.COTS_feeding_rates,2),1) ;
            density_COTS=remaining_COTS(rf,:);
            temp1_density_COTS(2:end) = squeeze(density_COTS(1:(end-1))) ;% Increment the age of all COTS before mortality
            % (note this eradicates the oldest COTS)
            new_density_COTS(rf,:) = (1-META.COTS_mortality').*temp1_density_COTS' ;
        end
        
        mature_COTS_density = squeeze(sum(new_density_COTS(:,META.COTS_fecundity~=0),2)) ;
        fertilization_success = 0.14 * ((10^8)*mature_COTS_density/META.total_area_cm2).^0.61 ;
        fertilization_success(fertilization_success>0.9)=0.9;  %fertilization cap (max obtained by Babcock)
        COTS_densities = squeeze(new_density_COTS(:,:));
        new_COTS_fecundity(:,1) = sum(squeeze(COTS_densities(1:META.nb_reefs,:)).*...
            META.COTS_fecundity(ones(1,META.nb_reefs),:),2).*fertilization_success ;
        
        COTS_output_larvae = new_COTS_fecundity.*META.area_habitat;
        
        %multiply the laral output by the average export from the reef
        centr=1;%if 1 then weighted outdegree
        cots_exportpot = COTS_output_larvae .* connmetrics.num(:,centr);
        
        [ ~, sorted_indices ] = sort( cots_exportpot(:), 'descend' );
        criteria(:,1)=new_COTS_fecundity;
        criteria(:,2)=COTS_output_larvae;
        criteria(:,3)=connmetrics.num(:,centr);
        criteria(:,4)=cots_exportpot;
        
    case 8%estimated COTS larval export connectvity strategy; DOES NOT MULTIPLY BY REEF SIZE
        %uses future matrix, knowledge of future connectivity perfect
        %pick reef based on short-range connectivity, select only on connectivity
        %see how much COTS larvae are about the be exported; this should
        %already include the size of the COTS population and the potnetial
        %to export; which means that it does not count coral on home reef
        %more, and could ignore large COTS popualtions on sinks,
        
        %Looking into the future!!
        
        
        new_density_COTS=zeros(META.nb_reefs,16);
        for rf=1:META.nb_reefs
            temp1_density_COTS = zeros(size(META.COTS_feeding_rates,2),1) ;
            density_COTS=remaining_COTS(rf,:);
            temp1_density_COTS(2:end) = squeeze(density_COTS(1:(end-1))) ;% Increment the age of all COTS before mortality
            % (note this eradicates the oldest COTS)
            new_density_COTS(rf,:) = (1-META.COTS_mortality').*temp1_density_COTS' ;
        end
        
        mature_COTS_density = squeeze(sum(new_density_COTS(:,META.COTS_fecundity~=0),2)) ;
        fertilization_success = 0.14 * ((10^8)*mature_COTS_density/META.total_area_cm2).^0.61 ;
        fertilization_success(fertilization_success>0.9)=0.9;  %fertilization cap (max obtained by Babcock)
        COTS_densities = squeeze(new_density_COTS(:,:));
        new_COTS_fecundity(:,1) = sum(squeeze(COTS_densities(1:META.nb_reefs,:)).*...
            META.COTS_fecundity(ones(1,META.nb_reefs),:),2).*fertilization_success ;
        
        COTS_output_larvae = new_COTS_fecundity;
        
        %multiply the laral output by the average export from the reef
        centr=1;%if 1 then weighted outdegree
        cots_exportpot = COTS_output_larvae .* connmetrics.num(:,centr);
        
        [ ~, sorted_indices ] = sort( cots_exportpot(:), 'descend' );
        criteria(:,1)=new_COTS_fecundity;
        criteria(:,2)=COTS_output_larvae;
        criteria(:,3)=connmetrics.num(:,centr);
        criteria(:,4)=cots_exportpot;
        
    case 2%estimated COTS larval export connectvity strategy
        %uses mean matrix, pick reef based on short-range connectivity,
        %see how much COTS larvae are about the be exported; this should
        %already include the size of the COTS population and the potnetial
        %to export; which means that it does not count coral on home reef
        %more, and could ignore large COTS popualtions on sinks,
        
        %Looking not so much into the future
        
        COTS_out=zeros(META.nb_reefs,1);
        for rf1=1:META.nb_reefs
            [COTS_out(rf1)] = sum(squeeze(remaining_COTS(rfl,:).*META.COTS_fecundity'));%used to be f_COTS_reproduction( RESULT, META, reef, time )
        end
        
        %multiply the larval output by the average export - weighted outdegree - from the reef
        centr=1;%if 1 then weighted outdegree
        cots_exportpot = COTS_out .* connmetrics.num(:,centr);
        
        [ ~, sorted_indices ] = sort( cots_exportpot(:), 'descend' );
        criteria(:,1)=COTS_out;
        criteria(:,2)=connmetrics.num(:,centr);
        criteria(:,3)=cots_exportpot;
    case 3%pick reef based on short-range connectivity
        %imperfect knowledge of connectivity, probabilistic, more realistic, but sometimes wrong and sensitive to threshold
        %sort out reefs using centrality; if 1 then weighted outdegree
        centr=1;
        [ ~, sorted_indices ] = sort( connmetrics.ranks(:,centr), 'descend' );
        criteria(:,1)=connmetrics.num(:,centr);
    case 4%pick reef based on short-range connectivity, only manage picked reefs
        %imperfect knowledge of connectivity, probabilistic, more realistic, but sometimes wrong and sensitive to threshold
        
        %sort out reefs using centrality; if 2 then number of links
        centr=2;
        [ ~, sorted_indices ] = sort( connmetrics.ranks(:,centr), 'descend' );
        criteria(:,1)=connmetrics.num(:,centr);
    case 5%Kays code
        %Use cost/benefit to select reefs chose cheap reefs first
        %get coral population
        cover = sum(squeeze(RESULT.coral_total_fecundity(:,t,:)),2);
        %get cots population weighted by size distribution
        cotsPop = squeeze(RESULT.COTS_all_densities(:,t,:));
        sizeCOTS = sum(cotsPop.* repmat(META.COTS_diameters, [META.nb_reefs,1]),2);
        %risk score is how much coral you have to lose
        %multiplied by how likely you are to lose it
        riskScore = cover.*sizeCOTS;
        
        %isolation = Minimum distance to port * mean of closest 10 reefs
        load('isolation.mat', 'isolation');
        
        %cost = isolation * size
        habA = META.area_habitat;
        cost = isolation.*habA;
        
        costBen = [cost, riskScore];
        
        %rank reefs based on risk first then cheapest as a result
        %reefs with lots to lose that don't cost much are tackeled first
        [ ~, sorted_indices ] = sortrows(costBen, [2,1], {'descend', 'ascend'});
        
    case 6
        %Use cost/benefit to select reefs chose risky reefs first
        %get coral cover
        cover = sum(squeeze(RESULT.coral_total_fecundity(:,t,:)),2);
        %get cots population weighted by size distribution
        cotsPop = squeeze(RESULT.COTS_all_densities(:,t,:));
        sizeCOTS = sum(cotsPop.* repmat(META.COTS_diameters, [META.nb_reefs,1]),2);
        %risk score is how much coral you have multiplied by how likely you are to lose it
        riskScore = cover.*sizeCOTS;
        
        %isolation = Minimum distance to port * mean of closest 10 reefs
        load('isolation.mat', 'isolation');
        
        %cost = isolation * size
        habA = META.area_habitat;
        cost = isolation.*habA;
        
        costBen = [cost, riskScore];
        
        %rank reefs based on cost first then risk as a result
        %reefs that are cheap and risky are tackeled first
        [ ~, sorted_indices ] = sortrows(costBen, [1,2], {'ascend', 'descend'});
    case 7
        %do nto use; always load list of reefs from file
        %load a fixed list of reefs with known outbreaks, only to be used for inital state obtained with hindcasts for simulation 1
        %list has 43 reefs that, without any COTS control, have tow >0.22 at t=1, have mean tow over
        %20 timesteps above 0.22, and have tow above 0.22 for more than 10
        %of the first 20 timesteps - in other words, are known to have high densities
        %and long-lasting outbreaks in simulation 1
        sio=load('sure_initial_outbreaks_S1.mat');
        sorted_indices=sio.sure_initial_outbreaks_S1;
        criteria(:,1)=sorted_indices;
        if t>1
            sorted_indices=sorted_indices(randperm(length(sorted_indices)));
        end
        criteria(:,2)=sorted_indices;
    case 9
        %fixed list, thisoboatorder determines if it permutes the list
        %load a fixed list of reefs with known outbreaks, only to be used for inital state obtained with hindcasts for simulation 1
        %list has 43 reefs that, without any COTS control, have tow >0.22 at t=1, have mean tow over
        %20 timesteps above 0.22, and have tow above 0.22 for more than 10
        %of the first 20 timesteps - in other words, are known to have high densities
        %and long-lasting outbreaks in simulation 1
        if thisboatorder==1
            sorted_indices=RESULT.COTS_inital_outbreaks_list;
        else
            iol=RESULT.COTS_inital_outbreaks_list;
            sorted_indices=iol(randperm(length(iol)));
        end
        criteria(:,1)=sorted_indices;
    case 10
        %fixed list, thisoboatorder determines if it permutes the list, always go to all reefs
        %load a fixed list of reefs with known outbreaks, only to be used for inital state obtained with hindcasts for simulation 1
        %list has 43 reefs that, without any COTS control, have tow >0.22 at t=1, have mean tow over
        %20 timesteps above 0.22, and have tow above 0.22 for more than 10
        %of the first 20 timesteps - in other words, are known to have high densities
        %and long-lasting outbreaks in simulation 1
        if thisboatorder==1
            sorted_indices=RESULT.COTS_inital_outbreaks_list;
        else
            iol=RESULT.COTS_inital_outbreaks_list;
            sorted_indices=iol(randperm(length(iol)));
        end
        criteria(:,1)=sorted_indices;
    case 11
        %fixed list, thisoboatorder determines if it permutes the list, always go to all reefs
        %load a fixed list of reefs with known outbreaks, only to be used for inital state obtained with hindcasts for simulation 1
        %list has 43 reefs that, without any COTS control, have tow >0.22 at t=1, have mean tow over
        %20 timesteps above 0.22, and have tow above 0.22 for more than 10
        %of the first 20 timesteps - in other words, are known to have high densities
        %and long-lasting outbreaks in simulation 1
        if thisboatorder==1
            sorted_indices=RESULT.COTS_inital_outbreaks_list;
        else
            iol=RESULT.COTS_inital_outbreaks_list;
            sorted_indices=iol(randperm(length(iol)));
        end
        criteria(:,1)=sorted_indices;
    case 12
        %fixed list, thisoboatorder determines if it permutes the list, always go to all reefs
        %load a fixed list of reefs with known outbreaks, only to be used for inital state obtained with hindcasts for simulation 1
        %list has 43 reefs that, without any COTS control, have tow >0.22 at t=1, have mean tow over
        %20 timesteps above 0.22, and have tow above 0.22 for more than 10
        %of the first 20 timesteps - in other words, are known to have high densities
        %and long-lasting outbreaks in simulation 1
        if thisboatorder==1
            sorted_indices=RESULT.COTS_inital_outbreaks_list;
        else
            iol=RESULT.COTS_inital_outbreaks_list;
            sorted_indices=iol(randperm(length(iol)));
        end
        criteria(:,1)=sorted_indices;
    case 13
        %attempt No 1 at CSIRO way (but see case 14 for a full implementation), first go to the reef that was visited last in order to check if
        %it is above threshold, then use the priority list and randomly
        %permute it and append it to last reef, then create the lsit of reefs outside of the priority
        %list and append it to the (last_reef+priority) list
        priority=RESULT.COTS_priority_list;
        %NOTE: create a strategy without nonprioty reefs if we exclusively want
        %to consider reefs on the priority list
        nonpriority=transpose(1:3806);%change this as needed if we only consider a certain region, i.e. then we would not create non-priority list from all 3806 GBR reefs
        nonpriority(priority)=[];
        %         if RESULT.last_reef_COTScontrolled>0
        %             nonpriority(priority)=[];
        %         else
        %             nonpriority(vertcat(priority,RESULT.last_reef_COTScontrolled))=[];%remove prioty reefs and last controlled form the list
        %         end
        if thisboatorder==1
            if RESULT.last_reef_COTScontrolled==0
                sorted_indices_all=vertcat(priority,nonpriority);
            else
                if ismember(RESULT.last_reef_COTScontrolled,priority)
                    priority(find(priority==RESULT.last_reef_COTScontrolled))=[];
                elseif ismember(RESULT.last_reef_COTScontrolled,nonpriority)
                    nonpriority(find(nonpriority==RESULT.last_reef_COTScontrolled))=[];
                end
                sorted_indices_all=vertcat(RESULT.last_reef_COTScontrolled,priority,nonpriority);
            end
        else
            if RESULT.last_reef_COTScontrolled==0
                sorted_indices_all=vertcat(priority(randperm(length(priority))),nonpriority(randperm(length(nonpriority))));
            else
                if ismember(RESULT.last_reef_COTScontrolled,priority)
                    priority(find(priority==RESULT.last_reef_COTScontrolled))=[];
                elseif ismember(RESULT.last_reef_COTScontrolled,nonpriority)
                    nonpriority(find(nonpriority==RESULT.last_reef_COTScontrolled))=[];
                end
                sorted_indices_all=vertcat(RESULT.last_reef_COTScontrolled,priority(randperm(length(priority))),nonpriority(randperm(length(nonpriority))));
            end
        end
        current_pref_coral=squeeze(sum(RESULT.coral_pct2D(:,t,1:4),3));%coral cover of preferred groups, currently 1 to 4
        ET_COTSha=0.015*((20.*current_pref_coral)+4);%coral cover-dependent ecological threshold per Condie et al
        ET_COTStow=0.015*ET_COTSha;%convert the ecological threshold to COTS per tow
        %COTS_current_tow_density=RESULT.COTS_total_perceived_density(:,t);%COTS per tow
        COTS_current_tow_density=sum(remaining_COTS(:,META.COTS_adult_min_age:end).*META.COTS_detectability(META.COTS_adult_min_age:end),2);
        sorted_COTS_density=zeros(size(sorted_indices_all,1),1);
        sorted_ET_density=zeros(size(sorted_indices_all,1),1);
        sorted_pref_coral=zeros(size(sorted_indices_all,1),1);
        for i=1:size(sorted_indices_all,1)%this assigns the apporiate values to the permuted list
            sorted_COTS_density(i,1)=COTS_current_tow_density(sorted_indices_all(i,1));
            sorted_ET_density(i,1)=ET_COTStow(sorted_indices_all(i,1));
            sorted_pref_coral(i,1)=current_pref_coral(sorted_indices_all(i,1));
        end
        reefs_below_ET=find(sorted_COTS_density<sorted_ET_density);%reefs that have COTS below ET
        sorted_indices=sorted_indices_all;
        sorted_indices(reefs_below_ET)=[];%remove the reefs below ET from active list
        criteria=zeros(size(sorted_indices_all,1),8);
        if ~isempty(sorted_indices)
            criteria(1:size(sorted_indices,1),1)=sorted_indices;%reefs above ET; only created if those exist, otherwise all zeros
        end
        criteria(:,2)=sorted_indices_all;%all reefs after perumtation, first priority list then nonpriority
        criteria(:,3)=sorted_COTS_density;%COTS density per tow that need sto correspond to list of reefs in 2nd column
        criteria(:,4)=sorted_ET_density;%COTS ET density per tow that need sto correspond to list of reefs in 2nd column
        criteria(:,5)=sorted_pref_coral;%current coral cover that need sto correspond to list of reefs in 2nd column
        criteria(:,6)=COTS_current_tow_density;%unsorted COTS ET density per tow per reef to check that the sorting was not broken
        criteria(:,7)=ET_COTStow;%unsorted COTS ET density per tow per reef to check that the sorting was not broken
        criteria(:,8)=current_pref_coral;%unsorted current coral cover per reef to check that the sorting was not broken; check against coral cover 2D
        
        %         ET_density=RESULT.COTS_total_perceived_density(:,t);
        %         sorted_density=zeros(size(sorted_indices_all,2),1);
        %         for i=1:size(sorted_indices_all,2)
        %             sorted_density(i,1)=ET_density(sorted_indices_all(i,1));
        %         end
        %         reefs_below_ET=find(sorted_density<META.COTS_ecological_threshold);
    case 14
        %CSIRO strategies that can handle all the triggers; use this as a
        %template to read different priority lists from different files eg.
        %for regional focus etc
        
        %change which file is read below to change focal region for different strategies
        if t==(META.COTS_control_start+1)%load only once to save time
            csirolist=load('CSIROprioritylist.mat');
            priority=csirolist.CSIROlist;
            others=find(isnan(META.COTS_control_sites(:,2))==1); %% Updated--This is the index of reefs outside the MP
            RESULT.COTS_priority_list(1).list=priority;
            
            nonpriority=transpose(1:3806);%change this as needed if we only consider a certain region, i.e. then we would not create non-priority list from all 3806 GBR reefs
            nonpriority=setdiff(nonpriority,others); %%Updated--remove reefs outside MPA-- control doesnt operate there
            nonpriority=setdiff(nonpriority,priority); %%remove priority
            RESULT.COTS_nonpriority_list(1).list= nonpriority;
            RESULT.COTS_effort_distribution=1;%needs to be specified for all strategies that use CSIRO
        end
        priority=RESULT.COTS_priority_list(1).list;
        nonpriority=RESULT.COTS_nonpriority_list(1).list;
        %NOTE: create a strategy without nonprioty reefs if we exclusively want
        %to consider reefs on the priority list
        
        if thisboatorder==1%this uses fixed order, no permuting of priority list
            if RESULT.last_reef_COTScontrolled==0%if we never visited any reef before, start with pure list, otherwise add last reef vissted on top of list
                sorted_indices_all=vertcat(priority,nonpriority);
            else
                if ismember(RESULT.last_reef_COTScontrolled,priority)
                    priority(find(priority==RESULT.last_reef_COTScontrolled))=[];
                elseif ismember(RESULT.last_reef_COTScontrolled,nonpriority)
                    nonpriority(find(nonpriority==RESULT.last_reef_COTScontrolled))=[];
                end
                sorted_indices_all=vertcat(RESULT.last_reef_COTScontrolled,priority,nonpriority);
            end
        else%else permute priority list
            if RESULT.last_reef_COTScontrolled==0%same as above
                sorted_indices_all=vertcat(priority(randperm(length(priority))),nonpriority(randperm(length(nonpriority))));
            else
                if ismember(RESULT.last_reef_COTScontrolled,priority)
                    priority(find(priority==RESULT.last_reef_COTScontrolled))=[];
                elseif ismember(RESULT.last_reef_COTScontrolled,nonpriority)
                    nonpriority(find(nonpriority==RESULT.last_reef_COTScontrolled))=[];
                end
                sorted_indices_all=vertcat(RESULT.last_reef_COTScontrolled,priority(randperm(length(priority))),nonpriority(randperm(length(nonpriority))));
            end
        end
        
       
        COTS_current_tow_density=sum(remaining_COTS(:,META.COTS_adult_min_age:end).*META.COTS_detectability(META.COTS_adult_min_age:end),2);%per tow density
        COTS_ctowd=COTS_current_tow_density(idx); %%Updated--only controlled reefs
        current_pref_coral=squeeze(sum(RESULT.coral_pct2D(:,t,META.COTS_pref_corals),3));%coral cover of preferred groups, currently 1 to 4
        current_pcoral=current_pref_coral(idx); %%Updated--only controlled reefs
        
        sorted_COTS_density=zeros(size(sorted_indices_all,1),1);%see below what each varaible stores
        sorted_ET_density=zeros(size(sorted_indices_all,1),1);
        sorted_pref_coral=zeros(size(sorted_indices_all,1),1);

%         for i=1:size(sorted_indices_all,1)%this assigns the appropiate values to the permuted list
%             sorted_COTS_density(i,1)=COTS_current_tow_density(sorted_indices_all(i,1));
%             sorted_ET_density(i,1)=current_COTS_ET(sorted_indices_all(i,1));
%             sorted_pref_coral(i,1)=current_pref_coral(sorted_indices_all(i,1));
%         end

        %%Updated based on the new reef list 
        for i=1:size(sorted_indices_all,1)
            id=find(idx(:,1)==sorted_indices_all(i,1));%this assigns the appropiate values to the permuted list
            sorted_COTS_density(i,1)=COTS_ctowd(id);
            sorted_ET_density(i,1)=current_COTS_ET(id);
            sorted_pref_coral(i,1)=current_pcoral(id);
        end
        
        criteria1=zeros(size(sorted_indices_all,1),8);%this stores everything we need so it is onyl calcualted once
        sorted_indices=sorted_indices_all;%lazy way to not change everything else in here
        criteria1(:,1)=sorted_indices_all;%all reefs after permutation, first priority list then nonpriority
        criteria1(:,2)=sorted_COTS_density;%COTS density per tow that need sto correspond to list of reefs in 2nd column
        criteria1(:,3)=sorted_ET_density;%COTS ET density per tow that need sto correspond to list of reefs in 2nd column
        criteria1(:,4)=sorted_pref_coral;%current coral cover that need sto correspond to list of reefs in 2nd column
        %%Updated for new reef list
        COTS_ctd=COTS_current_tow_density;
        COTS_ctd(others)=[];%% remove non-controlled reefs
        criteria1(:,5)=COTS_ctd;%unsorted COTS ET density per tow per reef to check that the sorting was not broken
        criteria1(:,6)=current_COTS_ET;
        cprefc=current_pref_coral;
        cprefc(others)=[];%% remove non-controlled reefs
        criteria1(:,7)=cprefc;
        
        %and finally check whether the strategy uses global trigger, if so check whether
        %condition is satisfied now that we know ET for each reef; note that this defaults to zero, but COTS_control
        %function only checks it if META.COTS_global_trigger=1
        global_trigger=0;
        if META.COTS_global_trigger==1%see if this strategy uses global, GBR-wide trigger
            priority_ET=current_COTS_ET(priority);
            priority_dens=COTS_current_tow_density(priority);
            priority_above_ET=priority_dens-priority_ET;
            if nnz(find(priority_above_ET>0))>round(numel(priority)*0.1)
                global_trigger=1;
            end
        end
        criteria.criteria=criteria1;
        %end of case 14 here, so copy up to this point to create new
        %strategies, change both which file a new strategy loads and how it handles nonpriority list
        
    case 15
        % Strategy to visit only North & Central. Added by Caro
        if t==(META.COTS_control_start+1)%load only once to save time
            csirolist=load('CCS2a_plist.mat');%%priority reefs only List for North & Central
            nonpriority=load('CCS2a_nplist.mat');%% Non-priority reefs for North & Central. I initially put this inside the previous loop when downloading csirolist,  but not sure that makes a difference?
            
            %others=find(isnan(META.COTS_control_sites(:,2))==1); %% This is the index of reefs outside the MP
            nonpriority= nonpriority.CCS2a_nplist;
            nonpriority=setdiff(nonpriority,others); %%remove reefs outside MPA-- control doesnt operate there
            
            priority= csirolist.CCS2a_plist;
            otherreefs=transpose(1:3806);
            otherreefs(vertcat(priority, nonpriority,others))=[];
            
            RESULT.COTS_priority_list(1).list=priority;
            RESULT.COTS_nonpriority_list(1).list= nonpriority;
            RESULT.COTS_otherreefs_list(1).list=otherreefs;
            RESULT.COTS_effort_distribution=1;%needs to be specified for all strategies that use CSIRO
        end
        priority=RESULT.COTS_priority_list(1).list;
        nonpriority=RESULT.COTS_nonpriority_list(1).list;
        otherreefs=RESULT.COTS_otherreefs_list(1).list;
        RESULT.COTS_N2consider_reefs=size(priority,1)+size(nonpriority,1);
        if thisboatorder==1%this uses fixed order, no permuting of priority list
            if RESULT.last_reef_COTScontrolled==0%if we never visited any reef before, start with pure list, otherwise add last reef vissted on top of list
                sorted_indices_all=vertcat(priority,nonpriority, otherreefs);
            else
                if ismember(RESULT.last_reef_COTScontrolled,priority)
                    priority(find(priority==RESULT.last_reef_COTScontrolled))=[];
                    sorted_indices_all=vertcat(RESULT.last_reef_COTScontrolled,priority,nonpriority, otherreefs);
                elseif ismember(RESULT.last_reef_COTScontrolled,nonpriority)
                    nonpriority(find(nonpriority==RESULT.last_reef_COTScontrolled))=[];
                    sorted_indices_all=vertcat(RESULT.last_reef_COTScontrolled,priority,nonpriority, otherreefs);
                else
                    sorted_indices_all=vertcat(priority,nonpriority, otherreefs);
                end
            end
        else%else permute priority list
            if RESULT.last_reef_COTScontrolled==0%same as above
                sorted_indices_all=vertcat(priority(randperm(length(priority))),nonpriority(randperm(length(nonpriority))),otherreefs(randperm(length(otherreefs))));
            else
                if ismember(RESULT.last_reef_COTScontrolled,priority)
                    priority(find(priority==RESULT.last_reef_COTScontrolled))=[];
                    sorted_indices_all=vertcat(RESULT.last_reef_COTScontrolled,priority(randperm(length(priority))),nonpriority(randperm(length(nonpriority))),otherreefs(randperm(length(otherreefs))));
                elseif ismember(RESULT.last_reef_COTScontrolled,nonpriority)
                    nonpriority(find(nonpriority==RESULT.last_reef_COTScontrolled))=[];
                    sorted_indices_all=vertcat(RESULT.last_reef_COTScontrolled,priority(randperm(length(priority))),nonpriority(randperm(length(nonpriority))),otherreefs(randperm(length(otherreefs))));
                else
                    sorted_indices_all=vertcat(priority(randperm(length(priority))),nonpriority(randperm(length(nonpriority))),otherreefs(randperm(length(otherreefs))));
                end
            end
        end
        
        COTS_current_tow_density=sum(remaining_COTS(:,META.COTS_adult_min_age:end).*META.COTS_detectability(META.COTS_adult_min_age:end),2);%per tow density
        COTS_ctowd=COTS_current_tow_density(idx); %%Updated--only controlled reefs
        current_pref_coral=squeeze(sum(RESULT.coral_pct2D(:,t,META.COTS_pref_corals),3));%coral cover of preferred groups, currently 1 to 4
        current_pcoral=current_pref_coral(idx); %%Updated--only controlled reefs

        sorted_COTS_density=zeros(size(sorted_indices_all,1),1);%see below what each varaible stores
        sorted_ET_density=zeros(size(sorted_indices_all,1),1);
        sorted_pref_coral=zeros(size(sorted_indices_all,1),1);
%         for i=1:size(sorted_indices_all,1)%this assigns the appropiate values to the permuted list
%             sorted_COTS_density(i,1)=COTS_current_tow_density(sorted_indices_all(i,1));
%             sorted_ET_density(i,1)=current_COTS_ET(sorted_indices_all(i,1));
%             sorted_pref_coral(i,1)=current_pref_coral(sorted_indices_all(i,1));
%         end

        %%Updated based on the new reef list 
        for i=1:size(sorted_indices_all,1)
            id=find(idx(:,1)==sorted_indices_all(i,1));%this assigns the appropiate values to the permuted list
            sorted_COTS_density(i,1)=COTS_ctowd(id);
            sorted_ET_density(i,1)=current_COTS_ET(id);
            sorted_pref_coral(i,1)=current_pcoral(id);
        end

        criteria1=zeros(size(sorted_indices_all,1),8);%this stores everythign we need so it is onyl calcualted once
        sorted_indices=sorted_indices_all;%lazy way to not change everythign else in here
        criteria1(:,1)=sorted_indices_all;%all reefs after permutation, first priority list then nonpriority
        criteria1(:,2)=sorted_COTS_density;%COTS density per tow that need sto correspond to list of reefs in 2nd column
        criteria1(:,3)=sorted_ET_density;%COTS ET density per tow that need sto correspond to list of reefs in 2nd column
        criteria1(:,4)=sorted_pref_coral;%current coral cover that need sto correspond to list of reefs in 2nd column
        %%Updated for new reef list
        COTS_ctd=COTS_current_tow_density;
        COTS_ctd(others)=[];%% remove non-controlled reefs
        criteria1(:,5)=COTS_ctd;%unsorted COTS ET density per tow per reef to check that the sorting was not broken
        criteria1(:,6)=current_COTS_ET;
        cprefc=current_pref_coral;
        cprefc(others)=[];%% remove non-controlled reefs
        criteria1(:,7)=cprefc;
        
        %and finally check whether the strategy uses global trigger, if so check whether
        %condition is satisfied now that we know ET for each reef; note that this defaults to zero, but COTS_control
        %function only checks it if META.COTS_global_trigger=1
        global_trigger=0;
        if META.COTS_global_trigger==1%see if this strategy uses global, GBR-wide trigger
            priority_ET=current_COTS_ET(priority);
            priority_dens=COTS_current_tow_density(priority);
            priority_above_ET=priority_dens-priority_ET;
            if nnz(find(priority_above_ET>0))>round(numel(priority)*0.1)
                global_trigger=1;
            end
        end
        criteria.criteria=criteria1;
        
    case 16
        %CSIRO strategy that goes to priority reefs, then also goes to 1' lat from priority reefs with outbreaks
        %this one is for the whole GBR
        
        %change which file is read below to change focal region for different strategies
        if t==(META.COTS_control_start+1)%load only once to save time
            csirolist=load('CSIROprioritylist.mat');
            priority=csirolist.CSIROlist;
            RESULT.COTS_priority_list(1).list=priority;
            
            RESULT.COTS_effort_distribution=1;%needs to be specified for all strategies that use CSIRO
        end
        priority=RESULT.COTS_priority_list(1).list;
        
        %first, find all priority reefs with sites that are above ET
        priority_outbreaks=[];
        for rfs=1:size(RESULT.COTS_priority_list(1).list,1)
            this_reef=RESULT.COTS_priority_list(1).list(rfs);
            this_reef_sites=COTS_site_densities{this_reef,1};
            this_sites_tow=zeros(size(this_reef_sites,1),1);
            this_sites_over_ET=zeros(size(this_reef_sites,1),1);
            for sts=1:size(this_reef_sites,1)%convert COTS at sites to COTS per tow, also check which ones are over ET
                this_sites_tow(sts,1)=sum(this_reef_sites(sts,META.COTS_adult_min_age:end).*META.COTS_detectability(META.COTS_adult_min_age:end),2);
                if this_sites_tow(sts,1)>current_COTS_ET(this_reef,1)
                    this_sites_over_ET(sts,1)=1;
                end
            end
            sites_over_ET=find(this_sites_over_ET);
            if ~isempty(sites_over_ET)
                priority_outbreaks=vertcat(priority_outbreaks,this_reef);
            end
        end
        
        %then, find all reefs that are within 1' latititude N-S from outbreak reefs
        %Update-- exclude reefs outside MPA
        META.reef_lat(others)=[];
        
        nonpriority=[];
        if ~isempty(priority_outbreaks)
            for pob=1:size(priority_outbreaks,1)
                this_pri_ob=priority_outbreaks(pob,1);
                this_lat=META.reef_lat(this_pri_ob,1);
                this_np=find(abs(this_lat-META.reef_lat(:,1))<1);
                if ~isempty(this_np)
                    nonpriority=vertcat(nonpriority,this_np);
                end
            end
        end
        nonpriority=unique(nonpriority);%clean it up
        nonpriority=setdiff(nonpriority,priority);%remove the reefs that are already in priority list if they were within 1' of another priority reef with outbreak
        others=find(isnan(META.COTS_control_sites(:,2))==1); %% This is the index of reefs outside the MP
        nonpriority=setdiff(nonpriority,others); %%remove reefs outside MPA-- control doesnt operate there
        
        otherreefs=transpose(1:3806);
        pnp=vertcat(priority, nonpriority);
        otherreefs(pnp)=[];
        otherreefs=setdiff(otherreefs,others);%%remove reefs outside MPA-- control doesnt operate ther
        
        RESULT.COTS_nonpriority_list(1).list= nonpriority;
        RESULT.COTS_dynamic_nonpriority_list(t).npl=nonpriority;
        
        RESULT.COTS_N2consider_reefs=size(priority,1)+size(nonpriority,1);
        if thisboatorder==1%this uses fixed order, no permuting of priority list
            if RESULT.last_reef_COTScontrolled==0%if we never visited any reef before, start with pure list, otherwise add last reef vissted on top of list
                sorted_indices_all=vertcat(priority,nonpriority, otherreefs);
            else
                if ismember(RESULT.last_reef_COTScontrolled,priority)
                    priority(find(priority==RESULT.last_reef_COTScontrolled))=[];
                    sorted_indices_all=vertcat(RESULT.last_reef_COTScontrolled,priority,nonpriority, otherreefs);
                elseif ismember(RESULT.last_reef_COTScontrolled,nonpriority)
                    nonpriority(find(nonpriority==RESULT.last_reef_COTScontrolled))=[];
                    sorted_indices_all=vertcat(RESULT.last_reef_COTScontrolled,priority,nonpriority, otherreefs);
                else
                    sorted_indices_all=vertcat(priority,nonpriority, otherreefs);
                end
            end
        else%else permute priority list
            if RESULT.last_reef_COTScontrolled==0%same as above
                sorted_indices_all=vertcat(priority(randperm(length(priority))),nonpriority(randperm(length(nonpriority))),otherreefs(randperm(length(otherreefs))));
            else
                if ismember(RESULT.last_reef_COTScontrolled,priority)
                    priority(find(priority==RESULT.last_reef_COTScontrolled))=[];
                    sorted_indices_all=vertcat(RESULT.last_reef_COTScontrolled,priority(randperm(length(priority))),nonpriority(randperm(length(nonpriority))),otherreefs(randperm(length(otherreefs))));
                elseif ismember(RESULT.last_reef_COTScontrolled,nonpriority)
                    nonpriority(find(nonpriority==RESULT.last_reef_COTScontrolled))=[];
                    sorted_indices_all=vertcat(RESULT.last_reef_COTScontrolled,priority(randperm(length(priority))),nonpriority(randperm(length(nonpriority))),otherreefs(randperm(length(otherreefs))));
                else
                    sorted_indices_all=vertcat(priority(randperm(length(priority))),nonpriority(randperm(length(nonpriority))),otherreefs(randperm(length(otherreefs))));
                end
            end
        end
        
        COTS_current_tow_density=sum(remaining_COTS(:,META.COTS_adult_min_age:end).*META.COTS_detectability(META.COTS_adult_min_age:end),2);%per tow density
        COTS_ctowd=COTS_current_tow_density(idx); %%Updated--only controlled reefs
        current_pref_coral=squeeze(sum(RESULT.coral_pct2D(:,t,META.COTS_pref_corals),3));%coral cover of preferred groups, currently 1 to 4
        current_pcoral=current_pref_coral(idx); %%Updated--only controlled reefs

        sorted_COTS_density=zeros(size(sorted_indices_all,1),1);%see below what each varaible stores
        sorted_ET_density=zeros(size(sorted_indices_all,1),1);
        sorted_pref_coral=zeros(size(sorted_indices_all,1),1);
  
        %         for i=1:size(sorted_indices_all,1)%this assigns the appropiate values to the permuted list
%             sorted_COTS_density(i,1)=COTS_current_tow_density(sorted_indices_all(i,1));
%             sorted_ET_density(i,1)=current_COTS_ET(sorted_indices_all(i,1));
%             sorted_pref_coral(i,1)=current_pref_coral(sorted_indices_all(i,1));
%         end

        %%Updated based on the new reef list 
        for i=1:size(sorted_indices_all,1)
            id=find(idx(:,1)==sorted_indices_all(i,1));%this assigns the appropiate values to the permuted list
            sorted_COTS_density(i,1)=COTS_ctowd(id);
            sorted_ET_density(i,1)=current_COTS_ET(id);
            sorted_pref_coral(i,1)=current_pcoral(id);
        end
        
        criteria1=zeros(size(sorted_indices_all,1),8);%this stores everythign we need so it is onyl calcualted once
        sorted_indices=sorted_indices_all;%lazy way to not change everythign else in here
        criteria1(:,1)=sorted_indices_all;%all reefs after perumtation, first priority list then nonpriority
        criteria1(:,2)=sorted_COTS_density;%COTS density per tow that need sto correspond to list of reefs in 2nd column
        criteria1(:,3)=sorted_ET_density;%COTS ET density per tow that need sto correspond to list of reefs in 2nd column
        criteria1(:,4)=sorted_pref_coral;%current coral cover that need sto correspond to list of reefs in 2nd column
         %%Updated for new reef list
        COTS_ctd=COTS_current_tow_density;
        COTS_ctd(others)=[];%% remove non-controlled reefs
        criteria1(:,5)=COTS_ctd;%unsorted COTS ET density per tow per reef to check that the sorting was not broken
        criteria1(:,6)=current_COTS_ET;
        cprefc=current_pref_coral;
        cprefc(others)=[];%% remove non-controlled reefs
        criteria1(:,7)=cprefc;
        
        %and finally check whether the strategy uses global trigger, if so check whether
        %condition is satisfied now that we know ET for each reef; note that this defaults to zero, but COTS_control
        %function only checks it if META.COTS_global_trigger=1
        global_trigger=0;
        if META.COTS_global_trigger==1%see if this strategy uses global, GBR-wide trigger
            priority_ET=current_COTS_ET(priority);
            priority_dens=COTS_current_tow_density(priority);
            priority_above_ET=priority_dens-priority_ET;
            if nnz(find(priority_above_ET>0))>round(numel(priority)*0.1)
                global_trigger=1;
            end
        end
        criteria.criteria=criteria1;
        
    case 17
        %CSIRO strategy that goes to priority reefs, then also goes to 1' lat from priority reefs with outbreaks
        %this one is for the different parts of the GBR
        
        
        if t==(META.COTS_control_start+1)%load only once to save time
            northcentr=load('CCS2a_plist.mat');%%priority reefs only List for North & Central
            northcentr=northcentr.CCS2a_plist;
            RESULT.COTS_priority_list(1).list=northcentr;
            southcentr=load('CCS2c_plist.mat');%%priority reefs only List for North & Central
            southcentr=southcentr.CCS2c_plist;
            south=setdiff(southcentr,northcentr);
            RESULT.COTS_priority_list(2).list=south;
           
            RESULT.COTS_effort_distribution=[0.75 0.25];%needs to be specified for all strategies that use CSIRO
        end
        for edb=1:length(RESULT.COTS_effort_distribution)
            priority=RESULT.COTS_priority_list(edb).list;
            
            %first, find all priority reefs with sites that are above ET
            priority_outbreaks=[];
            for rfs=1:size(priority,1)
                this_reef=priority(rfs);
                this_reef_sites=COTS_site_densities{this_reef,1};
                this_sites_tow=zeros(size(this_reef_sites,1),1);
                this_sites_over_ET=zeros(size(this_reef_sites,1),1);
                for sts=1:size(this_reef_sites,1)%convert COTS at sites to COTS per tow, also check which ones are over ET
                    this_sites_tow(sts,1)=sum(this_reef_sites(sts,META.COTS_adult_min_age:end).*META.COTS_detectability(META.COTS_adult_min_age:end),2);
                    if this_sites_tow(sts,1)>current_COTS_ET(this_reef,1)
                        this_sites_over_ET(sts,1)=1;
                    end
                end
                sites_over_ET=find(this_sites_over_ET);
                if ~isempty(sites_over_ET)
                    priority_outbreaks=vertcat(priority_outbreaks,this_reef);
                end
            end
            
            %then, find all reefs that are within 1' latititude N-S from outbreak reefs
            nonpriority=[];
            if ~isempty(priority_outbreaks)
                for pob=1:size(priority_outbreaks,1)
                    this_pri_ob=priority_outbreaks(pob,1);
                    this_lat=META.reef_lat(this_pri_ob,1);
                    this_np=find(abs(this_lat-META.reef_lat(:,1))<1);
                    if ~isempty(this_np)
                        nonpriority=vertcat(nonpriority,this_np);
                    end
                end
            end
            nonpriority=unique(nonpriority);%clean it up
            nonpriority=setdiff(nonpriority,priority);%remove the reefs that are already in priority list if they were within 1' of another priority reef with outbreak
            others=find(isnan(META.COTS_control_sites(:,2))==1); %% This is the index of reefs outside the MP
            nonpriority=setdiff(nonpriority,others); %%remove reefs outside MPA-- control doesnt operate ther
            
            otherreefs=transpose(1:3806);
            pnp=vertcat(priority, nonpriority);
            otherreefs(pnp)=[];
            otherreefs=setdiff(otherreefs,others);%%remove reefs outside MPA-- control doesnt operate ther
            
            RESULT.COTS_nonpriority_list(edb).list= nonpriority;
            %RESULT.COTS_dynamic_nonpriority_list(t).npl=nonpriority;
            RESULT.COTS_N2consider_reefs=size(priority,1)+size(nonpriority,1);
            if thisboatorder==1%this uses fixed order, no permuting of priority list
                if RESULT.last_reef_COTScontrolled==0%if we never visited any reef before, start with pure list, otherwise add last reef vissted on top of list
                    sorted_indices_all=vertcat(priority,nonpriority, otherreefs);
                else
                    if ismember(RESULT.last_reef_COTScontrolled,priority)
                        priority(find(priority==RESULT.last_reef_COTScontrolled))=[];
                        sorted_indices_all=vertcat(RESULT.last_reef_COTScontrolled,priority,nonpriority, otherreefs);
                    elseif ismember(RESULT.last_reef_COTScontrolled,nonpriority)
                        nonpriority(find(nonpriority==RESULT.last_reef_COTScontrolled))=[];
                        sorted_indices_all=vertcat(RESULT.last_reef_COTScontrolled,priority,nonpriority, otherreefs);
                    else
                        sorted_indices_all=vertcat(priority,nonpriority, otherreefs);
                    end
                end
            else%else permute priority list
                if RESULT.last_reef_COTScontrolled==0%same as above
                    sorted_indices_all=vertcat(priority(randperm(length(priority))),nonpriority(randperm(length(nonpriority))),otherreefs(randperm(length(otherreefs))));
                else
                    if ismember(RESULT.last_reef_COTScontrolled,priority)
                        priority(find(priority==RESULT.last_reef_COTScontrolled))=[];
                        sorted_indices_all=vertcat(RESULT.last_reef_COTScontrolled,priority(randperm(length(priority))),nonpriority(randperm(length(nonpriority))),otherreefs(randperm(length(otherreefs))));
                    elseif ismember(RESULT.last_reef_COTScontrolled,nonpriority)
                        nonpriority(find(nonpriority==RESULT.last_reef_COTScontrolled))=[];
                        sorted_indices_all=vertcat(RESULT.last_reef_COTScontrolled,priority(randperm(length(priority))),nonpriority(randperm(length(nonpriority))),otherreefs(randperm(length(otherreefs))));
                    else
                        sorted_indices_all=vertcat(priority(randperm(length(priority))),nonpriority(randperm(length(nonpriority))),otherreefs(randperm(length(otherreefs))));
                    end
                end
            end
            COTS_current_tow_density=sum(remaining_COTS(:,META.COTS_adult_min_age:end).*META.COTS_detectability(META.COTS_adult_min_age:end),2);%per tow density
            current_pref_coral=squeeze(sum(RESULT.coral_pct2D(:,t,META.COTS_pref_corals),3));%coral cover of preferred groups, currently 1 to 4
            sorted_COTS_density=zeros(size(sorted_indices_all,1),1);%see below what each varaible stores
            sorted_ET_density=zeros(size(sorted_indices_all,1),1);
            sorted_pref_coral=zeros(size(sorted_indices_all,1),1);
            for i=1:size(sorted_indices_all,1)%this assigns the appropiate values to the permuted list
                sorted_COTS_density(i,1)=COTS_current_tow_density(sorted_indices_all(i,1));
                sorted_ET_density(i,1)=current_COTS_ET(sorted_indices_all(i,1));
                sorted_pref_coral(i,1)=current_pref_coral(sorted_indices_all(i,1));
            end
            criteria1=zeros(size(sorted_indices_all,1),8);%this stores everythign we need so it is onyl calcualted once
            sorted_indices=sorted_indices_all;%lazy way to not change everythign else in here
            criteria1(:,1)=sorted_indices_all;%all reefs after perumtation, first priority list then nonpriority
            criteria1(:,2)=sorted_COTS_density;%COTS density per tow that need sto correspond to list of reefs in 2nd column
            criteria1(:,3)=sorted_ET_density;%COTS ET density per tow that need sto correspond to list of reefs in 2nd column
            criteria1(:,4)=sorted_pref_coral;%current coral cover that need sto correspond to list of reefs in 2nd column
            criteria1(:,5)=COTS_current_tow_density;%unsorted COTS ET density per tow per reef to check that the sorting was not broken
            criteria1(:,6)=current_COTS_ET;%unsorted COTS ET density per tow per reef to check that the sorting was not broken
            criteria1(:,7)=current_pref_coral;%unsorted current coral cover per reef to check that the sorting was not broken; check against coral cover 2D
            criteria(edb).criteria=criteria1;
        end
        %and finally check whether the strategy uses global trigger, if so check whether
        %condition is satisfied now that we know ET for each reef; note that this defaults to zero, but COTS_control
        %function only checks it if META.COTS_global_trigger=1
        global_trigger=0;
        if META.COTS_global_trigger==1%see if this strategy uses global, GBR-wide trigger
            priority_ET=current_COTS_ET(priority);
            priority_dens=COTS_current_tow_density(priority);
            priority_above_ET=priority_dens-priority_ET;
            if nnz(find(priority_above_ET>0))>round(numel(priority)*0.1)
                global_trigger=1;
            end
        end
    case 18
        % Strategy to visit only  Central (CCS-2b)
        if t==(META.COTS_control_start+1)%load only once to save time
            csirolist=load('CCS2b_plist.mat');%%priority reefs only List forCentral
            nonpriority=load('CCS2b_nplist.mat');%% Non-priority reefs forCentral. 
            priority= csirolist.CCS2b_plist;
            nonpriority= nonpriority.CCS2b_nplist;
            others=find(isnan(META.COTS_control_sites(:,2))==1); %% This is the index of reefs outside the MP
            nonpriority=setdiff(nonpriority,others); %%remove reefs outside MPA-- control doesnt operate ther
            
            
            otherreefs=transpose(1:3806);
            otherreefs(vertcat(priority, nonpriority))=[];
            otherreefs=setdiff(otherreefs,others); %%remove reefs outside MPA-- control doesnt operate ther
            
            
            RESULT.COTS_priority_list(1).list=priority;
            RESULT.COTS_nonpriority_list(1).list= nonpriority;
            RESULT.COTS_otherreefs_list(1).list=otherreefs;
            RESULT.COTS_effort_distribution=1;%needs to be specified for all strategies that use CSIRO
        end
        priority=RESULT.COTS_priority_list(1).list;
        nonpriority=RESULT.COTS_nonpriority_list(1).list;
        otherreefs=RESULT.COTS_otherreefs_list(1).list;
        RESULT.COTS_N2consider_reefs=size(priority,1)+size(nonpriority,1);
        
        if thisboatorder==1%this uses fixed order, no permuting of priority list
            if RESULT.last_reef_COTScontrolled==0%if we never visited any reef before, start with pure list, otherwise add last reef vissted on top of list
                sorted_indices_all=vertcat(priority,nonpriority, otherreefs);
            else
                if ismember(RESULT.last_reef_COTScontrolled,priority)
                    priority(find(priority==RESULT.last_reef_COTScontrolled))=[];
                    sorted_indices_all=vertcat(RESULT.last_reef_COTScontrolled,priority,nonpriority, otherreefs);
                elseif ismember(RESULT.last_reef_COTScontrolled,nonpriority)
                    nonpriority(find(nonpriority==RESULT.last_reef_COTScontrolled))=[];
                    sorted_indices_all=vertcat(RESULT.last_reef_COTScontrolled,priority,nonpriority, otherreefs);
                else
                    sorted_indices_all=vertcat(priority,nonpriority, otherreefs);
                end
            end
        else%else permute priority list
            if RESULT.last_reef_COTScontrolled==0%same as above
                sorted_indices_all=vertcat(priority(randperm(length(priority))),nonpriority(randperm(length(nonpriority))),otherreefs(randperm(length(otherreefs))));
            else
                if ismember(RESULT.last_reef_COTScontrolled,priority)
                    priority(find(priority==RESULT.last_reef_COTScontrolled))=[];
                    sorted_indices_all=vertcat(RESULT.last_reef_COTScontrolled,priority(randperm(length(priority))),nonpriority(randperm(length(nonpriority))),otherreefs(randperm(length(otherreefs))));
                elseif ismember(RESULT.last_reef_COTScontrolled,nonpriority)
                    nonpriority(find(nonpriority==RESULT.last_reef_COTScontrolled))=[];
                    sorted_indices_all=vertcat(RESULT.last_reef_COTScontrolled,priority(randperm(length(priority))),nonpriority(randperm(length(nonpriority))),otherreefs(randperm(length(otherreefs))));
                else
                    sorted_indices_all=vertcat(priority(randperm(length(priority))),nonpriority(randperm(length(nonpriority))),otherreefs(randperm(length(otherreefs))));
                end
            end
        end
        
        COTS_current_tow_density=sum(remaining_COTS(:,META.COTS_adult_min_age:end).*META.COTS_detectability(META.COTS_adult_min_age:end),2);%per tow density
        COTS_ctowd=COTS_current_tow_density(idx); %%Updated--only controlled reefs
        current_pref_coral=squeeze(sum(RESULT.coral_pct2D(:,t,META.COTS_pref_corals),3));%coral cover of preferred groups, currently 1 to 4
        current_pcoral=current_pref_coral(idx); %%Updated--only controlled reefs
        sorted_COTS_density=zeros(size(sorted_indices_all,1),1);%see below what each variable stores
        sorted_ET_density=zeros(size(sorted_indices_all,1),1);
        sorted_pref_coral=zeros(size(sorted_indices_all,1),1);
%         for i=1:size(sorted_indices_all,1)%this assigns the appropiate values to the permuted list
%             sorted_COTS_density(i,1)=COTS_current_tow_density(sorted_indices_all(i,1));
%             sorted_ET_density(i,1)=current_COTS_ET(sorted_indices_all(i,1));
%             sorted_pref_coral(i,1)=current_pref_coral(sorted_indices_all(i,1));
%         end

%%Updated based on the new reef list 
        for i=1:size(sorted_indices_all,1)
            id=find(idx(:,1)==sorted_indices_all(i,1));%this assigns the appropiate values to the permuted list
            sorted_COTS_density(i,1)=COTS_ctowd(id);
            sorted_ET_density(i,1)=current_COTS_ET(id);
            sorted_pref_coral(i,1)=current_pcoral(id);
        end
        
        criteria1=zeros(size(sorted_indices_all,1),8);%this stores everythign we need so it is onyl calcualted once
        sorted_indices=sorted_indices_all;%lazy way to not change everythign else in here
        criteria1(:,1)=sorted_indices_all;%all reefs after perumtation, first priority list then nonpriority
        criteria1(:,2)=sorted_COTS_density;%COTS density per tow that needs to correspond to list of reefs in 2nd column
        criteria1(:,3)=sorted_ET_density;%COTS ET density per tow that needs to correspond to list of reefs in 2nd column
        criteria1(:,4)=sorted_pref_coral;%current coral cover that need sto correspond to list of reefs in 2nd column
        COTS_ctd=COTS_current_tow_density;
        COTS_ctd(others)=[];%% remove non-controlled reefs
        criteria1(:,5)=COTS_ctd;%unsorted COTS ET density per tow per reef to check that the sorting was not broken
        criteria1(:,6)=current_COTS_ET;
        cprefc=current_pref_coral;
        cprefc(others)=[];%% remove non-controlled reefs
        criteria1(:,7)=cprefc;
        %and finally check whether the strategy uses global trigger, if so check whether
        %condition is satisfied now that we know ET for each reef; note that this defaults to zero, but COTS_control
        %function only checks it if META.COTS_global_trigger=1
        global_trigger=0;
        if META.COTS_global_trigger==1%see if this strategy uses global, GBR-wide trigger
            priority_ET=current_COTS_ET(priority);
            priority_dens=COTS_current_tow_density(priority);
            priority_above_ET=priority_dens-priority_ET;
            if nnz(find(priority_above_ET>0))>round(numel(priority)*0.1)
                global_trigger=1;
            end
        end
        criteria.criteria=criteria1;
        
    case 19
        % Strategy to visit only  Central & South (CCS-2c) Added by Caro
        if t==(META.COTS_control_start+1)%load only once to save time
            csirolist=load('CCS2c_plist.mat');%%priority reefs only List for Central & South
            nonpriority=load('CCS2c_nplist.mat');%% Non-priority reefs for Central & South
            priority= csirolist.CCS2c_plist;
            nonpriority= nonpriority.CCS2c_nplist;
            others=find(isnan(META.COTS_control_sites(:,2))==1); %% This is the index of reefs outside the MP
            nonpriority=setdiff(nonpriority,others); %%remove reefs outside MPA-- control doesnt operate there
            
            otherreefs=transpose(1:3806);
            otherreefs(vertcat(priority, nonpriority))=[];
            otherreefs=setdiff(otherreefs,others); %%remove reefs outside MPA-- control doesnt operate ther
            
            RESULT.COTS_priority_list(1).list=priority;
            RESULT.COTS_nonpriority_list(1).list= nonpriority;
            RESULT.COTS_otherreefs_list(1).list=otherreefs;
            RESULT.COTS_effort_distribution=1;%needs to be specified for all strategies that use CSIRO
        end
        priority=RESULT.COTS_priority_list(1).list;
        nonpriority=RESULT.COTS_nonpriority_list(1).list;
        otherreefs=RESULT.COTS_otherreefs_list(1).list;
        RESULT.COTS_N2consider_reefs=size(priority,1)+size(nonpriority,1);
        if thisboatorder==1%this uses fixed order, no permuting of priority list
            if RESULT.last_reef_COTScontrolled==0%if we never visited any reef before, start with pure list, otherwise add last reef vissted on top of list
                sorted_indices_all=vertcat(priority,nonpriority, otherreefs);
            else
                if ismember(RESULT.last_reef_COTScontrolled,priority)
                    priority(find(priority==RESULT.last_reef_COTScontrolled))=[];
                    sorted_indices_all=vertcat(RESULT.last_reef_COTScontrolled,priority,nonpriority, otherreefs);
                elseif ismember(RESULT.last_reef_COTScontrolled,nonpriority)
                    nonpriority(find(nonpriority==RESULT.last_reef_COTScontrolled))=[];
                    sorted_indices_all=vertcat(RESULT.last_reef_COTScontrolled,priority,nonpriority, otherreefs);
                else
                    sorted_indices_all=vertcat(priority,nonpriority, otherreefs);
                end
            end
        else%else permute priority list
            if RESULT.last_reef_COTScontrolled==0%same as above
                sorted_indices_all=vertcat(priority(randperm(length(priority))),nonpriority(randperm(length(nonpriority))),otherreefs(randperm(length(otherreefs))));
            else
                if ismember(RESULT.last_reef_COTScontrolled,priority)
                    priority(find(priority==RESULT.last_reef_COTScontrolled))=[];
                    sorted_indices_all=vertcat(RESULT.last_reef_COTScontrolled,priority(randperm(length(priority))),nonpriority(randperm(length(nonpriority))),otherreefs(randperm(length(otherreefs))));
                elseif ismember(RESULT.last_reef_COTScontrolled,nonpriority)
                    nonpriority(find(nonpriority==RESULT.last_reef_COTScontrolled))=[];
                    sorted_indices_all=vertcat(RESULT.last_reef_COTScontrolled,priority(randperm(length(priority))),nonpriority(randperm(length(nonpriority))),otherreefs(randperm(length(otherreefs))));
                else
                    sorted_indices_all=vertcat(priority(randperm(length(priority))),nonpriority(randperm(length(nonpriority))),otherreefs(randperm(length(otherreefs))));
                end
            end
        end
       COTS_current_tow_density=sum(remaining_COTS(:,META.COTS_adult_min_age:end).*META.COTS_detectability(META.COTS_adult_min_age:end),2);%per tow density
        COTS_ctowd=COTS_current_tow_density(idx); %%Updated--only controlled reefs
        current_pref_coral=squeeze(sum(RESULT.coral_pct2D(:,t,META.COTS_pref_corals),3));%coral cover of preferred groups, currently 1 to 4
        current_pcoral=current_pref_coral(idx); %%Updated--only controlled reefs
        
        sorted_COTS_density=zeros(size(sorted_indices_all,1),1);%see below what each varaible stores
        sorted_ET_density=zeros(size(sorted_indices_all,1),1);
        sorted_pref_coral=zeros(size(sorted_indices_all,1),1);
        
%         for i=1:size(sorted_indices_all,1)%this assigns the appropiate values to the permuted list
%             sorted_COTS_density(i,1)=COTS_current_tow_density(sorted_indices_all(i,1));
%             sorted_ET_density(i,1)=current_COTS_ET(sorted_indices_all(i,1));
%             sorted_pref_coral(i,1)=current_pref_coral(sorted_indices_all(i,1));
%         end

%%Updated based on the new reef list 
        for i=1:size(sorted_indices_all,1)
            id=find(idx(:,1)==sorted_indices_all(i,1));%this assigns the appropiate values to the permuted list
            sorted_COTS_density(i,1)=COTS_ctowd(id);
            sorted_ET_density(i,1)=current_COTS_ET(id);
            sorted_pref_coral(i,1)=current_pcoral(id);
        end

        criteria1=zeros(size(sorted_indices_all,1),8);%this stores everythign we need so it is onyl calcualted once
        sorted_indices=sorted_indices_all;%lazy way to not change everythign else in here
        criteria1(:,1)=sorted_indices_all;%all reefs after perumtation, first priority list then nonpriority
        criteria1(:,2)=sorted_COTS_density;%COTS density per tow that need sto correspond to list of reefs in 2nd column
        criteria1(:,3)=sorted_ET_density;%COTS ET density per tow that need sto correspond to list of reefs in 2nd column
        criteria1(:,4)=sorted_pref_coral;%current coral cover that need sto correspond to list of reefs in 2nd column
        COTS_ctd=COTS_current_tow_density;
        COTS_ctd(others)=[];%% remove non-controlled reefs
        criteria1(:,5)=COTS_ctd;%unsorted COTS ET density per tow per reef to check that the sorting was not broken
        criteria1(:,6)=current_COTS_ET;
        cprefc=current_pref_coral;
        cprefc(others)=[];%% remove non-controlled reefs
        criteria1(:,7)=cprefc;
        %and finally check whether the strategy uses global trigger, if so check whether
        %condition is satisfied now that we know ET for each reef; note that this defaults to zero, but COTS_control
        %function only checks it if META.COTS_global_trigger=1
        global_trigger=0;
        if META.COTS_global_trigger==1%see if this strategy uses global, GBR-wide trigger
            priority_ET=current_COTS_ET(priority);
            priority_dens=COTS_current_tow_density(priority);
            priority_above_ET=priority_dens-priority_ET;
            if nnz(find(priority_above_ET>0))>round(numel(priority)*0.1)
                global_trigger=1;
            end
        end
        criteria.criteria=criteria1;
end