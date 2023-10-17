function [sorted_indices, criteria, global_trigger, RESULT] = f_makeReefList_CCS(META, RESULT, t, connmetrics, remaining_COTS, thisboatorder, current_COTS_ET, COTS_site_densities)
%Sept2021: Updated code for Caro's paper (only relevant strategies here)
%Only site-trigger for GBR, Central, Central-South
%Applying new Habitat maps and re-defined GBR regions

%%Load the new regions
load('New_regions.mat')
new_regions=cell2table(nregions, 'VariableNames', {'Index' 'Priority' 'region'});
%criteria used

criteria=[];
others=find(isnan(META.COTS_control_sites(:,2))==1); %% This is the index of reefs outside the MP
idx=META.cntrl_reefID; %% Reefs within MPA

switch META.COTS_reefs2cull_strat
    case 14 %% No regional strategy: GBR-wide
        % YM: tentative changes below to initialise randomised priority and
        % non-priority lists and keep them unchanged for the rest of the  simulation (Caro's comparison with CoCoNet).
        % Requires turning 'COTS fixed list' from 0 to 1 in parsList2.csv to avoid any further randomisation (thisboatorder below).
        
        %change which file is read below to change focal region for different strategies
        if t==(META.COTS_control_start+1)%load only once to save time
            csirolist=load('CSIROprioritylist.mat');
            priority=csirolist.CSIROlist;
%             priority = csirolist.CSIROlist(randperm(length(csirolist.CSIROlist))) ;%YM: randomize the priority list from start.
            others=find(isnan(META.COTS_control_sites(:,2))==1); %% Updated--This is the index of reefs outside the MP
            RESULT.COTS_priority_list(1).list=priority;
            
            nonpriority=transpose(1:3806);%change this as needed if we only consider a certain region, i.e. then we would not create non-priority list from all 3806 GBR reefs
%             nonpriority = randperm(META.nb_reefs)';% YM: randomize the list of all simulated reefs from start
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
        COTS_ctowd=COTS_current_tow_density(find(idx==META.reef_ID)); %%Updated--only controlled reefs  %% YM: corrected for sizeable code
        current_pref_coral=squeeze(sum(RESULT.coral_pct2D(:,t,META.COTS_pref_corals),3));%coral cover of preferred groups, currently 1 to 4
        current_pcoral=current_pref_coral(find(idx==META.reef_ID)); %%Updated--only controlled reefs  %% YM: corrected for sizeable code
        
        sorted_COTS_density=zeros(size(sorted_indices_all,1),1);%see below what each variable stores
        sorted_ET_density=zeros(size(sorted_indices_all,1),1);
        sorted_pref_coral=zeros(size(sorted_indices_all,1),1);
        
        
        %%Updated based on the new reef list
        for i=1:size(sorted_indices_all,1)
            id=find(idx(:,1)==sorted_indices_all(i,1));%this assigns the appropiate values to the permuted list
            
            if isempty(id)==0 %% YM: corrected for sizeable code
                
                sorted_COTS_density(i,1)=COTS_ctowd(id);
                sorted_ET_density(i,1)=current_COTS_ET(id);
                sorted_pref_coral(i,1)=current_pcoral(id);
            end
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
        criteria1(find(idx==META.reef_ID),5)=COTS_ctd;%unsorted COTS ET density per tow per reef to check that the sorting was not broken %% YM: corrected for sizeable code
        criteria1(find(idx==META.reef_ID),6)=current_COTS_ET; %% YM: corrected for sizeable code
        cprefc=current_pref_coral;
        cprefc(others)=[];%% remove non-controlled reefs
        criteria1(find(idx==META.reef_ID),7)=cprefc; %% YM: corrected for sizeable code
        
        global_trigger=0;
        criteria.criteria=criteria1;
        %end of case 14 here
        
        
    case 18  % Strategy to visit only  Central: With new Regions, this
        % is the North and Central sectors
        
        if t==(META.COTS_control_start+1)%load only once to save time
            north=strcmp(new_regions.region,'N');
            north=find(north);%%Finds all non-zero (i.e. yes N)
            central=strcmp(new_regions.region,'C');
            central=find(central);%%Finds all non-zero (i.e. yes C)
            Central=cat(1,north,central);%% All reefs for 'Central' Scenario
            priority=strcmp(new_regions.Priority,'Y');
            priority=find(priority);
            nonpriority=strcmp(new_regions.Priority, 'N');
            nonpriority=find(nonpriority);
            
            others=find(isnan(META.COTS_control_sites(:,2))==1); %% This is the index of reefs outside the MP
            nonpriority=setdiff(nonpriority,others); %%remove reefs outside MPA-- control doesnt operate there
            %%Specific for Central Scenario
            priority=intersect(Central,priority);
            nonpriority=intersect(Central,nonpriority);
            
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
        COTS_ctowd=COTS_current_tow_density(find(idx==META.reef_ID)); %%Updated--only controlled reefs %% YM: corrected for sizeable code
        current_pref_coral=squeeze(sum(RESULT.coral_pct2D(:,t,META.COTS_pref_corals),3));%coral cover of preferred groups, currently 1 to 4
        current_pcoral=current_pref_coral(find(idx==META.reef_ID)); %%Updated--only controlled reefs %% YM: corrected for sizeable code
        sorted_COTS_density=zeros(size(sorted_indices_all,1),1);%see below what each variable stores
        sorted_ET_density=zeros(size(sorted_indices_all,1),1);
        sorted_pref_coral=zeros(size(sorted_indices_all,1),1);
        
        %%Updated based on the new reef list
        for i=1:size(sorted_indices_all,1)
            
            if isempty(id)==0 %% YM: corrected for sizeable code
                
                id=find(idx(:,1)==sorted_indices_all(i,1));%this assigns the appropiate values to the permuted list
                sorted_COTS_density(i,1)=COTS_ctowd(id);
                sorted_ET_density(i,1)=current_COTS_ET(id);
                sorted_pref_coral(i,1)=current_pcoral(id);
                
            end
        end
        
        criteria1=zeros(size(sorted_indices_all,1),8);%this stores everythign we need so it is onyl calcualted once
        sorted_indices=sorted_indices_all;%lazy way to not change everythign else in here
        criteria1(:,1)=sorted_indices_all;%all reefs after perumtation, first priority list then nonpriority
        criteria1(:,2)=sorted_COTS_density;%COTS density per tow that needs to correspond to list of reefs in 2nd column
        criteria1(:,3)=sorted_ET_density;%COTS ET density per tow that needs to correspond to list of reefs in 2nd column
        criteria1(:,4)=sorted_pref_coral;%current coral cover that need sto correspond to list of reefs in 2nd column
        COTS_ctd=COTS_current_tow_density;
        COTS_ctd(others)=[];%% remove non-controlled reefs
        criteria1(find(idx==META.reef_ID),5)=COTS_ctd;%unsorted COTS ET density per tow per reef to check that the sorting was not broken %% YM: corrected for sizeable code
        criteria1(find(idx==META.reef_ID),6)=current_COTS_ET; %% YM: corrected for sizeable code
        
        cprefc=current_pref_coral;
        cprefc(others)=[];%% remove non-controlled reefs
        criteria1(find(idx==META.reef_ID),7)=cprefc; %% YM: corrected for sizeable code
        
        global_trigger=0;
        criteria.criteria=criteria1;
        
    case 19  % Strategy to visit only  Central & South. With new Regions, this
        % are ALL except FN
        if t==(META.COTS_control_start+1)%load only once to save time
            
            Farnorth=strcmp(new_regions.region,'FN');
            Farnorth=find(Farnorth);%%Finds all non-zero (i.e. yes FN)
            otherreefs=transpose(1:3806);
            
            Centralsouth=setdiff(otherreefs,Farnorth);%% All reefs for 'Central-South' Scenario
            priority=strcmp(new_regions.Priority,'Y');
            priority=find(priority);
            nonpriority=strcmp(new_regions.Priority, 'N');
            nonpriority=find(nonpriority);
            others=find(isnan(META.COTS_control_sites(:,2))==1); %% This is the index of reefs outside the MP
            nonpriority=setdiff(nonpriority,others); %%remove reefs outside MPA-- control doesnt operate there
            %%Specific for Central-South Scenario
            priority=intersect(Centralsouth,priority);
            nonpriority=intersect(Centralsouth,nonpriority);
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
        COTS_ctowd=COTS_current_tow_density(find(idx==META.reef_ID)); %%Updated--only controlled reefs %% YM: corrected for sizeable code
        current_pref_coral=squeeze(sum(RESULT.coral_pct2D(:,t,META.COTS_pref_corals),3));%coral cover of preferred groups, currently 1 to 4
        current_pcoral=current_pref_coral(find(idx==META.reef_ID)); %%Updated--only controlled reefs %% YM: corrected for sizeable code
        
        sorted_COTS_density=zeros(size(sorted_indices_all,1),1);%see below what each varaible stores
        sorted_ET_density=zeros(size(sorted_indices_all,1),1);
        sorted_pref_coral=zeros(size(sorted_indices_all,1),1);
        
        %%Updated based on the new reef list
        for i=1:size(sorted_indices_all,1)
            
            if isempty(id)==0 %% YM: corrected for sizeable code
                
                id=find(idx(:,1)==sorted_indices_all(i,1));%this assigns the appropiate values to the permuted list
                sorted_COTS_density(i,1)=COTS_ctowd(id);
                sorted_ET_density(i,1)=current_COTS_ET(id);
                sorted_pref_coral(i,1)=current_pcoral(id);
                
            end
        end
        
        criteria1=zeros(size(sorted_indices_all,1),8);%this stores everythign we need so it is onyl calcualted once
        sorted_indices=sorted_indices_all;%lazy way to not change everythign else in here
        criteria1(:,1)=sorted_indices_all;%all reefs after perumtation, first priority list then nonpriority
        criteria1(:,2)=sorted_COTS_density;%COTS density per tow that need sto correspond to list of reefs in 2nd column
        criteria1(:,3)=sorted_ET_density;%COTS ET density per tow that need sto correspond to list of reefs in 2nd column
        criteria1(:,4)=sorted_pref_coral;%current coral cover that need sto correspond to list of reefs in 2nd column
        COTS_ctd=COTS_current_tow_density;
        COTS_ctd(others)=[];%% remove non-controlled reefs
        
        criteria1(find(idx==META.reef_ID),5)=COTS_ctd;%unsorted COTS ET density per tow per reef to check that the sorting was not broken %% YM: corrected for sizeable code
        criteria1(find(idx==META.reef_ID),6)=current_COTS_ET; %% YM: corrected for sizeable code
        
        cprefc=current_pref_coral;
        cprefc(others)=[];%% remove non-controlled reefs
        criteria1(find(idx==META.reef_ID),7)=cprefc; %% YM: corrected for sizeable code
        
        global_trigger=0;
        criteria.criteria=criteria1;
end