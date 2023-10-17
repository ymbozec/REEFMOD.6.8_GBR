function [ RESULT ] = f_COTS_control_CSIRO( META, RESULT, t, REEF_POP, CONNECT)
%F_COTS_CONTROL Summary of this function goes here

%COTS control function that uses CSIRO-like control effort and triggers

if t==META.COTS_control_start%set this to zero if this is the first time COTS control is run
    RESULT.last_reef_COTScontrolled=0;%no reefs controlled before, so set this to zero for makeReefList file
    %META.COTS_reefs2cull_strat;%enable to check if it is using the right strategy by writing it out in console...
    %     csirolist=load('CSIROprioritylist.mat');%this is now in makeReefList to enable different strategies to read different files
    %     RESULT.COTS_priority_list=csirolist.CSIROlist;
    %META.COTS_postcontrol_proportions=[ 0 0 (1-META.COTS_detectability(3:end))./(sum(1-META.COTS_detectability(3:end)))];
end

t=t+1;

% if t<META.nb_time_steps%if not the final time step, as the final time step creates problem with storing populations in future time steps...
if t<=META.nb_time_steps %otherwise stops controlling one step before the last step
    
    % keep track of: criteria, visited reef IDs, number of reefs visited, site distribution
    % before control, site distribution after control, dives per reef,
    COTS_records=struct('COTS_criteria',[],'control_record',[],'sites_initial',[],'sites_final',[],'COTS_initial',[],'COTS_final',[],'COTS_ET',[],'betarndN',[]);
    thisboatorder=META.boatProperties.fixedOrder(1);
    
    %calculate the full quota for timestep
    timestep_dives=0;
    
    for boat = 1:length(META.boatProperties.boatID)
        %calculate full timestep effort in terms of (team) dives
        %timestep_dives=timestep_dives+META.boatProperties.totalTeamDives(boat);
        timestep_dives=timestep_dives+META.boatProperties.totalInidvDives(boat);
    end
    
    %remaining_dives=timestep_dives;%full quota at the beginnign of a timestep; in the number of avaialble (team) dives
    
    %calculate the current number of COTS per tow; to be used in reef-level and season-level triggers
    current_COTS=reshape(RESULT.COTS_all_densities(:,t,1:META.COTS_maximum_age),META.nb_reefs,META.COTS_maximum_age);
    current_COTS_init = current_COTS ; %YM: added to record the number of CoTS killed
    %COTS_current_tow_density=sum(current_COTS(:,META.COTS_adult_min_age:end).*META.COTS_detectability(META.COTS_adult_min_age:end),2);
    %Newlist=[META.reef_ID' META.COTS_control_sites(:,2)]; %%Add index to track other metrics
    %Newlist=rmmissing(Newlist); %%Remove missing reefs
    idx=META.cntrl_reefID;
    
    
    COTS_records.COTS_initial=current_COTS(find(idx==META.reef_ID),:); %%Only track controlled reefs within MPA %% YM: corrected for sizeable code
    
    %calculate ecological threshold (ET) for each reef - dynamically for a current level of coral cover; used in both reef-level
    %and site-level calculations since coral tracked on a reef-level
    
    %current_COTS_ET=repmat(0.075,[3806 1]);
    
    current_COTS_ET=repmat(0.075,[size(idx) 1]);
    COTS_records.COTS_ET=current_COTS_ET;
    
    % redistribute reef-level COTS density to individual control sites; keep
    % track of random numbers used to check whether rng is working
    
    %Exclude Reefs outside MPark
    current_COTS_forControl=current_COTS(find(idx==META.reef_ID),:); %% YM: corrected for sizeable code
    [COTS_site_densities, stored_betarndN]=f_redistrib_COTS(current_COTS_forControl, META);
    
    COTS_records.sites_initial=COTS_site_densities;
    
    % make a list of reefs that determines the order in which they are visited;
    % only do this once per timestep, then go down the list checking the individual triggers until out of quota
    %[~, criteriaS, global_trigger, RESULT] = f_makeReefList(META, RESULT, t, 0, current_COTS, thisboatorder, current_COTS_ET, COTS_site_densities);
    [~, criteriaS, global_trigger, RESULT] = f_makeReefList_CCS(META, RESULT, t, 0, current_COTS, thisboatorder, current_COTS_ET, COTS_site_densities);
    
    
    
    global_control_trigger=1;%always do control, unless need to check global trigger; if checking global trigger, set control to global trigger (from f_makeReefList) and proceed
    if META.COTS_global_trigger==1
        global_control_trigger=global_trigger;
    end
    culled_reefs=0;%counter to keep track of how many reefs were controlled
    
    % keep track of: criteria, visited reef IDs, number of reefs visited, site distribution
    % before control, site distribution after control, dives per reef,
    control_records=struct('reef',[],'sites_dives',[],'size_classes_pre_control',[],'size_classes_post_control',[],'COTS_post_control_densities',[],'COTS_post_control_proportions',[],'COTS_post_control_overall_ET',[],'sites_over_ET',[],'this_sites_over_ET',[]);
    cnt1=1;
    if global_control_trigger==1
        %how to distribute effort using potentially different priority
        %lists; if there is onyl one list, i.e. effort is evenly dsitrubted
        %across focal regions, then this will be just a single list and the
        %for loop with efd iterator will run only once; as a result, criteria are now a
        %structure as well as priority and nonpriority lists
        effort_distribution=RESULT.COTS_effort_distribution;
        for efd=1:length(effort_distribution)
            criteria=criteriaS(efd).criteria;
            % YM: below code to resize the reef list and nly retain the simulated reef IDs
            I = ismember(criteria(:,1),META.reef_ID);
            mycriteria = criteria.*I(:,ones(1,size(criteria,2)));
            criteria = mycriteria(find(sum(mycriteria,2)~=0),:);
            % end of YM's code
            current_reef=1;%we will start with the first reef on the list
            remaining_dives=timestep_dives*effort_distribution(efd);%full quota at the beginning of a timestep fro this region; in the number of avaialble (team) dives
            while remaining_dives>0 && current_reef<size(criteria,1) % YM April 10 2022 - add constraint for nb of simulated reefs
                %while there are dives remaining
                %general description: get first reef ID from the list
                %check if reef-level trigger is used, then use it
                %check if site-level trigger is valid at any sites
                %go only to those sites, recheck after every run to see if trigger
                %still valid, if not, increase the current_reef by one until out of quota
                this_reef=criteria(current_reef,1);
                %check if reef level trigger is used; if yes, check whether this reef satisifies the criterion; if not, advance to the next reef on list
                if META.COTS_reef_trigger==1
                    treat_this_reef=0;
                    while treat_this_reef==0
                        %if criteria(this_reef,2)<0.22
                        if criteria(current_reef,2)<0.22 %%Updated to match new reef list
                            current_reef=current_reef+1;
                            %if current_reef>META.nb_reefs%if we came to the end of the list, stop searching
                            if current_reef>length(criteria)
                                break;
                            end
                            %newid=find(criteria(:,1) == current_reef);
                            newid=find(criteria(:,1) == this_reef);
                            %this_reef=criteria(current_reef,1);
                            this_reef=criteria(newid,1);
                        else
                            treat_this_reef=1;
                        end
                    end
                end
                %if current_reef>META.nb_reefs%safeguard in case reef trigger iterated through all the reefs and none of them triggered, so that the loop does not simply proceed with the last value of this_reef
                if current_reef>length(criteria)
                    break;
                end
                %check whether this reefs is either on priority or nonpriority list, ignore if it is not
                %should now work with both case 14 and case 15 from makereeflist
                %if ismember(this_reef, RESULT.COTS_priority_list(efd).list) || ismember(this_reef, RESULT.COTS_nonpriority_list(efd).list)
                
                %this_reef_sites=COTS_site_densities{this_reef,1};
                newid=find(META.cntrl_reefID(:,1) == this_reef);%%Update to match new reef list
                
                if isempty(newid)==0 %% YM: corrected for sizeable code
                    
                    this_reef_sites=COTS_site_densities{newid,1};%%Update to match new reef list
                    this_sites_tow=zeros(size(this_reef_sites,1),1);
                    this_sites_over_ET=zeros(size(this_reef_sites,1),1);
                    for sts=1:size(this_reef_sites,1)%convert COTS at sites to COTS per tow, also check which ones are over ET
                        this_sites_tow(sts,1)=sum(this_reef_sites(sts,META.COTS_adult_min_age:end).*META.COTS_detectability(META.COTS_adult_min_age:end),2);
                        newid=find(criteria(:,1) == this_reef);
                        if this_sites_tow(sts,1)>criteria(newid,3)
                            this_sites_over_ET(sts,1)=1;
                        end
                    end
                    sites_over_ET=find(this_sites_over_ET);
                    if ~isempty(sites_over_ET) && remaining_dives>0%if there are sites to treat, and dives remaining, do control
                        control_records(cnt1).sites_over_ET=sites_over_ET;
                        control_records(cnt1).this_sites_over_ET=this_sites_over_ET;
                        RESULT.last_reef_COTScontrolled=this_reef;%set thsi reefs as the one that has been last controlled, to be used later; set here for every reef because we never know wwhen the control runs out of dives
                        record_site=zeros(length(sites_over_ET),2);%record which site on a reef was visited, and how many dives were perfomed on each site
                        %we need to rescale the COTS post-control
                        %proportions, or cpcp, to the fact that only
                        %every second size class has COTS in it in reefmod,
                        %and that they progress dynamically
                        this_cpcp=zeros(1,length(META.COTS_detectability));%use this as a temporary variable to calculate reduction later
                        this_cpcp(1,1:META.COTS_adult_min_age)=1;%control never affects the first two size classes
                        if rem(META.COTS_adult_min_age,2)==0%check if this is an even number, then start culls from approriate size class
                            if rem(t,2)~=0%if it is an odd time step (remember that in this file t=t+1 because we are modifying future populations), COTS will only be in even numbered columns
                                temp_cpcp=META.COTS_postcontrol_proportions(META.COTS_adult_min_age:2:end);%get even numbered cpcp, start from 4, since column 2 is not affected
                                temp_cpcp=temp_cpcp./(sum(temp_cpcp));%rescale it to only consider those columns; their relative proportions now add up to one
                                this_cpcp(META.COTS_adult_min_age:2:end)=temp_cpcp;%include it into
                            else%do the same thing, but now only consider even numbered columns
                                temp_cpcp=META.COTS_postcontrol_proportions((META.COTS_adult_min_age+1):2:end);
                                temp_cpcp=temp_cpcp./(sum(temp_cpcp));
                                this_cpcp((META.COTS_adult_min_age+1):2:end)=temp_cpcp;
                            end
                        else%else if odd number, do this
                            if rem(t,2)~=0%if it is an odd time step (remember that in this file t=t+1 because we are modifying future populations), COTS will only be in even numbered columns
                                temp_cpcp=META.COTS_postcontrol_proportions((META.COTS_adult_min_age+1):2:end);%get even numbered cpcp, start from 4, since column 2 is not affected
                                temp_cpcp=temp_cpcp./(sum(temp_cpcp));%rescale it to only consider those columns; their relative priortions now add up to one
                                this_cpcp((META.COTS_adult_min_age+1):2:end)=temp_cpcp;%include it into
                            else%do the same thing, but now only consider even numbered columns
                                temp_cpcp=META.COTS_postcontrol_proportions(META.COTS_adult_min_age:2:end);
                                temp_cpcp=temp_cpcp./(sum(temp_cpcp));
                                this_cpcp(META.COTS_adult_min_age:2:end)=temp_cpcp;
                            end
                        end
                        newid=find(criteria(:,1) == this_reef);
                        site_ETtow=criteria(newid,3);%the ecological threshold (ET) for this site/reef
                        
                        site=1;
                        record_site_pre_densities=zeros(length(sites_over_ET),16);
                        record_site_post_densities=zeros(length(sites_over_ET),16);
                        while site<=length(sites_over_ET) && remaining_dives>0%go throguh all sites
                            ctrl_site=sites_over_ET(site);
                            record_site(site,1)=ctrl_site;
                            this_site_tow=sum(this_reef_sites(ctrl_site,META.COTS_adult_min_age:end).*META.COTS_detectability(META.COTS_adult_min_age:end),2);
                            %bring it down to ET, but use detectability to reduce the classes
                            site_size_classes=this_reef_sites(ctrl_site,:);%list size classes
                            record_site_pre_densities(site,:)=this_reef_sites(ctrl_site,:);
                            for s=META.COTS_adult_min_age:length(site_size_classes)%reduce density per size class based on detectability so that it matches the relative abundance post control, which is clacualted from detectability
                                sizeclass_relative_abundance=site_size_classes(s);%current abundace of a COTS size class on this site
                                sizeclass_threshold_abundance=(site_ETtow)*this_cpcp(s);%what abundance should be if reef at ET, equation from COndie et al; this uses adjutesd cpcp, see above
                                if sizeclass_relative_abundance>sizeclass_threshold_abundance%if abundance greater than at what it should be at ET
                                    site_size_classes(s)=sizeclass_threshold_abundance;%reduce to abundance at ET
                                end
                            end
                            this_reef_sites(ctrl_site,:)=site_size_classes;
                            record_site_post_densities(site,:)=site_size_classes;
                            %here 8 is the standard number of divers, I use team dives per site, using different boats would require a rewrite
                            %also uses ceil to round, as the whoel team will finish the dive on the same site
                            control_dives=(ceil((4.18)*(this_site_tow/0.015)^0.667));
                            %control_dives=(ceil((4.18/META.boatProperties.divers(1))*(this_site_tow/0.015)^0.667));
                            remaining_dives=remaining_dives-control_dives;%subtract perfomed dives from total
                            if remaining_dives<0%break if no more dives left
                                remaining_dives=0;
                            end
                            newid=find(META.cntrl_reefID(:,1) == this_reef);%%Update to match new list
                            %COTS_site_densities(this_reef,1)={this_reef_sites};%update sites for this reef
                            COTS_site_densities(newid,1)={this_reef_sites};%update sites for this reef
                            record_site(site,2)=control_dives;
                            site=site+1;%up the iterator
                        end
                    
                    control_records(cnt1).reef=this_reef;%record reef
                    control_records(cnt1).sites_dives=record_site;%record sites and dives on sites
                    control_records(cnt1).size_classes_pre_control=record_site_pre_densities;%record size classes for each site that was controlled prior to actually doing the control; hopefully easier to keep track of pre and post control site densities; the rows here match with sites that were recorded in sites_dives
                    control_records(cnt1).size_classes_post_control=record_site_post_densities;%record size classes for each site that was controlled after the control; hopefully easier to keep track of pre and post control site densities; again, the rows here match with sites that were recorded in sites_dives
                    control_records(cnt1).COTS_post_control_densities=(site_ETtow)*this_cpcp;%record COTS_post_control_proportions i.e what the densities should maximally be for each size class when reef is at ET; if the site is treated, post-control density at specific size classes should be at (if it was higher) or below (stays the same) these values
                    control_records(cnt1).COTS_post_control_proportions=this_cpcp;%proportions of different size classes in terms of how they contribute to the overall ET; e.g. a value of 0.5 for a size class means that at maximum 50% of the density at ET should be in this size class; size classes below min adult age are not affected and are set to 1
                    control_records(cnt1).COTS_post_control_overall_ET=site_ETtow;%record overall ET for this reef/site; note that ET is decided at reef level the same as coral cover and is the same for all sites
                    cnt1=cnt1+1;
                    
                    end
                end
                %end
                current_reef=current_reef+1;%up the reef
                culled_reefs=culled_reefs+1;%up the count of controlled reefs
                %if current_reef>META.nb_reefs%if you run out of reefs, break
                if current_reef>length(criteria)
                    break;
                end
            end
        end
    end
    COTS_records.sites_final=COTS_site_densities;%record the final disturbtion of COTS across sites after control, to compare that control is working
    %now get reef-level COTS from individual, now treated, sites
    %for rfs=1:META.nb_reefs
    for rfs=1:length(META.cntrl_sites) %%Update to match new reef list
        recalc_COTS=zeros(1,META.COTS_maximum_age);
        temp_sites=COTS_site_densities{rfs,1};
        for sites=1:size(COTS_site_densities{rfs,1},1)
            recalc_COTS=recalc_COTS+temp_sites(sites,:);
        end
        %%Update to match new reef list
%         current_COTS_forControl(rfs,:)=recalc_COTS./META.COTS_control_sites(rfs,2);
        current_COTS_forControl(rfs,:)=recalc_COTS./size(COTS_site_densities{rfs,1},1); %% YM: corrected for sizeable code
        %current_COTS(rfs,:)=recalc_COTS./META.COTS_control_sites(rfs,2);%divide by number of sites to get reef average
        if isnan(current_COTS_forControl(rfs,:))
            current_COTS_forControl(rfs,:)=recalc_COTS;
        end
    end
    %%Update now for all other reefs
    current_COTS(find(idx==META.reef_ID),:)=current_COTS_forControl; %% YM: corrected for sizeable code
    COTS_records.COTS_final=current_COTS(find(idx==META.reef_ID),:);%Keep track only for controlled reefs %% YM: corrected for sizeable code
    
    
    %update things for the next time step
    RESULT.COTS_all_densities(:,t,1:META.COTS_maximum_age)=reshape(current_COTS,META.nb_reefs,1,META.COTS_maximum_age);
    RESULT.COTS_adult_densities(:,t)=sum(current_COTS(:,META.COTS_adult_min_age:end),2);
    % YM: commented code below as we just need to know which reefs where culled
%     RESULT.COTS_culled_reefs(1,t) = culled_reefs; %(YM: in fact they were the visited reefs, not necessarily culled)
     culled_reef_IDs = cat(1,control_records.reef);
     RESULT.COTS_culled_reefs(1:META.nb_reefs,t) = ismember(META.reef_ID,culled_reef_IDs);
    %YM: added to record the remaining control effort when simulating less then 3,806 reefs
    RESULT.COTS_control_remaining_dives(1,t) = remaining_dives;
    %YM: added to record the amount of CoTS killed
    RESULT.COTS_culled_density(1:META.nb_reefs,t) = round(sum(current_COTS_init(:,META.COTS_adult_min_age:end),2)...
        - sum(current_COTS(:,META.COTS_adult_min_age:end),2),4);

    %keep track of things that have been done
    COTS_records.COTS_criteria=criteriaS;
    COTS_records.control_record=control_records;
    COTS_records.betarndN=stored_betarndN;
    RESULT.COTS_records(1,t)=COTS_records;
    
    % NOTE YM: need to harmonize Karlo's t/t+1 with the storage in RESULT
    % 1) do we really need the above condition t<=META.nb_time_steps oter
    % control stops before the last time step
    % 2) RECORD of culling should be made at t-1 then
    
end%if not last time step

end %function