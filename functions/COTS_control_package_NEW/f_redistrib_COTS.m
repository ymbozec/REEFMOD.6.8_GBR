function [COTS_sites_density, stored_betarndN]=f_redistrib_COTS(current_COTS_forControl, META)


% function that redistributes reef-level COTS population to individual
% control sites using beta infalted random numbers at the moment
mu = 0.3324884;
sigma = 0.2749365;
nu = 2.717073;
tau = 0.02195122;

%COTS_sites_density=cell(META.nb_reefs,1);
%stored_betarndN=cell(META.nb_reefs,1);

% for i=1:3806
%     if META.COTS_control_sites(i,2)==1
%         COTS_sites_density(i,1)={current_COTS(i,:)};
%     else
%         check_zero=0;
%         while check_zero==0%to avoid all random numbers beign zero and then not redistributing anything
%             %betarndN=betainfrnd(META.COTS_control_sites(i,2),mu, sigma, nu, tau);
%             tot=1;
%             Nsites=META.COTS_control_sites(i,2);
%             betarndN=nan(1,Nsites);
%             for s=1:(Nsites-1)
%                 r=betainfrnd(1, mu, sigma, nu, tau);%%Using random-intercept model
%                 betarndN(1,s)=tot*r;
%                 tot=tot-betarndN(1,s);
%                 if tot<0%if for some reason this becomes negative, break
%                     betarndN(1,s:(Nsites-1))=0;
%                     break;
%                 end
%             end
%             betarndN(1,Nsites) = 1-nansum(betarndN(1,:));
%             if sum(betarndN)>0%rpeat if all rndaon numbers zero
%                 check_zero=1;
%             end
%         end
        
            %now use the random numbers to redistribute reef-level COTS population; all
            %size classes are redistributed usign the same number; would be interesting to know whether outbreaks have different size class distribution from background populations
            %this_reef_sites=zeros(META.COTS_control_sites(i,2),META.COTS_maximum_age);
            
%             for j=1:META.COTS_control_sites(i,2)
%                 this_reef_cots_4redistrib=current_COTS(i,:)*META.COTS_control_sites(i,2);%note that this redistributes the same size class distrubtion to all sites and in all circumstances
%                 this_reef_sites(j,:)=this_reef_cots_4redistrib.*betarndN(j);%used to be current_COTS(i,:) but this just divided the reef-level COTS to all sites
%             end
%    COTS_sites_density(i,1)={this_reef_sites};
%             stored_betarndN(i,1)={betarndN};
%         end
%     end
%     
% end



%%A new list of reefs has been created to exclude reefs outside the Marine
%%Park for control, so remove those reefs and loop only over reefs within
%%the MP

%Newlist=[META.reef_ID' META.COTS_control_sites(:,2)]; %%Add index to track other metrics
%Newlist=rmmissing(Newlist); %%Remove missing reefs
% META.cntrl_reefID=[Newlist(:,1)];%% get the reefID for reefs within MPA


% others=find(isnan(META.COTS_control_sites(:,2))==1); %% This is the position of the nan values
% META.COTS_control_sites(isnan(META.COTS_control_sites))=0;

nbreefs = length(META.cntrl_sites);

COTS_sites_density=cell([nbreefs,1]);
stored_betarndN=cell([nbreefs,1]);

for i=1:length(META.cntrl_sites)
    %idx=[Newlist(i,1)];%% get the reefID
    if META.cntrl_sites(i,2)==1
        %COTS_sites_density(i,1)={current_COTS(idx,:)};
        COTS_sites_density(i,1)={current_COTS_forControl(i,:)};
    %elseif isnan(Newlist(i,2))
     %   COTS_sites_density(i,1)={current_COTS(idx,:)};
    else
        check_zero=0;
        while check_zero==0%to avoid all random numbers beign zero and then not redistributing anything
            tot=1;
            Nsites=META.cntrl_sites(i,2);
            betarndN=nan(1,Nsites);
            for s=1:(Nsites-1)
                r=betainfrnd(1, mu, sigma, nu, tau);%%Using random-intercept model
                betarndN(1,s)=tot*r;
                tot=tot-betarndN(1,s);
                if tot<0%if for some reason this becomes negative, break
                    betarndN(1,s:(Nsites-1))=0;
                    break;
                end
            end
            betarndN(1,Nsites) = 1-nansum(betarndN(1,:));
            
            if sum(betarndN)>0%repeat if all random numbers zero
                check_zero=1;
            end
        end
        
        %now use the random numbers to redistribute reef-level COTS population; all
        %size classes are redistributed usign the same number; would be interesting to know whether outbreaks have different size class distribution from background populations
        %this_reef_sites=zeros(META.COTS_control_sites(i,2),META.COTS_maximum_age);
        this_reef_sites=zeros(META.cntrl_sites(i,2),META.COTS_maximum_age);
        
        for j=1:META.cntrl_sites(i,2)
            
            %this_reef_cots_4redistrib=current_COTS(idx,:)*Newlist(i,2);%note that this redistributes the same size class distrubtion to all sites and in all circumstances
            this_reef_cots_4redistrib=current_COTS_forControl(i,:)*META.cntrl_sites(i,2);%note that this redistributes the same size class distrubtion to all sites and in all circumstances
            this_reef_sites(j,:)=this_reef_cots_4redistrib.*betarndN(j);%used to be current_COTS(i,:) but this just divided the reef-level COTS to all sites
        end
        COTS_sites_density(i,1)={this_reef_sites};
        stored_betarndN(i,1)={betarndN};
    end
end

end