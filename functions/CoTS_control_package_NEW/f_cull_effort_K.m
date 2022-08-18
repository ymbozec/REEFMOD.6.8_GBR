function [ META ] = f_cull_effort_K( META )
%F_CULL_EFFORT Summary of this function goes here
%   Detailed explanation goes here

boats = length(META.boatProperties.boatID);


%if META.COTS_boat_cull_effort==0%if total cull effort in km2 is not pre-specified, use calculations from number of days
    for i=1:boats %for each boat
        %average distance of 20 min swim in m; calculated from timed swims, assume AMPTO moves at same speed, although culls probably slower
        swimd=480;
        %average width covered during swim in m; manual reach; no changes due to habitat complexity etc, no slowdown due to high densities
        swimw=3;
        %average area covered by a diver per hour in m2; 4320m2, or 66x66m
        swima=swimd*swimw*3;
        %number of divers per boat; AMPTO
        divers=META.boatProperties.divers(i);
        %number of hours dived per dive; AMPTO
        divet=2/3;%was 2/3 originally, maybe still is?
        %number of dives per diver per day; AMPTO; 2.67h per diver, 26.67h per boat
        diven=4;
        %area covered by boat per day; 115200m2, or 340x340m, or 0.1152km2
        area_day=swima*divers*divet*diven;
        META.boatProperties.areaPerDay(i) = area_day;
        %boat days at sea; AMPTO 220-240/year and some lost due to weather, so 100 days/6months; some also lost to travel, but not for now
        boat_days=META.boatProperties.boatDays(i);
        %effort quota, or area cleaned by boat per 6 months; 11.52km2 per 6 months
        META.boatProperties.effortQuota(i) = area_day*boat_days;
        META.boatProperties.totalTeamDives(i)=META.boatProperties.boatDays(i)*diven;%total dives a team can make; each should be on a new site
        META.boatProperties.totalInidvDives(i)=META.boatProperties.boatDays(i)*diven*META.boatProperties.divers(i);%tdives of individual divers; note that this assumes divers fromt eh same boat can all be on different sites on the same reef durign the same dive whcih si probably unrealistic
        
    end   
end

%end

