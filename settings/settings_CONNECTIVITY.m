% REEFMOD-GBR connectivity settings
%
% Yves-Marie Bozec, y.bozec@uq.edu.au, 09/2019
%__________________________________________________________________________

%%%% Coral (Acropora) and COTS connectivity for the GBR (Hock et al. 2017)
load('GBR_CONNECT_7years.mat'); % Connectivity matrices for Acropora and CoTS (7 spawning seasons for each)
% We've got connectivity matrices for summers 2008-2009, 2010-11, 2011-12, 2012-13, 2014-15, 2015-16, 2016-2017
% Corresponding Reefmod summer-year names are: 2009, 2011, 2012, 2013, 2015, 2016, 2017 (i.e., 2014 is missing)
% Need to exclude 2009 and 2017 to match WQ layers + needs to fill matrix

id1 = randi(7); % randomly select one of the 7 matrices for 2013-14
id2 = randi(7); % randomly select one of the 7 matrices for 2017-18

self1 = sparse(diag(META.coral_min_selfseed * ones(1,META.nb_reefs)));

CONNECT(1).ACROPORA = self1 + GBR_CONNECT(2).ACROPORA(META.reef_ID,META.reef_ID) ; % 2010-11
CONNECT(2).ACROPORA = self1 + GBR_CONNECT(3).ACROPORA(META.reef_ID,META.reef_ID) ; % 2011-12
CONNECT(3).ACROPORA = self1 + GBR_CONNECT(4).ACROPORA(META.reef_ID,META.reef_ID) ; % 2012-13
CONNECT(4).ACROPORA = self1 + GBR_CONNECT(id1).ACROPORA(META.reef_ID,META.reef_ID) ; % 2013-14 fill with 2010-11 according to Karlo
CONNECT(5).ACROPORA = self1 + GBR_CONNECT(5).ACROPORA(META.reef_ID,META.reef_ID) ; % 2014-15
CONNECT(6).ACROPORA = self1 + GBR_CONNECT(6).ACROPORA(META.reef_ID,META.reef_ID) ; % 2015-16
CONNECT(7).ACROPORA = self1 + GBR_CONNECT(7).ACROPORA(META.reef_ID,META.reef_ID) ; % 2016-17
CONNECT(8).ACROPORA = self1 + GBR_CONNECT(id2).ACROPORA(META.reef_ID,META.reef_ID) ; % 2017-18


if META.doing_COTS == 1
    
    self2 = sparse(diag(META.COTS_min_selfseed * ones(1,META.nb_reefs,1)));
    
    CONNECT(1).COTS = self2 + GBR_CONNECT(2).COTS(META.reef_ID,META.reef_ID) ;
    CONNECT(2).COTS = self2 + GBR_CONNECT(3).COTS(META.reef_ID,META.reef_ID) ;
    CONNECT(3).COTS = self2 + GBR_CONNECT(4).COTS(META.reef_ID,META.reef_ID) ;
    CONNECT(4).COTS = self2 + GBR_CONNECT(id1).COTS(META.reef_ID,META.reef_ID) ; % fill 2014 with 2011 according to Karlo
    CONNECT(5).COTS = self2 + GBR_CONNECT(5).COTS(META.reef_ID,META.reef_ID) ;
    CONNECT(6).COTS = self2 + GBR_CONNECT(6).COTS(META.reef_ID,META.reef_ID) ;
    CONNECT(7).COTS = self2 + GBR_CONNECT(7).COTS(META.reef_ID,META.reef_ID) ;
    CONNECT(8).COTS = self2 + GBR_CONNECT(id2).COTS(META.reef_ID,META.reef_ID) ;
    
end


%% Jan 2022: for all reefs of the focal region, determine in-degree from reefs outside of that region
% Returns a vector (length = nb of inside reefs) of summed in-degree probabilitites of larval dispersal from outside reef
if isempty(META.outside_reef_ID)==0
    
    CONNECT(1).ACROPORA_ext = sum(GBR_CONNECT(2).ACROPORA(META.outside_reef_ID,META.reef_ID),1) ; % 2010-11
    CONNECT(2).ACROPORA_ext = sum(GBR_CONNECT(3).ACROPORA(META.outside_reef_ID,META.reef_ID),1) ; % 2011-12
    CONNECT(3).ACROPORA_ext = sum(GBR_CONNECT(4).ACROPORA(META.outside_reef_ID,META.reef_ID),1) ; % 2012-13
    CONNECT(4).ACROPORA_ext = sum(GBR_CONNECT(id1).ACROPORA(META.outside_reef_ID,META.reef_ID),1) ; % 2013-14 fill with 2010-11 according to Karlo
    CONNECT(5).ACROPORA_ext = sum(GBR_CONNECT(5).ACROPORA(META.outside_reef_ID,META.reef_ID),1) ; % 2014-15
    CONNECT(6).ACROPORA_ext = sum(GBR_CONNECT(6).ACROPORA(META.outside_reef_ID,META.reef_ID),1) ; % 2015-16
    CONNECT(7).ACROPORA_ext = sum(GBR_CONNECT(7).ACROPORA(META.outside_reef_ID,META.reef_ID),1) ; % 2016-17
    CONNECT(8).ACROPORA_ext = sum(GBR_CONNECT(id2).ACROPORA(META.outside_reef_ID,META.reef_ID),1) ; % 2017-18

%% We don't do it for COTS due to their very patchy distribution
%     if META.doing_COTS == 1
%         
%         CONNECT(1).COTS_ext = sum(GBR_CONNECT(2).COTS(META.outside_reef_ID,META.reef_ID),1) ; % 2010-11
%         CONNECT(2).COTS_ext = sum(GBR_CONNECT(3).COTS(META.outside_reef_ID,META.reef_ID),1) ; % 2011-12
%         CONNECT(3).COTS_ext = sum(GBR_CONNECT(4).COTS(META.outside_reef_ID,META.reef_ID),1) ; % 2012-13
%         CONNECT(4).COTS_ext = sum(GBR_CONNECT(id1).COTS(META.outside_reef_ID,META.reef_ID),1) ; % 2013-14 fill with 2010-11 according to Karlo
%         CONNECT(5).COTS_ext = sum(GBR_CONNECT(5).COTS(META.outside_reef_ID,META.reef_ID),1) ; % 2014-15
%         CONNECT(6).COTS_ext = sum(GBR_CONNECT(6).COTS(META.outside_reef_ID,META.reef_ID),1) ; % 2015-16
%         CONNECT(7).COTS_ext = sum(GBR_CONNECT(7).COTS(META.outside_reef_ID,META.reef_ID),1) ; % 2016-17
%         CONNECT(8).COTS_ext = sum(GBR_CONNECT(id2).COTS(META.outside_reef_ID,META.reef_ID),1) ; % 2017-18
%     end
end

clear GBR_CONNECT self1 self2

