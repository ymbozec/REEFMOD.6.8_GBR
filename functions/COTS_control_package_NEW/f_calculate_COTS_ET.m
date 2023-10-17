function [ET_COTStow]= f_calculate_COTS_ET(RESULT,t, META)

% fucntion that checks whether cover of preffered coral groups is above ecological threshold (ET) for
% COTS control

current_pref_coral=squeeze(sum(RESULT.coral_pct2D(:,t,META.COTS_pref_corals),3));%coral cover of preferred groups, currently 1 to 4

%%OLD equations used in NESP report
%ET_COTSha=0.015*((20.*current_pref_coral)+4);% coral cover-dependent ecological threshold per Condie et al:
%ET_COTStow=0.015*ET_COTSha;%Used in convert the ecological threshold to COTS per tow which is originally from Moran and De'ath (see below)

%% ET = alpha * (20*(C/K) + 4). (alpha=0.015 is convertion factor as per condie to convert to manta tow calculated from Moran and De'ath 92) 
%% K is the maximum area of available habitat and C is the cover of corals preferred by CoTS, both expressed in the same arbitrary units. 
%% C/K is therefore the proportional cover of corals preferred by CoTS (for Reef MOd using percentage should be (C/(K/100))

%%REVISITED for paper after input from Scott
ET_COTStow=0.015*((20.*current_pref_coral/100)+4);% 

end

%%Note that Ecological Threshold (0.04 COTS / minute bottom time CPUE, is 
%%equivalent to 0.075 COTS per manta tow in the model.

