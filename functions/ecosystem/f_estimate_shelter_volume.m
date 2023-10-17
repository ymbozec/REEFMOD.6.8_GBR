% Estimate total shelter volume provided by all colonies of a given species
%
% Yves-Marie Bozec, y.bozec@uq.edu.au, 03/2023
%__________________________________________________________________________


function [total_shelter_volume_per_taxa] = f_estimate_shelter_volume(living_planar_coral_cover_cm2, SV_a, SV_b, SV_rand)

% SV_a and SV_b are the slope and intercept of the log-log relationship between colony area 
% and shelter volume (in dm3) provided by that colony (following Urbina-Barreto et al, 2021):
% log(S) = b + a*log(PA)
% SV parameters are specific to coral morphology (see Urbina-Barreto et al 2021) and defined in PARAMETERS.m.
% The function sums all colony-based shelter volumes to give an estimate of the total shelter volume (dm3)
% provided by the coral group for the reef grid.
% SV_rand adds uncertainty in the prediction (95% prediction intervals, calculated by Ryan Heneghan from the published data)

% living_planar_coral_cover_cm2 = vector of all colony 2D areas of the coral group

if living_planar_coral_cover_cm2 ~= 0
    
    S = exp(SV_a*log(living_planar_coral_cover_cm2) + SV_b + normrnd(0, SV_rand, length(living_planar_coral_cover_cm2),1));    
    total_shelter_volume_per_taxa = sum(S);
    
else
    total_shelter_volume_per_taxa = 0;
end