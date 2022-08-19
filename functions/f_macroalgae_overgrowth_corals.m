%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Y.-M. Bozec, MSEL, created Jun 2013.
%
% Implements the overgrowth of corals by any macroalgae
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% a = 3 -> Lobophora
% a = 2 -> Dictyota
% local_Malg_pct is the current cover of macroalgae over the 5 cells
% the function operates for all the colonies of every species within a cell (vectorization)

% // in the refrence code this is implemented as a series of code blocks one for each coral type
% // here each coral is processed by the same code, since the only difference appears to be in
% // the value returned by CORAL.lobophora_reduce_rate which formerly was ALGAL_CORAL_BRD_RESIST,
% // ALGAL_CORAL_SPW_RESIST, ALGAL_CORAL_BRD_VUL, ALGAL_CORAL_SPW_VUL and effectively zero for Acropora
% // since this did not have a code block. Now instead of having explicit code for each coral type just
% // set the value of CORAL.lobophora_reduce_rate for each coral type (zero for no effect, eg, acropora)

function [coral_cm2, algal_cm2] = ...
    f_macroalgae_overgrowth_corals(coral_cm2, algal_cm2, species_ID, a, alg_env_prop, reduce_rate, nb_coral_types)

id1=spones(coral_cm2); % faster than ones, here
id1(coral_cm2<0)=0 ; % exclude dead colonies (negatives)

[ncell,ncol]=size(coral_cm2); % ncell = number of cells, ncol = number of colonies per species

alg_env_prop(alg_env_prop<0.01)=0; % ignore algal environment with less than 1% algal cover

perim = 2*pi*sqrt(coral_cm2.*id1/pi) ;
perim_alg = alg_env_prop(:,ones(1,ncol)).* perim ;  % arc length (cm) in contact with alga a

reduce_rates = zeros(ncell,ncol);

col_start = 1;
col_stop = 0;

for s = 1:nb_coral_types
    
    col_stop = col_stop + species_ID(s) ; 
    reduce_rates(:,col_start:col_stop)=reduce_rate(s) ;
    col_start = col_start + species_ID(s) ;
    
end

area_lost = reduce_rates.*perim_alg ; % loss of coral in cm2 (sparse matrix
area_lost = floor(area_lost) ;

area_lost(area_lost > coral_cm2.*id1) = coral_cm2(area_lost > coral_cm2.*id1) ; % cannot loose more than the colony area

coral_cm2 = coral_cm2 - area_lost ; % update coral cover according to area_lost
algal_cm2(:,a) = algal_cm2(:,a) + sum(area_lost,2) ;