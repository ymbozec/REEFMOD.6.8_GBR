% Create demographic parameters of coral outplants as duplicate of their native counterparts
% (ie, native corals of the same morphological group)
% This is just a duplication, specific parameters can be refined in
% settings_RESTORATION.m
%
% Yves-Marie Bozec, y.bozec@uq.edu.au, 02/2021
%__________________________________________________________________________


function [CORAL,nb_coral_types] = f_create_outplant_demographics(CORAL, outplanted_species, nb_coral_types)


lengths = structfun(@(x) size(x,1), CORAL);

I=find(lengths==nb_coral_types);
fields = fieldnames(CORAL);
select_fields = fields(I);

for s = 1:length(outplanted_species)
    
    for f=1:length(select_fields)
        
        V = extractfield(CORAL,cell2mat(select_fields(f))) ;
        newV = [V  V(outplanted_species(s))];
        CORAL = setfield(CORAL,cell2mat(select_fields(f)),newV') ;
   % Note that prop_settlers (proportion of each coral type in maximum settlement) summed to 1 initially
   % Now above 1 but probably not a big deal as an addition of a new population?
    end
    
end

nb_coral_types = nb_coral_types + length(outplanted_species);

%% TEMP
% CORAL.BH_alpha(7:8) = 0