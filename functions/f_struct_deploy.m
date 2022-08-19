% -------------------------------------------------------------------------
% Y.-M. Bozec, MSEL, created Aug 2015.
% For optimization
% -------------------------------------------------------------------------

function [coral_cm2, surface_cm2, volume_cm3, clade, colony_ID, species_ID] = f_struct_deploy(coral)

% This extract coral matrices from the 1 x s structure array 'coral' 
% This is necessary because operations are not possible (?) on structure
% arrays and we want to avoid looping on coral species when operating each
% process

nb_coral_types = size(coral,2) ;
% Extract cover of every species within a sparse matrix
coral_cm2 = [coral.cover_cm2];

% Do the same for:
surface_cm2 = [coral.surface_cm2];
volume_cm3 = [coral.volume_cm3];
colony_ID = [coral.colony_ID];
clade = [coral.clade];

% Finally build a matrix storing the dmensions of each coral matrix
% This allows tracking the id of each species in coral_cm2
species_ID = zeros(1,nb_coral_types);

for s = 1:nb_coral_types
    
    species_ID(1,s)=size(coral(s).cover_cm2,2);
    
end