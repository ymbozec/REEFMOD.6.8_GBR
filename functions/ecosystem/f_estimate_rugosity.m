function [rugosity, SI_reef] = f_estimate_rugosity(REEF, META, CORAL, coral_sp_pct2D)

% faire d'abord un test sans les ratio species-specific
SI_reef= sum(REEF.floor_SA_cm2)/META.total_area_cm2;
rugosity = 0.88*SI_reef;

% sinon:
% LiveCoralRugosity = sum(CORAL.CI_SI_ratio.*coral_sp_pct2D)/100 
% DeadCoralRugosity = 0.88*(100-sum(coral_sp_pct2D))/100
% rugosity = SI_reef*(LiveCoralRugosity + DeadCoralRugosity)