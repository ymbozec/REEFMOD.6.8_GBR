% Determine the balance amount of grazing among available algae
%
% Yves-Marie Bozec, y.bozec@uq.edu.au, created Nov 2011
%
% 20/09/2016: This is a new version, with code optimized and new rules of algal consumption:
% Preference order for switching is in ALGAL.feeding_prefs
%__________________________________________________________________________

function ALGALREMOVAL2 = f_algal_removal(algal_cm2, ALGALREMOVAL, algal_prefs, total_area_cm2)


total_algal_cm2 = sum(algal_cm2,1) ; % total cover of each algae over the whole grid

ALGALinit = (total_algal_cm2/total_area_cm2)';
algal_balance = ALGALinit - ALGALREMOVAL ; % >0 excess, <0 shortfall

shortfall = find(algal_balance<0);
excess = find(algal_balance>0);

excess(excess==3)=[]; % Pessimistic scenario whereby LOB is never eaten to compensate shortfalls in other algae

if isempty(shortfall)==0 && isempty(excess)==0
    
    for i = 1:length(shortfall) %1:length(algal_balance)  % check balance for each alga
        
        a = shortfall(i) ; % a designates the algae in shortfall
        [ordered_pref,I] = sort(algal_prefs(excess),'ascend');
        
        for k = 1:length(excess)
        
%         while algal_balance(a)<0 && isempty(excess)==0
%             k=k+1 ;
            available = algal_balance(excess(I(k))) + algal_balance(a) ;
            
            if available >= 0 % the other algae still in excess
                algal_balance(a) = 0 ; % the shorfall of a is offset
                algal_balance(excess(I(k))) = available ;
                
            else % shortfall of the other algae                
                algal_balance(a) = algal_balance(a) - algal_balance(excess(I(k))) ;
                algal_balance(excess(I(k))) = 0 ;
            end
            
            if algal_balance(a)>=0
                break
            end
        end
    end
end
algal_balance(algal_balance<0) = 0; % TEST

ALGALREMOVAL2 = (ALGALinit - algal_balance)';
