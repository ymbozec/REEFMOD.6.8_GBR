%  /////////////////////////////////////////////////////////////////////////////
%  // 3. FINALLY ADJUST THE TOTALS FOR CORAL GROWTH AND MORTALITY //////////////
%  /////////////////////////////////////////////////////////////////////////////

function [algal_cm2, dictlob_cm2] = f_adjust_total_cover (algal_cm2, dictlob_cm2, overshoot)


% algal_order = randperm(3); % process coral overgrowth of macro in random order
algal_order = randperm(2); % now limited to Lob and dict if Turf treated before


for o=1:2
    
    if (algal_order(o)+1 == 2) && (algal_cm2(algal_order(o)+1) > 0) %if look at Dictyota
        

        t12 = algal_cm2(2) - dictlob_cm2; %only consider dict without underlying lob here
        
        if (t12 > overshoot)% don't need to worry about dictlob_cm2
            algal_cm2(2) = algal_cm2(2) - overshoot; 
            break %Done by simply reducing Dictyota
            
        elseif t12 > 0; %need to use up remaining dict that's free of lob
            algal_cm2(2) = algal_cm2(2) - t12;
            overshoot = overshoot - t12;
        end
        
    elseif (algal_order(o)+1 == 3) && (algal_cm2(algal_order(o)+1) > 0) %if look at Lobophora
        
        % then look at Lobophora
        m6 = algal_cm2(3);
        
        if (m6 > overshoot) %can absorb full increase in coral by loss of lob
            if (algal_cm2(3)-overshoot) >= dictlob_cm2; %if coral excess can overgrow lob without affecting dict on top
                algal_cm2(3) = m6 - overshoot;
                break
                
            else % going to have to remove some dict too
                dictlost = overshoot - algal_cm2(3) + dictlob_cm2;%amount of dict to lose
                dictlob_cm2 = dictlob_cm2 - dictlost;
                algal_cm2(2) = algal_cm2(2) - dictlost;
                algal_cm2(3) = m6 - overshoot;
                break
            end
        else %going to lose all lob with coral overgrowth but still some expansion needed
            algal_cm2(3) = 0;
            algal_cm2(2) = algal_cm2(2)-dictlob_cm2;
            dictlob_cm2 = 0;% none left
            overshoot = overshoot - m6;
        end
        
    elseif (algal_order(o)+1 == 4) && (algal_cm2(algal_order(o)+1) > 0)
        % then look at any other kind of algae
        m7 = algal_cm2(4);
        % // if can absorb whole of overshoot with this type
        if m7 > overshoot
            algal_cm2(4) = m7 - overshoot;
            break
        else
            algal_cm2(4) = 0;
            overshoot = overshoot - m7;
        end
    end % end of if
end % end of for

