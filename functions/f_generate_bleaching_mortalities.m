% -------------------------------------------------------------------------
% Y.-M. Bozec, MSEL, created 09/2013.
%
% GENERATES BLEACHING MORTALITIES FROM DHW
% Now based on Hughes et al 2018 relationship between DHW and initial
% mortality
% -------------------------------------------------------------------------

function mortalities = f_generate_bleaching_mortalities(BleachingLinearModel,predicted_DHWs,deterministic_bleaching)

%  New function to generate whole-colony mortality values
%  from Hughes et al. (2018) relationship between initial mortality and
%  DHW. The relationship is a linear model contained in 'BleachingLinearModel'
%  which is stored in 'BleachingModelHughes.mat'

%  ADAPTATION: predicted_DHWs will be a matrix (for each Topt), not just a row
[l,c] = size(predicted_DHWs);
mortalities = zeros(l,c);

if deterministic_bleaching == 1
    % always same mortality from a given DHW from Hughes model 
    for j=1:c        
        mortalities(:,j) = exp(predict(BleachingLinearModel,predicted_DHWs(:,j)))-1;
    end
    
else
    % generate random mortalities from Hughes model (Gaussian random
    % noise with SD of model error)
    for j=1:c        
        mortalities(:,j) = exp(random(BleachingLinearModel,predicted_DHWs(:,j)))-1;
    end
    
end

mortalities = mortalities/100; % convert as proportion (note this can exceed 1)
mortalities(mortalities<0) = 0; % (can be negative when randomly generated from DHW=3)
% mortalities(mortalities>1) = 1; % now capped in f_bleaching_new


%% TESTING SECTION
% test = zeros(1,1000)
% DHW = 10*rand(1,1000)
% 
% for i=1:1000 
%     
%     mortalities(i) = exp(random(BleachingLinearModel,DHW(i)))-1;
%     
% end
% 
% mortalities = mortalities/100;
% mortalities(mortalities<0) = 0;
% mortalities(mortalities>1) = 1;
% mortalities = 1-(1- mortalities).^mortality_offset ;
% figure; plot(DHW,mortalities,'o')
% line( [4 4], [0 1])
% xlabel('DHW')
% ylabel('Initial mortality')
% 
% pred = zeros(1,10);
% for k=1:10
%     
%     pred(1,k) = exp(predict(BleachingLinearModel,k))-1;
% end
% 
% figure;plot(1:10,pred,'.')