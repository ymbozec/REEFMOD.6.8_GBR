function [brnd]=betainfrnd(N, mu, sigma, nu, tau)

a=mu*(1-sigma^2)/(sigma^2);
b=a*(1-mu)/mu;
brnd=zeros(N,1);
if rand>=((1+nu)/(1+nu+tau))
    brnd(1,1)=1;
else
    for n=1:N
        if rand<=(nu/(1+nu+tau))
            brnd(n,1)=0;
        else
            brnd(n,1)=betarnd(a,b);
        end
    end
end
  


end


%If one of them is one, then the others cannot be used; have to sum up to one

% check if it is 1
% check how many are zero
% for the rest, use beta distrubtion to fill them out
%     
%     
%     find how many reefs with mro ethan one site have 1
%     detemrnine proprtion of reefs with 1
%     probability of zero, given the number of sites

%does the probability of having 1, i.e. all COTS in one site, decrease as
%we have more sites? if so, we cannot use the same probabiltiy for all
%number of sites