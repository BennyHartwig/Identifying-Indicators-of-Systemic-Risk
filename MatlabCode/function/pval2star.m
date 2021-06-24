%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Identifiy Indicators of Systemic Risk (2020)
% Benny Hartwig, Christoph Meinering, Yves Schueler
% Deutsche Bundesbank
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [star  ] = pval2star(pval)

for ii = 1 : length(pval)
    if pval(ii) < .01
        star(ii,1) = {'***'};
    elseif pval(ii) <.05
        star(ii,1) = {'**'};
        
    elseif pval(ii) <.10
        star(ii,1) = {'*'};
    else
        star(ii,1) = {''};
    end
end