%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Identifiy Indicators of Systemic Risk (2020)
% Benny Hartwig, Christoph Meinering, Yves Schueler
% Deutsche Bundesbank
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ordinal  ] = pval2ordinal(pval)

for ii = 1 : length(pval)
    if pval(ii) < .01
        ordinal(ii) = 3;
    elseif pval(ii) <.05
        ordinal(ii) = 2;
        
    elseif pval(ii) <.10
        ordinal(ii) = 1;
    else
        ordinal(ii) = 0;
    end
end