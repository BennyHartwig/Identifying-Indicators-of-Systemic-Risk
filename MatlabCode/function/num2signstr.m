%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Identifiy Indicators of Systemic Risk (2020)
% Benny Hartwig, Christoph Meinering, Yves Schueler
% Deutsche Bundesbank
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [signstr ] = num2signstr(num)

signnum = sign(num);

for ii = 1 : length(signnum)
    if signnum(ii) == 1
        signstr(ii,1) = {'+'};
    elseif signnum(ii) == -1
        signstr(ii,1) = {'-'};
    elseif signnum(ii) == 0
        signstr(ii,1) = {'0'};
    else 
        signstr(ii,1) = {''};

    end
end
