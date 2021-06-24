%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Identifiy Indicators of Systemic Risk (2020)
% Benny Hartwig, Christoph Meinering, Yves Schueler
% Deutsche Bundesbank
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [bracnum ] = num2bracnum(num)

for ii = 1 : length(num)
    bracnum(ii,1) = {strcat('[',num2str(num(ii)),']')};
end