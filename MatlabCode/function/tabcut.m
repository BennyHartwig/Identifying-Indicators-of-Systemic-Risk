%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Identifiy Indicators of Systemic Risk (2020)
% Benny Hartwig, Christoph Meinering, Yves Schueler
% Deutsche Bundesbank
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [tab ] = tabcut(tab,idxStart,idxEnd)

vnames = fields(tab);
if ~isempty(idxStart) || ~isempty(idxEnd)
    for ii = 1 : length(vnames)
        tab.(vnames{ii})(idxStart:idxEnd,:) = [];
    end
end