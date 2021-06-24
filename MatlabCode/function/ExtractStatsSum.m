function output = ExtractStatsSum(struc,regnames,fieldname,rows,cols);

for ii = 1 : length(regnames)
    aux(ii,1) = sum(struc.(regnames{ii}).(fieldname)(1:end-1));
end

output = reshape(aux,rows,cols);
