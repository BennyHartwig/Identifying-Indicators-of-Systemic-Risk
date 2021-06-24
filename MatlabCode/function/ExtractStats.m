function output = ExtractStats(struc,regnames,fieldname,index,rows,cols);
kk = index(1); ll = index(2);
for ii = 1 : length(regnames)
    aux(ii,1) = struc.(regnames{ii}).(fieldname)(kk,ll);
end

output = reshape(aux,rows,cols);
