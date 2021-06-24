%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Identifiy Indicators of Systemic Risk (2020)
% Benny Hartwig, Christoph Meinering, Yves Schueler
% Deutsche Bundesbank
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; clc; close all
addpath(genpath('function'));

Folder = cd;


all_data_freq = {'data'};%{'data_semiannual_avg','data_quarterly'};
all_dep1_name ={'LVRR_crises_final'};%'RomerRomer','LaevenValencia','ReinhartRogoff','ESRB'};
all_STDname = {'R2','H2','R2Istar','H2Istar','R2star','H2star'};

for aa = 1 : length(all_data_freq)
    for bb = 1 : length(all_dep1_name)
        data_freq = all_data_freq{aa}; % #data_quarterly  #data_semiannual_avg #data_semiannual_last
        dep1_name = all_dep1_name{bb}; % #LaevenValencia #RomerRomer #ReinhartRogoff #ESRB
        logit_type = 'Binary'; % #Binary #Ordered
        
        
        
        
               vnames = {'cqgap','cqgap_ham','Mend_cyc','Mend_cyc_1600','Mend_cyc_rt',...
                'reer','cabgdp','spreads',...
                'FCycle','FCyc_dreh','FCyc_dreh_rt',...
                'd1lrtcredit','d4lrtcredit','d8lrtcredit','d12lrtcredit',...
                'd1lRPPI','d4lRPPI','d8lRPPI','d12lRPPI',...
                'd1lrstock','d4lrstock','d8lrstock','d12lrstock',...
                'd1lcabgdp','d4lcabgdp','d8lcabgdp','d12lcabgdp',...
                'd1lreer','d4lreer','d8lreer','d12lreer',...
                'd1lrbondpr','d4lrbondpr'};
               vnames_label = {'cqgap','cqgap_ham','Mend_cyc','Mend_cyc_1600','Mend_cyc_rt',...
                'reer','cabgdp','spreads',...
                'FCycle','FCyc_dreh','FCyc_dreh_rt',...
                'd1lrtcredit','d4lrtcredit','d8lrtcredit','d12lrtcredit',...
                'd1lRPPI','d4lRPPI','d8lRPPI','d12lRPPI',...
                'd1lrstock','d4lrstock','d8lrstock','d12lrstock',...
                'd1lcabgdp','d4lcabgdp','d8lcabgdp','d12lcabgdp',...
                'd1lreer','d4lreer','d8lreer','d12lreer',...
                'd1lrbondpr','d4lrbondpr'};
                pOPT =0; mOPT = 4;
                horizon_legend = num2cell([0:20]);      
        
        
        
        
 cnames = {'AUS','AUT','BEL','CAN','DNK','FIN','FRA','DEU','GRC','IRL','ITA','JPN','NLD','NZL','NOR','PRT','ESP','SWE','CHE','GBR','USA',...
    'BRA','BGR','CHL','CHN','COL','CZE','HKG','HUN','IND','IDN','ISR','KOR','LVA',...
    'MYS','MEX','PHL','POL','RUS','SGP','SVK','SVN','ZAF','THA','TUR'}; %
  
        
        for cc = 1 : length(cnames)
            country = cnames{cc};
            
            load(strcat('results/smpl1970/results_',country,'_LogitQReg05_',dep1_name,'_d1lRGDP'))
            load(strcat('results/smpl1970/results_',country,'_LogitOLS_',dep1_name,'_d1lRGDP'))
            results.LogitOLS = results_LogitOLS;
            results.LogitQReg05 = results_LogitQReg05;
            
            
            %% Inference 1. Stage
            for zz = 1 : length(all_STDname)
                STDname = all_STDname{zz};
                
                depname = fieldnames(results.LogitOLS.m0_p0);
                regname = fieldnames(indexStruc);
                m = mOPT;
                p = pOPT;
                p_crit = 0.1;
                name_optIC = 'BIC';
                % options for first stage
                rnames1st = {'sum_theta','star',name_optIC};
                
                cols =size(vnames,2); % number of candidate variables
                rows = length(regname)/cols;
                rowsOPT = 1;
                
                col_lag = repmat([[0 :rows-1 ]]',length(rnames1st),1);
                
                rnames1st = vec(repmat(rnames1st,rowsOPT,1));
                numROUND = 2;
                
                % options for second stage
                rnames2nd = {name_optIC,'coeff','se','tstat','pval',strcat('pval<',num2str(p_crit))};
                col_lag = repmat([[0 :rows-1 ]]',length(rnames2nd),1);
                
                rnames2nd = vec(repmat(rnames2nd,rowsOPT,1));
                
                switch STDname
                    case 'R2'; reg_name_stats = {'theta','stdR2','tstatR2'};
                    case 'R2Istar'; reg_name_stats = {'theta','stdR2Istar','tstatR2Istar'};
                    case 'R2star'; reg_name_stats = {'theta','stdR2star','tstatR2star'};
                    case 'H2'; reg_name_stats = {'theta','stdH2','tstatH2'};
                    case 'H2Istar'; reg_name_stats = {'theta','stdH2Istar','tstatH2Istar'};
                    case 'H2star'; reg_name_stats = {'theta','stdH2star','tstatH2star'};
                end
                qreg_quantile = {'05'};
                
                Hname = strcat('m',num2str(m),'_p',num2str(p));
                index1st = [1 1 ];
                index2nd = [m+2 1 ];
                ML2step_name = {'LogitOLS','LogitQReg05'};
                vnames_row_label = {'OLS';'OLS (corrected)';'QReg'; 'QReg (corrected)'};
                idx_excel_select =  [1;2;3;14;12;11;10;9];% [ 6;10;14;18;3 ;2;1;4] ;
                presc = 2;
                for ii = 1 : length(depname)
                    clear out
                    
                    %% Extract Optimal Logit
                    for jj = 1 : length(ML2step_name)
                        struc.(ML2step_name{jj}) = results.(ML2step_name{jj}).(Hname).(depname{ii});
                    end
                    
                    T = struc.LogitOLS.(regname{1}).nobs;
                    K = struc.LogitOLS.(regname{1}).K2;
                    
                    
                    rows = length(regname)/cols;
                    optIC = ExtractStats(struc.LogitOLS,regname,strcat(name_optIC,'1'),index1st,rows,cols);
                    [val lagg]  = min(optIC );
                    for  kk = 1 : length(val)
                        if isempty(find(optIC(:,kk) == val(kk)))
                            idx(kk) = lagg(kk)+(kk-1)*rows;
                        else
                            idx(kk) = find(optIC(:,kk) == val(kk))+(kk-1)*rows;
                        end
                    end
                    regnameOPT =regname(idx);
                    
                    NaNmat = NaN(length(regnameOPT),1);
                    NaNcell = repmat({' '} ,length(regnameOPT),1);
                    
                    out.Stage1= [ ...
                        cell2table(vnames_label'  ,'VariableNames',{strcat(depname{ii},'_Indicator')}), ...
                        array2table( round(ExtractStatsSum(struc.LogitOLS,regnameOPT,'theta1',rowsOPT,cols)',presc)  ,'VariableNames',{strcat(depname{ii},'_sum_theta1')}), ...
                        cell2table( pval2star(ExtractStats(struc.LogitOLS,regnameOPT,'pvalLR1',index1st,rowsOPT,cols)),'VariableNames',{strcat(depname{ii},'_star1')}), ...
                        cell2table( num2bracnum([lagg-1 ]'),'VariableNames',{strcat(depname{ii},'_',name_optIC)}),...
                        ] ;
                    
                    %% Extract OLS based on Optimal Logit
                    T = struc.LogitOLS.(regname{1}).nobs;
                    K = struc.LogitOLS.(regname{1}).K2;
                    
                    for jj = 1 : length(ML2step_name)
                        gamma_2 =round( (ExtractStats(struc.(ML2step_name{jj}),regnameOPT,reg_name_stats{:,1},index2nd,rowsOPT,cols))',presc);
                        std2 = round( (ExtractStats(struc.(ML2step_name{jj}),regnameOPT,reg_name_stats{:,2},index2nd,rowsOPT,cols))',presc);
                        tval2 = ExtractStats(struc.(ML2step_name{jj}),regnameOPT,reg_name_stats{:,3},index2nd,rowsOPT,cols)';
                        
                      try  
                        out.(ML2step_name{jj}) = [...
                            array2table( gamma_2 ,'VariableNames',{strcat(depname{ii},'_',ML2step_name{jj},'_gamma2_coeff')}), ...
                            cell2table(  pval2star(tcdf(tval2(:,:,1),T-K))  ,'VariableNames',{strcat(depname{ii},'_',ML2step_name{jj},'_gamma2_pval')}) , ...
                            cell2table(  NaNcell,'VariableNames',{strcat(depname{ii},'_',ML2step_name{jj},'_Empty')}), ...
                            ] ;
                      catch
                        out.(ML2step_name{jj}) = [...
                            array2table( gamma_2 ,'VariableNames',{strcat(depname{ii},'_',ML2step_name{jj},'_gamma2_coeff')}), ...
                            cell2table(  pval2star(tcdf(zeros(size(tval2(:,:,1))),T-K))  ,'VariableNames',{strcat(depname{ii},'_',ML2step_name{jj},'_gamma2_pval')}) , ...
                            cell2table(  NaNcell,'VariableNames',{strcat(depname{ii},'_',ML2step_name{jj},'_Empty')}), ...
                            ] ;
                      end
                    end
                    output.(country).(Hname).(STDname).(depname{ii}).ML2step = [out.Stage1(idx_excel_select,:) out.LogitOLS(idx_excel_select,:)  out.LogitQReg05(idx_excel_select,:)  ];
                    
                    
                    
                end
            end
        end
         cnames = {'AUS','AUT','BEL','CAN','DNK','FIN','FRA','DEU','GRC','IRL','ITA','JPN','NLD','NZL','NOR','PRT','ESP','SWE','CHE','GBR','USA',...
    'BRA','BGR','CHL','CHN','COL','CZE','HKG','HUN','IND','IDN','ISR','KOR','LVA',...
    'MYS','MEX','PHL','POL','RUS','SGP','SVK','SVN','ZAF','THA','TUR'};

        tmp = transposeTable( output.(country).(Hname).(STDname).(depname{ii}).ML2step);
        vnames_indicator = tmp.Properties.VariableNames(2:end);

        
        for zz = 1 : length(all_STDname)
            STDname = all_STDname{zz};
            
            
            
            for jj = 1 : length(vnames_indicator)
                cellBigRows = [];
                for ii = 1 : length(depname)
                    cellBigColumns = [];
                    for cc = 1 : length(cnames)
                        country = cnames{cc};
                        if isempty( find(strcmp(fields(output),cnames{cc})==1))
                            celltmp = repmat({' '},3,3);
                        else
                            tmp = transposeTable( output.(country).(Hname).(STDname).(depname{ii}).ML2step);
                            celltmp = (tmp.(vnames_indicator{jj}))';
                        end
                        cellBigColumns = [cellBigColumns; celltmp];
                    end
                    cellBigRows = [cellBigRows cellBigColumns];
                end
                
                vnamesBigColumns = {};
                
                for cc = 1 : length(cnames)
                    country = cnames{cc};
                    celltmp = country;%strcat(country,{'_coeff','_pval','_BIC'});
                    vnamesBigColumns = [vnamesBigColumns; celltmp];
                    
                end
                
                vnames_label_horizon = vec([strcat(depname,'_Logit1') strcat(depname,'_Logit2') strcat(depname,'_Logit3') strcat(depname,'_OLS') strcat(depname,'_OLS') strcat(depname,'_OLS') strcat(depname,'_QReg') strcat(depname,'_QReg') strcat(depname,'_QReg')]');
                out2xls = cell2table([vnames_label_horizon'; cellBigRows]','VariableNames',[{'horizon_type'};vnamesBigColumns]);
                
                
                fullFileName = fullfile(strcat(Folder,'\table\mainresult\',STDname), strcat(Hname,'_',data_freq(6:end),'_',dep1_name ,'.xlsx'));
                mkdir(strcat(Folder,'\table\mainresult\',STDname))
                writetable(out2xls,fullFileName,'Sheet',vnames_indicator{jj},'WriteVariableNames',true)
            end
        end
        clear output
        
    end
end