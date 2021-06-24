%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Identifiy Indicators of Systemic Risk (2020)
% Benny Hartwig, Christoph Meinering, Yves Schueler
% Deutsche Bundesbank
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; clc; close all
addpath(genpath('function'));



all_data_freq = {'data'}; %
all_dep1_name = {'LVRR_crises_final'};%
all_STDname = {'R2vsR2star','H2vsH2star','R2vsR2Istar','H2vsH2Istar'};

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
                    case 'R2vsR2Istar'; reg_name_stats = {'theta','stdR2','tstatR2','stdR2Istar','tstatR2Istar'};
                    case 'R2vsR2star'; reg_name_stats = {'theta','stdR2','tstatR2','stdR2star','tstatR2star'};
                    case 'H2vsH2Istar'; reg_name_stats = {'theta','stdH2','tstatH2','stdH2Istar','tstatH2Istar'};
                    case 'H2vsH2star'; reg_name_stats = {'theta','stdH2','tstatH2','stdH2star','tstatH2star'};
                end
                qreg_quantile = {'05'};%;'50'};
                
                Hname = strcat('m',num2str(m),'_p',num2str(p));
                index1st = [1 1 ];
                index2nd = [m+2 1 ];
                ML2step_name = {'LogitOLS','LogitQReg05'}; %
                vnames_row_label = {'OLS';'OLS (corrected)';'QReg'; 'QReg (corrected)'};
                idx_excel_start =  [ 1;2;3;9;8;7;6;5]*4-3 ;
                idx_excel_select =  vec([idx_excel_start, idx_excel_start+1, idx_excel_start+2,  idx_excel_start+3]');
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
                        cell2table( repmat(vnames_row_label,length(regnameOPT),1) ,'VariableNames',{strcat(depname{ii},'_LabelRow')}), ...
                        cell2table(vec([vnames_label ;NaNcell';NaNcell';NaNcell']) ,'VariableNames',{strcat(depname{ii},'_Indicator')}), ...
                        cell2table(  vec([num2signstr(ExtractStatsSum(struc.LogitOLS,regnameOPT,'theta1',rowsOPT,cols)),NaNcell,NaNcell,NaNcell]')  ,'VariableNames',{strcat(depname{ii},'_sum_theta1')}), ...
                        cell2table( vec([pval2star(ExtractStats(struc.LogitOLS,regnameOPT,'pvalLR1',index1st,rowsOPT,cols)),NaNcell,NaNcell,NaNcell]'),'VariableNames',{strcat(depname{ii},'_star1')}), ...
                        ] ;
                    
                    %% Extract OLS based on Optimal Logit
                    T = struc.LogitOLS.(regname{1}).nobs;
                    K = struc.LogitOLS.(regname{1}).K2;
                    %   array2table( round(ExtractStats(strucOLS,regnameOPT,OLSreg_name_stats{:,1},index2nd,rowsOPT,cols)',numROUND),'VariableNames',{strcat(depname{ii},'_gamma2')}), ...
                    
                    for jj = 1 : length(ML2step_name)
                        gamma_2(:,:,jj) =round( (ExtractStats(struc.(ML2step_name{jj}),regnameOPT,reg_name_stats{:,1},index2nd,rowsOPT,cols))',presc);
                        
                        std2(:,:,jj) = round( (ExtractStats(struc.(ML2step_name{jj}),regnameOPT,reg_name_stats{:,2},index2nd,rowsOPT,cols))',presc);
                        RegCI2(:,:,jj) = round( gamma_2(:,:,jj)+ std2(:,:,jj)*[ tinv(p_crit,T-K), tinv(1-p_crit,T-K)],presc);
                        tval2(:,:,jj) = ExtractStats(struc.(ML2step_name{jj}),regnameOPT,reg_name_stats{:,3},index2nd,rowsOPT,cols)';
                        
                        std2corr(:,:,jj) =round(  (ExtractStats(struc.(ML2step_name{jj}),regnameOPT,reg_name_stats{:,4},index2nd,rowsOPT,cols))',presc);
                        RegCI2corr(:,:,jj) =round(  gamma_2(:,:,jj)+ std2corr(:,:,jj)*[ tinv(p_crit,T-K), tinv(1-p_crit,T-K)],presc);
                        tval2corr(:,:,jj) = ExtractStats(struc.(ML2step_name{jj}),regnameOPT,reg_name_stats{:,5},index2nd,rowsOPT,cols)';
                    end
                    try
                    out.Stage2 = [...
                        array2table( vec([ gamma_2(:,:,1) NaNmat , gamma_2(:,:,2) NaNmat]') ,'VariableNames',{strcat(depname{ii},'_gamma2_coeff')}), ...
                        array2table( vec([ std2(:,:,1) std2corr(:,:,1), std2(:,:,2) std2corr(:,:,2)]') ,'VariableNames',{strcat(depname{ii},'_gamma2_std')}), ...
                        cell2table( vec([ pval2star(tcdf(tval2(:,:,1),T-K)) ,pval2star(tcdf(tval2corr(:,:,1),T-K)), pval2star(tcdf(tval2(:,:,2),T-K)) ,pval2star(tcdf(tval2corr(:,:,2),T-K))]')   ,'VariableNames',{strcat(depname{ii},'_gamma2_pval')}) , ...
                        array2table( [vec([RegCI2(:,1,1)'; RegCI2corr(:,1,1)'; RegCI2(:,1,2)'; RegCI2corr(:,1,2)']) vec([RegCI2(:,2,1)'; RegCI2corr(:,2,1)';RegCI2(:,2,2)'; RegCI2corr(:,2,2)'])] ,'VariableNames',{strcat(depname{ii},'_gamma2_CIlow'),strcat(depname{ii},'_gamma2_CIup')}), ...
                        cell2table(  repmat(NaNcell,4,1),'VariableNames',{strcat(depname{ii},'_Empty')}), ...
                        ] ;
                    catch
                           out.Stage2 = [...
                        array2table( vec([ gamma_2(:,:,1) NaNmat , gamma_2(:,:,2) NaNmat]') ,'VariableNames',{strcat(depname{ii},'_gamma2_coeff')}), ...
                        array2table( vec([ std2(:,:,1) std2corr(:,:,1), std2(:,:,2) std2corr(:,:,2)]') ,'VariableNames',{strcat(depname{ii},'_gamma2_std')}), ...
                        cell2table( vec([ pval2star(tcdf(zeros(size(tval2(:,:,1))),T-K)) ,pval2star(tcdf(zeros(size(tval2corr(:,:,1))),T-K)), pval2star(tcdf(zeros(size(tval2(:,:,2))),T-K)) ,pval2star(tcdf(zeros(size(tval2corr(:,:,2))),T-K))]')   ,'VariableNames',{strcat(depname{ii},'_gamma2_pval')}) , ...
                        array2table( [vec([RegCI2(:,1,1)'; RegCI2corr(:,1,1)'; RegCI2(:,1,2)'; RegCI2corr(:,1,2)']) vec([RegCI2(:,2,1)'; RegCI2corr(:,2,1)';RegCI2(:,2,2)'; RegCI2corr(:,2,2)'])] ,'VariableNames',{strcat(depname{ii},'_gamma2_CIlow'),strcat(depname{ii},'_gamma2_CIup')}), ...
                        cell2table(  repmat(NaNcell,4,1),'VariableNames',{strcat(depname{ii},'_Empty')}), ...
                        ] ;
                    end
                    
                    
                    output.(Hname).(STDname).(depname{ii}).ML2step = [out.Stage1(idx_excel_select,:) out.Stage2(idx_excel_select,:)  ];
                    
                    
                    
                end
                
                
                
                
                
                Folder = cd;
                
                
                switch data_freq
                    case {'data_semiannual_avg', 'data_semiannual_last'}
                        Hname = 'm2_p0';
                        
                        out2xls = [output.(Hname).(STDname).F0.ML2step , ...
                            output.(Hname).(STDname).F1.ML2step(:,2:end),...
                            output.(Hname).(STDname).F2.ML2step(:,2:end),...
                            output.(Hname).(STDname).F3.ML2step(:,2:end),...
                            output.(Hname).(STDname).F4.ML2step(:,2:end),...
                            output.(Hname).(STDname).F6.ML2step(:,2:end)] ;
                    case 'data' %_quarterly
                        Hname = 'm4_p0';
                        out2xls = [output.(Hname).(STDname).F0.ML2step, ...
                            output.(Hname).(STDname).F1.ML2step(:,2:end),...
                            output.(Hname).(STDname).F2.ML2step(:,2:end),...
                            output.(Hname).(STDname).F4.ML2step(:,2:end),...
                            output.(Hname).(STDname).F6.ML2step(:,2:end),...
                            output.(Hname).(STDname).F8.ML2step(:,2:end),...
                            output.(Hname).(STDname).F12.ML2step(:,2:end)] ;
                        
                end
                fullFileName = fullfile(strcat(Folder,'\table\robustness\',STDname), strcat(Hname,'_',data_freq(6:end),'_',dep1_name ,'.xlsx'));
                mkdir(strcat(Folder,'\table\robustness\',STDname))
                writetable(out2xls,fullFileName,'Sheet',country,'WriteVariableNames',true)
            end
        end
        
    end
end