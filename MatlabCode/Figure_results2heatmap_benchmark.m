%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Identifiy Indicators of Systemic Risk (2020)
% Benny Hartwig, Christoph Meinering, Yves Schueler
% Deutsche Bundesbank
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear all; clc; close all
%% 
addpath(genpath('function'));
set(groot,'defaulttextinterpreter','latex');  
set(groot, 'defaultAxesTickLabelInterpreter','latex');  
set(groot, 'defaultLegendInterpreter','latex'); 

Folder = cd;




RegType =2; RegTypeName1 = 'QReg05'; RegTypeName2= 'OLS';


p_crit = 0.1;


data_freq = 'data'; %
dep1_name = 'LVRR_crises_final'; % 
dep2_name = 'd1lRGDP';

all_STDname = {'R2star','R2Istar','R2'}; % 
% Benchmark: R2star

all_TestType = {'Hierarchical','Joint','JointCounterfactual'};
% Benchmark: Hierarchical

           
logit_type = 'Binary'; % 
            
            
            
load(data_freq)

                 vnames_dict = {'cqgap','Basel III credit-to-GDP gap';...
                     'cqgap_ham','Hamilton-filtered credit-to-GDP ratio';...
                     'Mend_cyc','Credit-boom indicator (Mendoza and Terrones (2008,2014))';...
                     'Mend_cyc_1600','Credit-boom indicator with 1600 (Mendoza and Terrones (2008,2014))';...
                     'Mend_cyc_rt','Credit-boom indicator (Mendoza and Terrones (2008,2014)), real time';...
                     'reer','Real effective exchange rate';...
                     'cabgdp','Current account to GDP ratio';...
                     'spreads','Country risk spreads';...
                     'FCycle','Composite financial cycle (Schüler et al. (2020))';...
                     'FCyc_dreh','Composite financial cycle (Drehmann et al. (2012))';...
                     'FCyc_dreh_rt','Composite financial cycle (Drehmann et al. (2012)), real time';...
                     'd1lrtcredit','Real credit growth (qoq)';...
                     'd4lrtcredit','Real credit growth (yoy)';...
                     'd8lrtcredit','Real credit growth (2yo2y)';...
                     'd12lrtcredit','Real credit growth (3yo3y)';...
                     'd1lRPPI','Real house price growth (qoq)';...
                     'd4lRPPI','Real house price growth (yoy)';...
                     'd8lRPPI','Real house price growth (2yo2y)';...
                     'd12lRPPI','Real house price growth (3yo3y)';...
                     'd1lrstock','Real stock price growth (qoq)';...
                     'd4lrstock','Real stock price growth (yoy)';...
                     'd8lrstock','Real stock price growth (2yo2y)';...
                     'd12lrstock','Real stock price growth (3yo3y)';...
                     'd1lcabgdp','Current account balance (% of GDP) (qoq)';...
                     'd4lcabgdp','Current account balance (% of GDP) (yoy)';...
                     'd8lcabgdp','Current account balance (% of GDP) (2yo2y)';...
                     'd12lcabgdp','Current account balance (% of GDP) (3yo3y)';...
                     'd1lreer','REER (qoq)';...
                     'd4lreer','REER (yoy)';...
                     'd8lreer','REER (2yo2y)';...
                     'd12lreer','REER (3yo3y)';...
                     'd1lrbondpr','Real bond price growth (qoq)';...
                     'd4lrbondpr','Real bond price growth (yoy)'};
                   
vnames = vnames_dict(:,1); vnames_label = vnames_dict(:,2);

                    
                    pOPT =0; mOPT = 4;
                    horizon_legend = num2cell([0:20]); %0:20
                    xlabel_name = {'Quarters ahead'};
                    stage2_sign = -1; % test left tail

            cnames = {'AUS','AUT','BEL','CAN','DNK','FIN','FRA','DEU','GRC','IRL','ITA','JPN','NLD','NZL','NOR','PRT','ESP','SWE','CHE','GBR','USA',...
            'BRA','BGR','CHL','CHN','COL','CZE','HKG','HUN','IND','IDN','ISR','KOR','LVA',...
            'MYS','MEX','PHL','POL','RUS','SGP','SVK','SVN','ZAF','THA','TUR'};
                    
            
            
            
     

            
            
            for cc = 1 :  length(cnames)
                country = cnames{cc};
%                 load(strcat('../results/results_',data_freq(6:end),'_',country,'_',dep1_name ,'_',dep2_name))
                
                % new import
                ML2step_name = [];
                mname = strcat('Logit',RegTypeName1); ML2step_name =[ ML2step_name , {mname}]; 
                load(strcat('results/smpl1970/results_',country,'_',mname,'_',dep1_name ,'_',dep2_name))
                eval(['results.(mname)  = results_',mname,' ;'] );
                clear(strcat('results_',mname))
                mname = strcat('Logit',RegTypeName2);  ML2step_name =[ ML2step_name , {mname}];
                load(strcat('results/smpl1970/results_',country,'_',mname,'_',dep1_name ,'_',dep2_name))
                eval(['results.(mname)  = results_',mname,' ;'] );
                clear(strcat('results_',mname))

                
                
                for yy = 1 : length(all_TestType)
                    for zz = 1 : length(all_STDname)
                        TestType = all_TestType{yy};
                        STDname = all_STDname{zz};
                        
                        
                        
                        
                        %% Inference 1. Stage
                        
                        depname = fieldnames(results.LogitOLS.m0_p0);
                        regname = fieldnames(indexStruc);
                        name_optIC = 'BIC';
                        % options for first stage
                        
                        cols =size(vnames,1); % number of candidate variables
                        rows = length(regname)/cols;
                        rowsOPT = 1;
                        
                        
                        numROUND = 2;
                        
                        % options for second stage
                        
                        
                        switch STDname
                            case 'R2Istar'; reg_name_stats = {'theta','stdR2Istar','tstatR2Istar'};
                            case 'R2star'; reg_name_stats = {'theta','stdR2star','tstatR2star'};
                            case 'R2'; reg_name_stats = {'theta','stdR2','tstatR2'};
                            case 'H2Istar'; reg_name_stats = {'theta','stdH2Istar','tstatH2Istar'};
                            case 'H2star'; reg_name_stats = {'theta','stdH2star','tstatH2star'};
                            case 'H2'; reg_name_stats = {'theta','stdH2','tstatH2'};
                        end
                        
                        
                        
                        m = mOPT;
                        p = pOPT;
                        Hname = strcat('m',num2str(m),'_p',num2str(p));
                        index1st = [1 1 ];
                        index2nd = [m+2 1 ];
                        % old import
%                         ML2step_name = {'LogitOLS','LogitQReg05','LogitQReg95'};
                        for ii = 1 : length(depname)
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
                            
                            for kif = 1:length(vnames)
                            save_gfc.(cnames{cc}).(depname{ii}).(vnames{kif}) = results.LogitOLS.m4_p0.(depname{ii}).(regnameOPT{kif});
                            end
                            
                            LogitSign=  sign(ExtractStatsSum(struc.LogitOLS,regnameOPT,'theta1',rowsOPT,cols))'   ;
                            LogitPval = double(pval2ordinal(ExtractStats(struc.LogitOLS,regnameOPT,'pvalLR1',index1st,rowsOPT,cols))'>0) ;
                            LogitPval(isnan(LogitSign)) = NaN;

                            %% Extract OLS based on Optimal Logit
                            switch stage2_sign
                                case -1
                                    switch RegType
                                        case 1
                                            TypeName1 = strcat('Logit',RegTypeName1);
                                            output_name = strcat('ML2step',RegTypeName1);
                                            
                                            RegPval =  (tcdf(ExtractStats(struc.(TypeName1),regnameOPT,reg_name_stats{:,3},index2nd,rowsOPT,cols),T-K))'<p_crit ;
                                            
                                            
                                        case 2
                                            TypeName1  = strcat('Logit',RegTypeName1);
                                            TypeName2  = strcat('Logit',RegTypeName2);
                                            output_name = strcat('ML2step',RegTypeName1,'vs',RegTypeName2);
                                            
                                             try
                                            RegPval_1 = double( (tcdf(ExtractStats(struc.(TypeName1),regnameOPT,reg_name_stats{:,3},index2nd,rowsOPT,cols),T-K))'<p_crit );
                                             catch
                                            disp(country)
                                            RegPval_1 =  zeros(cols,1);    
                                            end
                                            RegTheta_1 =   (ExtractStats(struc.(TypeName1),regnameOPT,reg_name_stats{:,1},index2nd,rowsOPT,cols)');
                                            RegStd_1 =   (ExtractStats(struc.(TypeName1),regnameOPT,reg_name_stats{:,2},index2nd,rowsOPT,cols)');
                                            RegCI_1 = RegTheta_1+ [ RegStd_1*tinv(p_crit,T-K), RegStd_1*tinv(1-p_crit,T-K)];
                                            
                                            try
                                            RegPval_2 =  double((tcdf((ExtractStats(struc.(TypeName2),regnameOPT,reg_name_stats{:,3},index2nd,rowsOPT,cols)),T-K))'<p_crit) ;
                                            catch
                                             disp(country)
                                            RegPval_2 =  zeros(cols,1);    
                                            end
                                            RegTheta_2 =   (ExtractStats(struc.(TypeName2),regnameOPT,reg_name_stats{:,1},index2nd,rowsOPT,cols)');
                                            RegStd_2 =   (ExtractStats(struc.(TypeName2),regnameOPT,reg_name_stats{:,2},index2nd,rowsOPT,cols)');
                                            RegCI_2 = RegTheta_2+ [ RegStd_2*tinv(p_crit,T-K), RegStd_2*tinv(1-p_crit,T-K)];
                                            
                                            % insert nan
                                            RegPval_1(isnan(RegTheta_1)) = NaN;
                                            RegPval_2(isnan(RegTheta_2)) = NaN;

                                            %
                                            
                                            idx_pval = find(RegPval_1 == 1); idx_ratio = find( RegCI_1(:,2)<RegCI_2(:,1));% find( RegTheta_1./RegTheta_2>1);
                                            idx_Q = intersect( idx_ratio, idx_pval );
                                            
                                            RegPval =  double((RegPval_1+RegPval_2)>0) ; RegPval(idx_Q) = 2;
                                            
                                            
                                            
                                            
                                    end
                                case 1 % Test positive sign
                                    switch RegType
                                        case 1 
                                            TypeName1 = strcat('Logit',RegTypeName1);
                                            output_name = strcat('ML2step',RegTypeName1);
                                            
                                            RegPval =  (1-(tcdf(ExtractStats(struc.(TypeName1),regnameOPT,reg_name_stats{:,3},index2nd,rowsOPT,cols),T-K))')<p_crit ;
                                            
                                            
                                        case 2
                                            TypeName1  = strcat('Logit',RegTypeName1);
                                            TypeName2  = strcat('Logit',RegTypeName2);
                                            output_name = strcat('ML2step',RegTypeName1,'vs',RegTypeName2);
                                            
                                            RegPval_1 =  double((1-(tcdf(ExtractStats(struc.(TypeName1),regnameOPT,reg_name_stats{:,3},index2nd,rowsOPT,cols),T-K))')<p_crit) ;
                                            RegTheta_1 =   (ExtractStats(struc.(TypeName1),regnameOPT,reg_name_stats{:,1},index2nd,rowsOPT,cols)');
                                            RegStd_1 =   (ExtractStats(struc.(TypeName1),regnameOPT,reg_name_stats{:,2},index2nd,rowsOPT,cols)');
                                            RegCI_1 = RegTheta_1+ [ RegStd_1*tinv(p_crit,T-K), RegStd_1*tinv(1-p_crit,T-K)];
                                            
                                            
                                            RegPval_2 = double((1- (tcdf(ExtractStats(struc.(TypeName2),regnameOPT,reg_name_stats{:,3},index2nd,rowsOPT,cols),T-K))')<p_crit );
                                            RegTheta_2 =   (ExtractStats(struc.(TypeName2),regnameOPT,reg_name_stats{:,1},index2nd,rowsOPT,cols)');
                                            RegStd_2 =   (ExtractStats(struc.(TypeName2),regnameOPT,reg_name_stats{:,2},index2nd,rowsOPT,cols)');
                                            RegCI_2 = RegTheta_2+ [ RegStd_2*tinv(p_crit,T-K), RegStd_2*tinv(1-p_crit,T-K)];
                                            
                                            % insert nan
                                            RegPval_1(isnan(RegTheta_1)) = NaN;
                                            RegPval_2(isnan(RegTheta_2)) = NaN;
                                            
                                            
                                            idx_pval = find(RegPval_1 == 1); idx_ratio = find( RegCI_2(:,2)<RegCI_1(:,1));% OLS_UPPER < QREG_LOWER 
                                            idx_Q = intersect( idx_ratio, idx_pval );
                                            
                                            RegPval =  double((RegPval_1+RegPval_2)>0) ; RegPval(idx_Q) = 2;
                                            
                                            
                                            
                                            
                                    end
                            end
                            
                            
                            
                            
                            clear idx_logit idx_logitReg  LogitSignReg  LogitReg
                            idx_logit = find(LogitPval >0);
                            idx_logitReg = find(RegPval>0);
                            idx_logitRegJoint = find(LogitPval.*RegPval>0);
                            LogitSignReg = LogitSign.*RegPval;
                            
                            LogitReg = zeros(length(LogitPval),1);
                            LogitReg(isnan(LogitPval)) = NaN;
                            switch TestType
                                case 'Hierarchical'
                                    LogitReg(idx_logit) = -0.5;
                                    LogitReg(idx_logitRegJoint) = LogitSignReg(idx_logitRegJoint);
                                case 'Joint'
                                    LogitReg(idx_logit) = -0.5;
                                    LogitReg(idx_logitReg) = LogitSignReg(idx_logitReg);
                                    
                                case 'JointCounterfactual'
                                    LogitReg(idx_logitReg) = 0.5;
                                    LogitReg(idx_logit) = -0.5;
                                    LogitReg(idx_logitRegJoint) = LogitSignReg(idx_logitRegJoint);
                                    
                            end
                            
                            for kif = 1:length(vnames)
                            save_gfc_color.(cnames{cc}).(depname{ii}).(vnames{kif}) = LogitReg(kif);
                            end
                            
                            
                            output.(TestType).(STDname).(country).(Hname).(output_name)(:,ii) = [LogitReg ];
                            
                            
                            
                            
                        end
                        
                        
                        
                    end
                end
            end
            

cnames_all = cnames;

            switch RegType
                case 1
                    map = [0.7 0.7 1; ...
                        .8 .8 .8 ; 1 1 1;    0.4 1 0.4; ...
                        1 0.7 0.7];
                    ColorLimit = [-1 ; 1];
                    
                case 2
                    map = [  0 0 1  ; 0.7 0.7 1; ...
                        .8 .8 .8 ; 1 1 1;    0.4 1 0.4;...
                        1 0.7 0.7 ;1 0 0];
                    ColorLimit = [-2 ; 2];
                    
            end
            width = 8/1.5 ; height = 8;
%HEATMAP
            for yy = 1 : length(all_TestType)
                for zz = 1 : length(all_STDname)
                    TestType = all_TestType{yy};
                    STDname = all_STDname{zz};
                    
                    
                    
                    for ii = 1 : length(vnames)
                        for cc = 1 : length(cnames_all)
                            if isempty( find(strcmp(fields(output.(TestType).(STDname)),cnames_all{cc})==1))
                                output_fig.(TestType).(STDname).(Hname).(output_name).(vnames{ii})(cc,:) = NaN(1,length(horizon_legend));
                            else
                                output_fig.(TestType).(STDname).(Hname).(output_name).(vnames{ii})(cc,:) = output.(TestType).(STDname).(cnames_all{cc}).(Hname).(output_name)(ii,:);
                            end
                        end
                        h = figure('Renderer', 'painters','units','inch', 'Position', [1 1 width height],'color','w')
                        heatmap(horizon_legend,cnames_all,output_fig.(TestType).(STDname).(Hname).(output_name).(vnames{ii}),...
                            'ColorLimits',ColorLimit,'ColorbarVisible','off','GridVisible','off',...
                            'Colormap',map,'CellLabelColor','none','MissingDataColor',[0 0 0],'FontSize',12)
                        xlabel(xlabel_name)
                        %title(vnames_label{ii})
                        FileName= strcat('fig_',vnames{ii});
                        fullFolderName = strcat(Folder,'\figures\smpl1970\',output_name,'\',TestType,'\',STDname,'\',data_freq,'_',dep1_name,'_',dep2_name,'\');
                     
                        mkdir(fullFolderName)
                        fullFileName = fullfile(fullFolderName, strcat(FileName,'.png'));
                        print(h,fullFileName,'-dpng')
                         fullFileName = fullfile(fullFolderName, strcat(FileName,'.eps'));
                         print(h,fullFileName,'-depsc')    
                    end
                    close all
                    
                    
                end
            end
            clear output output_fig
            
%%
% GFC Case Study
%%%%%%%

%1. Countries that underwent GFC
cnames_gfc = {'AUT','BEL','DNK','FRA','DEU','GRC','IRL','ITA','NLD','PRT','ESP','SWE','CHE','GBR','USA',...
            'HUN','LVA',...
            'RUS','SVN'};
        
%2. Loop
for jiv = 1:length(vnames) %Indicator
    temp1 = []; 
    temp2 = [];
   for jih = 1:length(depname) %Horizon
       temp = [];
       temp_color = [];
       for  jic = 1:length(cnames_gfc) %Countries
           if jic<=13 || jic >= 16
               idx = 155;
           elseif jic == 14
               idx = 151;
           elseif jic == 15
               idx = 152;
           end
           temp = [temp; save_gfc.(cnames_gfc{jic}).(depname{jih}).(vnames{jiv}).Y1hat(idx,1)];
           temp_color = [temp_color; save_gfc_color.(cnames_gfc{jic}).(depname{jih}).(vnames{jiv})];
       end
       temp1 = [temp1 temp]; 
       temp2 = [temp2 temp_color];
   end
           boxi.(vnames{jiv}) = temp1; %Horizon 1=0h, 21 = 20h
           boxi_color.(vnames{jiv}) = temp2;

end

width = 4.5 ; height = 2;
%3. Plot Scatterplot with colors
for jiv = 1:length(vnames)
    h = figure('Renderer', 'painters','units','inch', 'Position', [1 1 width height],'color','w')

   % boxplot(fliplr(boxi.(vnames{jiv})),'Labels',{'20','19','18','17','16','15','14','13','12','11','10','9','8','7','6','5','4','3','2','1','0'});
    hold on
    x = (repmat([1:21],19,1));
    y = fliplr(boxi.(vnames{jiv}));
    y_color =  fliplr(boxi_color.(vnames{jiv}));
    y_color_vec = y_color(:);
    

   
    y_color_rgb = [y_color(:) y_color(:) y_color(:)];
    y_color_rgb(y_color_vec == 0,:) = repmat([1 1 1],sum(y_color_vec == 0),1);
    y_color_rgb(y_color_vec == -0.5,:) = repmat([0.8 0.8 0.8],sum(y_color_vec == -0.5),1);
    y_color_rgb(y_color_vec == 2,:) = repmat([1 0 0],sum(y_color_vec == 2),1);
    y_color_rgb(y_color_vec == -2,:) = repmat([0 0 1],sum(y_color_vec == -2),1);
    y_color_rgb(y_color_vec == 1,:) = repmat([1 0.7 0.7],sum(y_color_vec == 1),1);
    y_color_rgb(y_color_vec == -1,:) = repmat([0.7 0.7 1],sum(y_color_vec == -1),1);



    
    scatter(x(:),y(:),50,y_color_rgb,'filled','MarkerEdgeColor',[0 0 0])
      b1 = plot(x(1,:),nanmedian(y),'-k','LineWidth',4);

      
    if jiv==1 || jiv == 7
        legend([b1],{'Median'},'Location','northwest')
    end
    ylim([0,1])
    xlim([1,21])
    xticks([1 5 9 13 17 21]);
    xticklabels({'20','16','12','8','4','0'})
    ylabel('Predicted probability')
    xlabel('Horizon')
    
    FileName= strcat('fig_',vnames{jiv});
    fullFolderName = strcat(Folder,'\figures\GFC_case\');
    mkdir(fullFolderName)
    fullFileName = fullfile(fullFolderName, strcat(FileName,'.png'));
    print(h,fullFileName,'-dpng')
    fullFileName = fullfile(fullFolderName, strcat(FileName,'.eps'));
    print(h,fullFileName,'-depsc')   
end
close all 
