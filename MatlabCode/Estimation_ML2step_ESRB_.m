%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Identifiy Indicators of Systemic Risk (2020)
% Benny Hartwig, Christoph Meinering, Yves Schueler
% Deutsche Bundesbank
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear all; clc;

% Start parallel pool
parfevalOnAll(gcp(), @warning, 0, 'off')

%Prelim
addpath(genpath('function'));
spl = '1970';
data_freq = 'data'; % 
dep1_name = 'ESRB';% 
dep2_name = 'd1lRGDP'; %

            
            % Specify settings
            logit_type = 'Binary'; % #Binary #Ordered
            
            
            load(data_freq)

            switch spl
                case '1970'; [tab ] = tabcut(tab,1,80); % start 1970
                case '1985'; [tab ] = tabcut(tab,1,140); % start 1985
            end
            mkdir(['results/smpl',spl])
            vnames = {'cqgap','cqgap_ham','Mend_cyc','Mend_cyc_1600','Mend_cyc_rt',...
                'reer','cabgdp','spreads',...
                'FCycle','FCyc_dreh','FCyc_dreh_rt',...
                'd1lrtcredit','d4lrtcredit','d8lrtcredit','d12lrtcredit',...
                'd1lRPPI','d4lRPPI','d8lRPPI','d12lRPPI',...
                'd1lrstock','d4lrstock','d8lrstock','d12lrstock',...
                'd1lcabgdp','d4lcabgdp','d8lcabgdp','d12lcabgdp',...
                'd1lreer','d4lreer','d8lreer','d12lreer',...
                'd1lrbondpr','d4lrbondpr'};
            horizon= 20; % forecast horizon %20
            nlag = 12; % maximum number of lagged regressors %12
            startIND = 1;  % 1970/Q1 
            endIND = 188; %2017Q4 
            cnames = tab.(dep1_name).Properties.VariableNames(2:end);
            idx = max(table2array(tab.(dep1_name)(:,2:end)));
            cnames(isnan(idx)) = [];
            
            
          
            
            
            for ff =  1 : length(cnames)
                country = cnames{ff};
                
                dep1 = tab.(dep1_name).(country); dep1(dep1>0) = 1;
                dep2 = tab.(dep2_name).(country);
                
                
                mat1st = [];
                INDEX = [1 : nlag];
                for i = 1 : length(vnames)
                    
                    for lag = 0 : nlag
                        mat_.(vnames{i})(:,lag+1) = lagmatrix(tab.(vnames{i}).(country) , lag) ;
                        indexStruc.(strcat('c_',vnames{i},'_',num2str(lag))) = [ (1+(length(INDEX)+1)*(i-1) :lag+1+(length(INDEX)+1)*(i-1)) ];
                        
                    end
                    mat1st =[  mat1st ,  mat_.(vnames{i}) ];
                    
                end
                
                regname = fieldnames(indexStruc);
                
                mat2nd = lagmatrix(dep2,[1:4]);
                
                % set up options
                opt.tau = 0.05;
                opt.horizon = horizon;
                opt.pMIN= 0;
                opt.pMAX = 0;
                opt.startIND = startIND;
                opt.endIND = endIND;
                
                display(country)
                opt.method = 'LogitOLS';
                %%  1. Stage Estimation
                display('Run 2 Step ML: LogitOLS')
                mname = 'LogitOLS';
                results.(mname) = ML_2Step(dep1,mat1st,regname,indexStruc,dep2,mat2nd,opt);
                eval(['results_',mname, ' =  results.(mname);'] );
                save(strcat('results/smpl',spl,'/results_',country,'_',mname,'_',dep1_name,'_',dep2_name),strcat('results_',mname),'indexStruc')
                clear('results',strcat('results_',mname))
                
                
                opt.method = 'LogitQReg';
                tau = [ 0.05];  %specificy quantile(s) 
                tau_pct = round(tau*100,0);
                for ii = 1 : length(tau)
                    mname =strcat('LogitQReg',sprintf('%02d',tau_pct(ii)));
                    display(strcat('Run 2 Step ML: LogitQReg',sprintf('%02d',tau_pct(ii))))
                    opt.tau = tau(ii);
                    results.(mname) = ML_2Step(dep1,mat1st,regname,indexStruc,dep2,mat2nd,opt);
                    eval(['results_',mname, ' =  results.(mname);'] );
                    save(strcat('results/smpl',spl,'/results_',country,'_',mname,'_',dep1_name,'_',dep2_name),strcat('results_',mname),'indexStruc')
                    clear('results',strcat('results_',mname))

                end
                
                
                clear mat_ INDEX  mat1st results indexStruc;
            end


