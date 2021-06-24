%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Identifiy Indicators of Systemic Risk (2020)
% Benny Hartwig, Christoph Meinering, Yves Schueler
% Deutsche Bundesbank
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function results = ML_2Step(...
    dep1,mat1st,regname,indexStruc,...
    dep2,mat2nd,...
    opt)

% Extract options
pMIN = opt.pMIN;
pMAX = opt.pMAX;
horizon = opt.horizon;
tau = opt.tau;
method = opt.method;
startIND = opt.startIND;
endIND = opt.endIND;


Tfull = length(dep1);
for p = pMIN : pMAX
    for m = 0 : 4
        if m == 0
            indexREG2 = [];
        elseif m > 0
            
            indexREG2 = [  1:m];
        end
        Cname = strcat('m',num2str(m),'_p',num2str(p));
        for lagg = 0 : horizon
            lag = lagg;
            Hname = strcat('F',num2str(lag));
            results.(Cname).(Hname) = [];
            startIND = opt.startIND+lag+p;
            
            clear sfields
            output = struct;
            parfor ll = 1 : length(regname)
                indexREG1 = indexStruc.(regname{ll});
                
                % REGRESSION: \Delta(1) y(t+h) = y(t+h) - y(t+h-1), where y(t) is in log
                Hname = strcat('F',num2str(lag));
                y1 = dep1;
                x1 = [lagmatrix(mat1st(:,indexREG1),lag) ones(Tfull,1)];
                y2 = dep2;
                if lag == 0
                    x2 = [ ones(Tfull,1) mat2nd(:,indexREG2) ];
                else
                    x2 = [ ones(Tfull,1) lagmatrix(mat2nd(:,indexREG2),lag-1) ];
                end
                startIND =max([ sum(isnan(y1),1) sum(isnan(x1),1) sum(isnan(y2),1)  sum(isnan(x2),1)])+1;
                % check for all 1's or all 0's
                tmp = find(y1(startIND:end) ==1);
                chk = length(tmp);
                [nobs junk] = size(y1(startIND:end));
                % make NaN regression: 
                if startIND > endIND; startIND = 1 ; % when there is no value for the indicator
                elseif chk == nobs || chk == 0; startIND = 1 ;  % when there is no crises in the sample
                end
                
                %                 Switch Estimation Method
                switch method
                    case 'LogitOLS'
                        output = ML_2step_LogitOLS( y1, x1,y2,x2, pMIN,p,startIND,endIND);
                    case 'LogitQReg'
                        output= ML_2step_LogitQReg( y1, x1,y2,x2, pMIN,p,tau,startIND,endIND);
                        
                end
                %output.Y1 = [];
                output.X1 = [];
                %output.Y1hat = [];
                output.U1 = [];
                output.Y2 = [];
                output.X2 = [];
                output.Z = [];
                output.Y2hat = [];
                output.U2 = [];

                sfields{ll,1}= output;
                
                
                %                 results.(Cname).(Hname).(regname{ll}) = output;
                
                
                
                
                
            end
            for ll = 1 : length(regname)
                results.(Cname).(Hname).(regname{ll}) = sfields{ll,1};
            end
        end
    end
end

results.INFO.dep1 = dep1;
results.INFO.mat1st = mat1st;
results.INFO.regname = regname;
results.INFO.indexStruc = indexStruc;
results.INFO.dep2 = dep2;
results.INFO.mat2nd = mat2nd;
results.INFO.opt = opt;