%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Identifiy Indicators of Systemic Risk (2020)
% Benny Hartwig, Christoph Meinering, Yves Schueler
% Deutsche Bundesbank
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ bhat ] = qregressMatlab( y, x, tau )
%   bhat are the estimates
%   y is a vector of outcomes
%   x is a matrix with columns of explanatory variables
%   tau is a scalar for choosing the conditional quantile to be estimated
options = optimoptions('linprog','Algorithm','dual-simplex','Display','none');
n=size(x,1);
m=size(x,2);
% vectors and matrices for linprog
f=[tau*ones(n,1);(1-tau)*ones(n,1);zeros(m,1)];
Aeq=[eye(n),-eye(n),x];
beq=y;
lb=[zeros(n,1);zeros(n,1);-inf*ones(m,1)];
ub=inf*ones(m+2*n,1);

% Solve the linear programme
try
[bhat,~,~]=linprog(f,[],[],Aeq,beq,lb,ub,options);
 catch
 bhat = nan(size(Aeq,2),1);
 end


% Pick out betas from (u,v,beta)-vector.
bhat=bhat(end-m+1:end);

end