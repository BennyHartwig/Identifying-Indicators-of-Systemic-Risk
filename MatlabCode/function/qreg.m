function [ output] = qreg(y, X,tau, vce )

% check inputs
if nargin < 4
    vce.bwidth = 'hsheather';
    vce.method = 'iid';
    vce.type = 'HendricksKoenker';
    vce.level = 0.95;
end

% check missing inputs for vce
if isempty(vce.type), vce.type = 'HendricksKoenker'; end
if isempty(vce.bwidth), vce.bwidth = 'hsheather'; end
if isempty(vce.method), vce.method = 'iid'; end
if isempty(vce.level), vce.level = 0.95; end

T = length(y);
K = size(X,2);
% Perform quantile regression
[ beta ] = qregressMatlab( y, X, tau );
yhat= X*beta;
u =  y- yhat;
check_loss = tau-(u<0);
scale = 1/T*diag(u'* check_loss)';

% Density: location and scale form of asymmetric lapace
% f(u) = tau*(1-tau)*exp( -rho_tau[u]/sigma )/sigma
% where u = y-x*beta(tau)
% rho_tau(u) = u(tau-I(u<0)) is error check function
% and sigma =1/T*sum[rho_tau(u)] is ML estimate
logLik = T*log(tau.*(1-tau)) -T*log(scale) -T ;
% Information Criterion
AIC = -logLik + K;
BIC  = -logLik + 1/2*K*log(T);

%% Estimation of Bandwidth Parameters
alpha = vce.level; % confidence niveau
switch vce.bwidth
    case 'hsheather' % Hall and Sheather (1988)
        h = T^(-1/3)*abs(norminv((1-alpha)/2))^(2/3) * (3/2 * normpdf(norminv(tau))^2/(2*norminv(tau)^2+1))^(1/3);
        
    case 'bofinger' % Bofinger (1975)
        h = T^(-1/5)*((9/2*(normpdf(2*norminv(tau)))^4)/(2*norminv(tau)^2+1)^2 )^(1/5);
        
    case 'chamberlain' % Chamberlain (1994)
        h = abs(norminv((1-alpha)/2))*sqrt(tau*(1-tau)/T);
end

%% Estimation of Variance Covariance Matrix
% Hendricks and Koenker (1991)
tau_high = min(0.99,tau+h); tau_low = max(tau-h,0.01);
beta_high = qregressMatlab(y,X, tau_high);
beta_low = qregressMatlab(y,X, tau_low);

% Used in Powell (1986,1991), Kim and White (2003)
kappa = median(abs(u-median(u)));  % median absolute deviation of residual
c = kappa*(norminv(tau_high)-norminv(tau_low));


switch vce.method
    case 'iid'
        switch vce.type
            case 'HendricksKoenker' % Hendricks and Koenker (1991)
                f = (2*h)./(mean(X)*(beta_high-beta_low));
                covb = tau*(1-tau)*f^-2*(X'*X)^-1;
            case {'Powell','KimWhite'}
                f = 1/(2*T*c)*sum((abs(u)<c)); % empirical histogram
                covb = tau*(1-tau)*f^-2*(X'*X)^-1;

        end
        
        
    case 'robust'
        switch vce.type
            case 'HendricksKoenker' % Hendricks and Koenker (1991)
                e= 10^-5; % offset constant
                f = max(0,(2*h)./(X*(beta_high-beta_low) -e)); % avoid zero residuals
                D1 = (X'*sparse(1:T,1:T,f)*X);
                D0 = X'*X;
                covb = tau*(1-tau)*D1^-1*D0*D1^-1;

            case 'Powell' % Powell (1986,1991)
                f = 1/(2*c)* (abs(u)<c);
                D1 =X'*sparse(1:T,1:T,f)*X;
                D0 = X'*X;
                covb = tau*(1-tau)*D1^-1*D0*D1^-1;
            case 'KimWhite'    
                f = 1/(2*c)* (abs(u)<c);
                check_loss = (tau-(u<0));
                D1 =X'*sparse(1:T,1:T,f)*X;
                D0 = X'*sparse(1:T,1:T,check_loss.^2)*X;
                covb = D1^-1*D0*D1^-1;
        end
        
        
        
end
se = sqrt(diag(covb));






%% Store output
output.y = y;
output.X = X;
output.tau = tau;
output.vce = vce;

output.u = u;
output.yhat = yhat;
output.beta = beta;
output.scale  = scale;
output.covb = covb;
output.se = se;
output.t = beta./se;
output.pval = (1-tcdf(abs(beta./se),T-K))*2;
output.nobs = T;
output.K = K;
output.lik = logLik;
output.AIC = AIC;
output.BIC = BIC;
output.f = f;
output.h = h;
output.beta_high = beta_high;
output.beta_low = beta_low;
output.kappa = kappa;
output.c = c;