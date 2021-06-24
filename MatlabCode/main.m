%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Identifiy Indicators of Systemic Risk (2020)
% Benny Hartwig, Christoph Meinering, Yves Schueler
% Deutsche Bundesbank
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; clc

%%%%%%%%%%%%%%%%%%%%%%%%
%% BENCHMARK ANALYSIS %%
%%%%%%%%%%%%%%%%%%%%%%%%


%% Estimation of the two step model using Maximum Likelihood
Estimation_ML2step_benchmark
% Estimates 2 Step model
% - OLS regression
% - Quantile regression for tau = [0.05]


%% Evaluation of results 
%% Create Figure
% Construct Heatmap and GFC case study
Figure_results2heatmap_benchmark
% Builds heatmap for three decision rules
% (1) Hierachical
% - 1st Stage no pass: white
% - 1st Stage pass: grey2
% - 2nd Stage pass: overwrite only grey by either redish or blueish
%
% (2) Joint
% - 2nd Stage no pass: white
% - 2nd Stage pass: redish or blueish
%
% (3) JointCounterfactual
% - 1st Stage no pass: white
% - 1st Stage pass: grey
% - 2nd Stage pass: overwrite all by green
% - 2nd Stage pass: overwrite all grey by either redish or blueish
%
% Color coding: 
% - redish indicates sum of coefficients on 1st stage is positive
% - bluish indicates sum of coefficients on 1st stage is negative


%% Create Table with standard errors
Table_results2xls_mainresult_benchmark
% Table format (header)
% - Type: OLS, Qreg
% - Indicator: Basel III credit-to-GDP gap, ...
% - 1st Stage: Sign
% - 1st Stage: p-val significance (stars)
% - 2nd Stage: Coefficient 
% - 2nd Stage: Standard Error
% - 2nd Stage: p-val significance (stars)


%% Comparison of standard errors (naiv s.e. vs 2Step s.e. under independece vs  2Step s.e. with cross term)
Table_results2xls_robustness_benchmark
% Table format (header)
% - Type: OLS (naiv s.e.), OLS (2Step s.e. under independence), OLS (2Step s.e. with cross terms), ...
% - Indicator: Credit growth, ...
% - 1st Stage: Sign
% - 1st Stage: p-val significance (stars)
% - 2nd Stage: Coefficient 
% - 2nd Stage: Standard Error
% - 2nd Stage: p-val significance (stars)
% - 2nd Stage: confidence interval lower bound
% - 2nd Stage: confidence interval upper bound



%%%%%%%%%%%%%%%%%%%%%%%%
%% Subsample Analysis %%
%%%%%%%%%%%%%%%%%%%%%%%%
Estimation_ML2step_benchmark_subsample
Figure_results2heatmap_benchmark_subsample

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Alternative financial disruptions %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%RomerRomer
Estimation_ML2step_RomerRomer_
Figure_results2heatmap_RomerRomer_

%European financial crises database
Estimation_ML2step_ESRB_
Figure_results2heatmap_ESRB_


