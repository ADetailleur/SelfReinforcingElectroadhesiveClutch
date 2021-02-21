clc
clear
close all

%% Script Information
% This script stores all of variables related to the sensitivity analysis
% in one struct

% Define the absolute radius, needed to determine variance introduced in
% constrained optimization using arms
SensitivityAnalysis.Parameters.r            = 0.06;     %m

% Fraction of the arm length that is used to calculate a 3 sigma interval
SensitivityAnalysis.Variance.fraction_a_0   = 0.1;
SensitivityAnalysis.Variance.fraction_a_N   = 0.1;
SensitivityAnalysis.Variance.fraction_a_F   = 0.1;
SensitivityAnalysis.Variance.fraction_mu    = 0.1;

% Weights to examine [c1, c2]^T * Var([xi, mu*q_2]) * [c1, c2]
SensitivityAnalysis.Weights.c1  = 1;
SensitivityAnalysis.Weights.c2  = 0;

% Constraints to take into account
SensitivityAnalysis.Constraint.phi_mindiff      = (2 * pi) / 10;
SensitivityAnalysis.Constraint.r_fraction_max   = 0.925;
SensitivityAnalysis.Constraint.normXi           = 2.7;
SensitivityAnalysis.Constraint.normMu           = 0.63;
% note - normXi = normAmpFactor

% Manual configuration to examine
SensitivityAnalysis.ManualConfig.mu          = 0.63;
SensitivityAnalysis.ManualConfig.phi_1       = 1.1684;
SensitivityAnalysis.ManualConfig.phi_2       = 1.7968;
SensitivityAnalysis.ManualConfig.r_fraction  = 0.925;

% Starting configuratoin for Optimization
SensitivityAnalysis.Optimization.startstate     = ...
    [0.63, pi/4, 3*pi/4, 0.75];            %mu, phi1, phi2, r_frac

% Confidence interval for covariance ellipse
SensitivityAnalysis.ConfidenceInterval      = sqrt(9.21); 
% The eucledian distance from the average, after transformation by the
% inverse covariance matrix is given by a 2 DOF chi-squared distribution.
% Thus, this distribution should be used to determine the confidence
% interval
% sqrt(5.991) = 95% confidence ellipse
% sqrt(9.21)  = 99% confidence ellipse

% Save the configured struct
save('SensitivityAnalysis.mat', 'SensitivityAnalysis')
disp("Sensitivity Analysis Structure Successfully Updated")