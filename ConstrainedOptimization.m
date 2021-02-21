clc
clear
close all
%% Script Information
% This script attempts to apply constrained optimization, to find the
% combination of parameters that is 'most robust'. This is done by
% minimizing the weighted 'spread' in both variables (mu*q_2 and xi), where
% the radii a and b are understood to be the square roots of the eigenvalue
% matrix

% Make use of the first and second order taylor approximations

%% Apply constrained optimization
load("DerivedFactors.mat", "DerivedFactors");
load("SensitivityAnalysis.mat", "SensitivityAnalysis");

% Define symbolic variables
GeoVars = sym('GeoVars', [1 4]);              %mu, phi_1, phi_2, r_frac

% Entered weights to be minimized. Set c2 to 0 if only xi is to be
% optimized
c_1     = SensitivityAnalysis.Weights.c1;
c_2     = SensitivityAnalysis.Weights.c2;

% Enforced conditions
global normAmpFactor AmpFactorfunc AmpFactorfunc_2Order
phi_mindiff     = SensitivityAnalysis.Constraint.phi_mindiff;
r_fraction_max  = SensitivityAnalysis.Constraint.r_fraction_max;
normAmpFactor   = SensitivityAnalysis.Constraint.normXi;

% Starting point for constrained optimization. May not be zero
StartPoint      = SensitivityAnalysis.Optimization.startstate;     %mu, phi1, phi2, r_frac

%% Determine the function used for the constraints
% Determine the first order estimate for the variance as a function of
% the variances
[~, index_mu, ~]        = find(DerivedFactors.symbolicvariables == 'mu');
[~, index_phi_1, ~]     = find(DerivedFactors.symbolicvariables == 'phi_1');
[~, index_phi_2, ~]     = find(DerivedFactors.symbolicvariables == 'phi_2');
[~, index_r_frac, ~]    = find(DerivedFactors.symbolicvariables == 'r_frac');

mu      = DerivedFactors.symbolicvariables(1, index_mu);
phi_1   = DerivedFactors.symbolicvariables(1, index_phi_1);
phi_2   = DerivedFactors.symbolicvariables(1, index_phi_2);
r_frac  = DerivedFactors.symbolicvariables(1, index_r_frac);

fmin_amp    = subs(DerivedFactors.AmpFactor_nondim, [mu, phi_1, phi_2, r_frac], GeoVars);
fmin_muQ    = subs(mu * DerivedFactors.q_1_nondim, [mu, phi_1, phi_2, r_frac], GeoVars); 

AmpFactorfunc       = matlabFunction(fmin_amp, 'Vars', {GeoVars});
q_1func             = matlabFunction(fmin_muQ, 'Vars', {GeoVars});

%% Determine the first and second order taylor approximations of xi and mu * q_1
% Assuming the 4 i.i.d. input variables result in a two dimensional vector
% output, calculate the 2x2 covariance matrix of a first order
% approximation

% First, import the variables necessary
ArmVars = sym('ArmVars', [1 4]);              %mu, a_0, a_N, a_F

[~, index_a_0, ~]       = find(DerivedFactors.symbolicvariables == 'a_0');
[~, index_a_N, ~]       = find(DerivedFactors.symbolicvariables == 'a_N');
[~, index_a_F, ~]       = find(DerivedFactors.symbolicvariables == 'a_F');

a_0     = DerivedFactors.symbolicvariables(1, index_a_0);
a_N     = DerivedFactors.symbolicvariables(1, index_a_N);
a_F     = DerivedFactors.symbolicvariables(1, index_a_F);

% Define the function to be optimized
f_optimAmpFull      = subs(DerivedFactors.AmpFactorArms, [mu, a_0, a_N, a_F], ArmVars);
f_optimMuQ_1Full    = subs(mu*DerivedFactors.q_1Arms, [mu, a_0, a_N, a_F], ArmVars);

% Determine the average values of all arms for a given geometric vector x
% Use this average value to determine the standard deviation
mu_bar          = mu;
mu_bar          = subs(mu_bar, [mu, phi_1, phi_2, r_frac], GeoVars);
sigma_mu        = (SensitivityAnalysis.Variance.fraction_mu / 3) * mu_bar;

a_0_bar         = r_frac * SensitivityAnalysis.Parameters.r * cos(pi/2 - ...
    DerivedFactors.alpha);
a_0_bar         = subs(a_0_bar, [mu, phi_1, phi_2, r_frac], GeoVars);
sigma_a_0       = (SensitivityAnalysis.Variance.fraction_a_0 / 3) * a_0_bar;

a_N_bar         = r_frac * SensitivityAnalysis.Parameters.r * cos(pi/2 - ...
    DerivedFactors.beta);
a_N_bar         = subs(a_N_bar, [mu, phi_1, phi_2, r_frac], GeoVars);
sigma_a_N       = (SensitivityAnalysis.Variance.fraction_a_N / 3) * a_N_bar;

a_F_bar         = SensitivityAnalysis.Parameters.r - ...
    r_frac * SensitivityAnalysis.Parameters.r * sin(pi/2 - DerivedFactors.beta);
a_F_bar         = subs(a_F_bar, [mu, phi_1, phi_2, r_frac], GeoVars);
sigma_a_F       = (SensitivityAnalysis.Variance.fraction_a_F / 3) * a_F_bar;

% Combine all standard deviations in one matrix. Perturbations are IID,
% hence the diagonal matrix
sigmaA          = diag([sigma_mu, sigma_a_0, sigma_a_N, sigma_a_F]).^2;

% Evaluate jacobian at a_0_bar, a_N_bar, a_F_bar, mu_bar, which are
% expressed in terms of geometric variables. Thus, the geometric variables
% are used to calculate the according nominal arms, whose values are used
% to compute the power of the noise accordingly. The variance and average
% of xi and mu * q_2 in terms of arm variables can thus be brought back to
% a function of the geometric variables.
Jacobian        = subs(jacobian([f_optimAmpFull, f_optimMuQ_1Full], ArmVars), ...
    ArmVars, [mu_bar, a_0_bar, a_N_bar, a_F_bar]);
Variance_1Order = Jacobian * sigmaA * ...
    transpose(Jacobian);
Average_1Order  = subs([f_optimAmpFull; f_optimMuQ_1Full], ArmVars, ...
    [mu_bar, a_0_bar, a_N_bar, a_F_bar]);

% Continue and calculate the second order variance. First, calculate the
% required hessians at a_0_bar, a_N_bar, a_F_bar, mu_bar
Hessian_1       = subs(hessian(f_optimAmpFull, ArmVars), ArmVars, ...
    [mu_bar, a_0_bar, a_N_bar, a_F_bar]);
Hessian_2       = subs(hessian(f_optimMuQ_1Full, ArmVars), ArmVars, ...
    [mu_bar, a_0_bar, a_N_bar, a_F_bar]);

Variance_2Order = Variance_1Order + 1/2 .* ...
    [trace(sigmaA * Hessian_1 * sigmaA * Hessian_1), ...
    trace(sigmaA * Hessian_1 * sigmaA * Hessian_2); ...
    trace(sigmaA * Hessian_2 * sigmaA * Hessian_1), ...
    trace(sigmaA * Hessian_2 * sigmaA * Hessian_2)];

Average_2Order = Average_1Order + 1/2 .* [trace(Hessian_1 * sigmaA); ...
    trace(Hessian_2 * sigmaA)];

AmpFactorfunc_2Order    = matlabFunction(Average_2Order(1,1), 'Vars', {GeoVars});
q_1func_2Order          = matlabFunction(Average_2Order(2,1), 'Vars', {GeoVars});

% First order with incorporation of the covariance
func_Var_combined_1Order    = matlabFunction(c_1^2 * Variance_1Order(1,1) + ...
    2 * c_1 * c_2 * Variance_1Order(1,2) + c_2^2 * Variance_1Order(2,2),'Vars', {GeoVars});

func_Var_combined_2Order    = matlabFunction(c_1^2 * Variance_2Order(1,1) + ...
    2 * c_1 * c_2 * Variance_2Order(1,2) + c_2^2 * Variance_2Order(2,2),'Vars', {GeoVars});

%% Set constraint matrices for optimization
A = [   [0, 1,-1, 0]
        [0, 0, 1, 0]
        [0,-1, 0, 0]
        [0, 0, 0,-1]
        [0, 0, 0, 1]];
b = [-phi_mindiff; pi; 0; 0; r_fraction_max];

Aeq = [1, 0, 0, 0];
beq = [SensitivityAnalysis.Constraint.normMu];

lb = [];
ub = [];

%% Minimise variances
nonlcon                         = @local_nonlcon_1order;
options                         = optimoptions('fmincon'); 
options.ConstraintTolerance     = 1e-10;
options.StepTolerance           = 1e-10;
options.OptimalityTolerance     = 1e-10;
options.Algorithm              = 'interior-point';
disp("First Order Taylor Approximation of Variations")
[optim, funValue]       = fmincon(func_Var_combined_1Order, StartPoint, ...
    A, b, Aeq, beq, lb, ub, nonlcon, options)

% Calculate the amplification factor, the value of q_1 and the value of 
% the linearised variances accompanying the found point
AmplificationFactor_1Order      = AmpFactorfunc(optim)
q_1_1Order                      = q_1func(optim)
q_1_check                       = vpa(subs(fmin_muQ, GeoVars, optim), 5)
Variance_xi                     = vpa(subs(Variance_1Order(1,1), GeoVars, optim), 5)
Covariance_xi_muq_1             = vpa(subs(Variance_1Order(2,1), GeoVars, optim), 5)
Variance_muq_1                  = vpa(subs(Variance_1Order(2,2), GeoVars, optim), 5)
a_0_bar_1                       = vpa(subs(a_0_bar, GeoVars, optim), 5)
a_N_bar_1                       = vpa(subs(a_N_bar, GeoVars, optim), 5)
a_F_bar_1                       = vpa(subs(a_F_bar, GeoVars, optim), 5)
optim_1                         = optim;

% Repeat the procedure with the second order approximation
% Set the start point to equal the solution of the first order
% approximation
disp("Second Order Taylor Approximation of Variations")
nonlcon                         = @local_nonlcon_2order;
[optim, funValue]       = fmincon(func_Var_combined_2Order, StartPoint, ...
    A, b, Aeq, beq, lb, ub, nonlcon, options)

AmplificationFactor_2Order      = AmpFactorfunc_2Order(optim)
q_1_2Order                      = q_1func_2Order(optim)
Variance_xi                     = vpa(subs(Variance_2Order(1,1), GeoVars, optim), 5)
Covariance_xi_muq_1             = vpa(subs(Variance_2Order(2,1), GeoVars, optim), 5)
Variance_muq_1                  = vpa(subs(Variance_2Order(2,2), GeoVars, optim), 5)
a_0_bar_2                       = vpa(subs(a_0_bar, GeoVars, optim), 5)
a_N_bar_2                       = vpa(subs(a_N_bar, GeoVars, optim), 5)
a_F_bar_2                       = vpa(subs(a_F_bar, GeoVars, optim), 5)
Average_correction_2Order       = vpa(subs(1/2 .* [trace(Hessian_1 * sigmaA); ...
    trace(Hessian_2 * sigmaA)], GeoVars, optim), 5)
optim_2                         = optim;

%% Declaration of nonlinear constraint
function [c, ceq] = local_nonlcon_1order(x)
    global normAmpFactor AmpFactorfunc AmpFactorfunc_2Order
    c = [];
    ceq = normAmpFactor - AmpFactorfunc(x);
end

function [c, ceq] = local_nonlcon_2order(x)
    global normAmpFactor AmpFactorfunc AmpFactorfunc_2Order
    c = [];
    ceq = normAmpFactor - AmpFactorfunc_2Order(x);
end