clc
clear
close all
%% Information
% This script calculates the values of the objective function as functions
% of phi_1 and phi_2
load("DerivedFactors.mat", "DerivedFactors");
load("SensitivityAnalysis.mat", "SensitivityAnalysis");

% Define symbolic variables
GeoVars = sym('GeoVars', [1 4]);              %mu, phi_1, phi_2, r_frac

% Entered weights to be minimized. Set c2 to 0 if only xi is to be
% optimized
c_1     = SensitivityAnalysis.Weights.c1;
c_2     = SensitivityAnalysis.Weights.c2;

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

%% Plot these functions for r_frac is 0.9 and mu is 0.63
% First, enter the desired value of the amplification factor, the
% corresponding value of mu, the maximum r_fraction value and the minimum
% size of the shoe
normAmpFactor   = SensitivityAnalysis.Constraint.normXi;
normmu          = SensitivityAnalysis.Constraint.normMu;
r_frac_max      = SensitivityAnalysis.Constraint.r_fraction_max;
phi_mindiff     = SensitivityAnalysis.Constraint.phi_mindiff;

% Next, express r_frac as a function of phi_1 and phi_2
% The following line is used to extract mu from the DerivedFactors
% structure
[~, index_mu, ~]        = find(DerivedFactors.symbolicvariables == 'mu');

AmpFactorConstraint     = DerivedFactors.AmpFactor_nondim == normAmpFactor;
AmpFactorConstraint     = isolate(subs(AmpFactorConstraint, ...
    DerivedFactors.symbolicvariables(index_mu), normmu), 'r_frac');

% Enter number of points to be displayed along one axis and possible
% minimum shoe size
n           = 50;

phi_1_array     = linspace(0, pi, n);
phi_2_array     = linspace(0, pi, n);

[Phi_1_array, Phi_2_array]  = meshgrid(phi_1_array, phi_2_array);

r_fracfunc              = matlabFunction(rhs(AmpFactorConstraint));
R_frac                  = r_fracfunc(Phi_1_array, Phi_2_array);

Mu_array        = normmu .* ones(n,n);
R_frac_array    = r_frac_max .* ones(n,n);

% Calclate value of objective functions
for i = 1:n
    for j = 1:n
        if (R_frac(i,j) < 0 || R_frac(i,j) > r_frac_max || j >= i)
            R_frac(i,j)     = NaN;
            ObjFuncVal(i,j) = NaN;
            A_0(i,j)        = NaN;
            A_N(i,j)        = NaN;
            A_F(i,j)        = NaN;
            Alpha(i,j)      = NaN;
            Beta(i,j)       = NaN;
        elseif (phi_1_array(j) + phi_mindiff > phi_2_array(i))
            R_frac(i,j)     = NaN;
            ObjFuncVal(i,j) = NaN;
            A_0(i,j)        = NaN;
            A_N(i,j)        = NaN;
            A_F(i,j)        = NaN;
            Alpha(i,j)      = NaN;
            Beta(i,j)       = NaN;
        else
            ObjFuncVal(i,j)  = func_Var_combined_2Order([normmu, Phi_1_array(i,j), Phi_2_array(i,j), R_frac(i,j)]);
            [A_0(i,j), A_N(i,j), A_F(i,j)]  = ArmLengthCalc(normmu, Phi_1_array(i,j), ...
            Phi_2_array(i,j), R_frac(i,j), DerivedFactors, SensitivityAnalysis);
            [Alpha(i,j), Beta(i,j)]         = AlphaBeta(Phi_1_array(i,j), Phi_2_array(i,j), DerivedFactors);
        end
    end
end

%% Calculate the gradient of the objective function
gridsize          = Phi_1_array(1,2) - Phi_1_array(1,1);
[J_phi_1, J_phi_2] = gradient(ObjFuncVal, gridsize);

%% Plot surface and contour plot of R_frac (manifold of solutions)
figure(1)
surf(Phi_1_array, Phi_2_array, R_frac)
xlabel('$\phi_1$', 'Interpreter', 'Latex', 'FontSize', 18)
ylabel('$\phi_2$', 'Interpreter', 'Latex', 'FontSize', 18)
xlim([0 pi])
ylim([0 pi])
zlabel('$r_{frac}$', 'Interpreter', 'Latex', 'FontSize', 18)
colorbar

levels  = [0.5, 0.6, 0.7, 0.8, 0.9];
figure(2)
contour(Phi_1_array, Phi_2_array, R_frac, levels)

% Objective function value for mu = 0.9 and r_frac = 0.9
%title('Objective Function value','FontSize', 18)
xlabel('$\phi_1$', 'Interpreter', 'Latex', 'FontSize', 18)
ylabel('$\phi_2$', 'Interpreter', 'Latex', 'FontSize', 18)
xlim([0 pi])
ylim([0 pi])
colorbar

%% Plot arm lengths within the solution space
figure(3)
surf(Phi_1_array, Phi_2_array, A_0)
xlabel('$\phi_1$', 'Interpreter', 'Latex', 'FontSize', 18)
ylabel('$\phi_2$', 'Interpreter', 'Latex', 'FontSize', 18)
xlim([0 pi])
ylim([0 pi])
zlabel('$a_{0}$', 'Interpreter', 'Latex', 'FontSize', 18)
colorbar

figure(4)
surf(Phi_1_array, Phi_2_array, A_N)
xlabel('$\phi_1$', 'Interpreter', 'Latex', 'FontSize', 18)
ylabel('$\phi_2$', 'Interpreter', 'Latex', 'FontSize', 18)
xlim([0 pi])
ylim([0 pi])
zlabel('$a_{N}$', 'Interpreter', 'Latex', 'FontSize', 18)
colorbar

figure(5)
surf(Phi_1_array, Phi_2_array, A_F)
xlabel('$\phi_1$', 'Interpreter', 'Latex', 'FontSize', 18)
ylabel('$\phi_2$', 'Interpreter', 'Latex', 'FontSize', 18)
xlim([0 pi])
ylim([0 pi])
zlabel('$a_{F}$', 'Interpreter', 'Latex', 'FontSize', 18)
colorbar

%% Plot alpha and beta
figure(6)
surf(Phi_1_array, Phi_2_array, Alpha)
xlabel('$\phi_1$', 'Interpreter', 'Latex', 'FontSize', 18)
ylabel('$\phi_2$', 'Interpreter', 'Latex', 'FontSize', 18)
xlim([0 pi])
ylim([0 pi])
zlabel('$\alpha$', 'Interpreter', 'Latex', 'FontSize', 18)
colorbar

figure(7)
surf(Phi_1_array, Phi_2_array, Beta)
xlabel('$\phi_1$', 'Interpreter', 'Latex', 'FontSize', 18)
ylabel('$\phi_2$', 'Interpreter', 'Latex', 'FontSize', 18)
xlim([0 pi])
ylim([0 pi])
zlabel('$\beta$', 'Interpreter', 'Latex', 'FontSize', 18)
colorbar

figure(8)
surf(Phi_1_array, Phi_2_array, Alpha - Beta)
xlabel('$\phi_1$', 'Interpreter', 'Latex', 'FontSize', 18)
ylabel('$\phi_2$', 'Interpreter', 'Latex', 'FontSize', 18)
xlim([0 pi])
ylim([0 pi])
zlabel('$\alpha - \beta$', 'Interpreter', 'Latex', 'FontSize', 18)
colorbar

%% Plot Objective Function Value and derivatives
figure(9)
surf(Phi_1_array, Phi_2_array, ObjFuncVal)
xlabel('$\phi_1$', 'Interpreter', 'Latex', 'FontSize', 18)
ylabel('$\phi_2$', 'Interpreter', 'Latex', 'FontSize', 18)
xlim([0 pi])
ylim([0 pi])
zlabel('Objective Function', 'Interpreter', 'Latex', 'FontSize', 18)
colorbar

levels  = [0.118, 0.12, 0.13, 0.15, 0.18, 0.2, 0.25, 0.3, 0.4, 0.5];
figure(10)
contour(Phi_1_array, Phi_2_array, ObjFuncVal, levels)

% Objective function value for mu = 0.9 and r_frac = 0.9
title('Objective Function value', 'Interpreter', 'Latex', 'FontSize', 18)
xlabel('$\phi_1$', 'Interpreter', 'Latex', 'FontSize', 14)
ylabel('$\phi_2$', 'Interpreter', 'Latex', 'FontSize', 14)
xlim([0 pi])
ylim([0 pi])
colorbar

figure(11)
surf(Phi_1_array, Phi_2_array, J_phi_1)
xlabel('$\phi_1$', 'Interpreter', 'Latex', 'FontSize', 18)
ylabel('$\phi_2$', 'Interpreter', 'Latex', 'FontSize', 18)
xlim([0 pi])
ylim([0 pi])
zlabel('$D_{\phi_1} J(\phi_1, \phi_2)$', 'Interpreter', 'Latex', 'FontSize', 18)
colorbar

figure(12)
surf(Phi_1_array, Phi_2_array, J_phi_2)
xlabel('$\phi_1$', 'Interpreter', 'Latex', 'FontSize', 18)
ylabel('$\phi_2$', 'Interpreter', 'Latex', 'FontSize', 18)
xlim([0 pi])
ylim([0 pi])
zlabel('$D_{\phi_2} J(\phi_1, \phi_2)$', 'Interpreter', 'Latex', 'FontSize', 18)
colorbar

%% Function to calculate arm lengths
function [a_0, a_N, a_F] = ArmLengthCalc(mu, phi_1, phi_2, r_frac, ...
                                            DerivedFactors, SensitivityAnalysis)
    
    alpha_func = matlabFunction(DerivedFactors.alpha);
    beta_func  = matlabFunction(DerivedFactors.beta);
    
    a_0         = r_frac * SensitivityAnalysis.Parameters.r * cos(pi/2 - ...
                    alpha_func(phi_1, phi_2));
    a_N         = r_frac * SensitivityAnalysis.Parameters.r * cos(pi/2 - ...
                    beta_func(phi_1, phi_2));
    a_F         = SensitivityAnalysis.Parameters.r - r_frac * ...
                    SensitivityAnalysis.Parameters.r * sin(pi/2 - ...
                    beta_func(phi_1, phi_2));
end

%% Function to calculate alpha and beta
function [alpha, beta]  = AlphaBeta(phi_1, phi_2, DerivedFactors)
    % Define the functions used to calculate alpha and beta
    alpha_func = matlabFunction(DerivedFactors.alpha);
    beta_func  = matlabFunction(DerivedFactors.beta);
    
    alpha   = alpha_func(phi_1, phi_2);
    beta    = beta_func(phi_1, phi_2);
end
