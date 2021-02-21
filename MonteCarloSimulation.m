clear; clc; close all;
%% Load factors and variables
load('DerivedFactors.mat', 'DerivedFactors');
load("SensitivityAnalysis.mat", "SensitivityAnalysis");

% Define geometric variables
GeoVars = sym('GeoVars', [1 4]);              %mu, phi_1, phi_2, r_frac

[~, index_mu, ~]        = find(DerivedFactors.symbolicvariables == 'mu');
[~, index_phi_1, ~]     = find(DerivedFactors.symbolicvariables == 'phi_1');
[~, index_phi_2, ~]     = find(DerivedFactors.symbolicvariables == 'phi_2');
[~, index_r_frac, ~]    = find(DerivedFactors.symbolicvariables == 'r_frac');

mu      = DerivedFactors.symbolicvariables(1, index_mu);
phi_1   = DerivedFactors.symbolicvariables(1, index_phi_1);
phi_2   = DerivedFactors.symbolicvariables(1, index_phi_2);
r_frac  = DerivedFactors.symbolicvariables(1, index_r_frac);

% Define the arm variables
ArmVars = sym('ArmVars', [1 4]);              %mu, a_0, a_N, a_F

[~, index_a_0, ~]       = find(DerivedFactors.symbolicvariables == 'a_0');
[~, index_a_N, ~]       = find(DerivedFactors.symbolicvariables == 'a_N');
[~, index_a_F, ~]       = find(DerivedFactors.symbolicvariables == 'a_F');

a_0     = DerivedFactors.symbolicvariables(1, index_a_0);
a_N     = DerivedFactors.symbolicvariables(1, index_a_N);
a_F     = DerivedFactors.symbolicvariables(1, index_a_F);

debug_plot  = true

%% Settings
nc  = 1e5;     %number of simulations per geometry

% Manual addition of configuration of interest. Set boolean to true if
% configuration should be examined
mu_manual       = SensitivityAnalysis.ManualConfig.mu;
phi_1_manual    = SensitivityAnalysis.ManualConfig.phi_1;
phi_2_manual    = SensitivityAnalysis.ManualConfig.phi_2;
r_frac_manual   = SensitivityAnalysis.ManualConfig.r_fraction;

% Group this information in a vector
geo_manual      = [mu_manual, phi_1_manual, phi_2_manual, r_frac_manual];

% Calculate the variances of this configuration by calculating the arm
% length. First calculate the value of alpha and beta
alpha_def       = DerivedFactors.alpha;
alpha           = subs(alpha_def, [mu, phi_1, phi_2, r_frac], geo_manual);
beta_def        = DerivedFactors.beta;
beta            = subs(beta_def, [mu, phi_1, phi_2, r_frac], geo_manual);

mu_bar          = mu_manual;
sigma_mu        = (SensitivityAnalysis.Variance.fraction_mu / 3) * mu_bar;

a_0_bar         = double(r_frac_manual * SensitivityAnalysis.Parameters.r * cos(pi/2 - alpha));
sigma_a_0       = (SensitivityAnalysis.Variance.fraction_a_0 / 3) * a_0_bar;

a_N_bar         = double(r_frac_manual * SensitivityAnalysis.Parameters.r * cos(pi/2 - beta));
sigma_a_N       = (SensitivityAnalysis.Variance.fraction_a_N / 3) * a_N_bar;

a_F_bar         = double(SensitivityAnalysis.Parameters.r - ...
    r_frac_manual * SensitivityAnalysis.Parameters.r * sin(pi/2 - beta));
sigma_a_F       = (SensitivityAnalysis.Variance.fraction_a_F / 3) * a_F_bar;

% Combine all standard deviations in one matrix. Perturbations are IID,
% hence the diagonal matrix
sigmaA          = diag([sigma_mu, sigma_a_0, sigma_a_N, sigma_a_F]).^2;

%% Take the viable configuration and add noise
% Define the function to calculate the amplification factor and mu * q_1
AmpFacFunc  = matlabFunction(subs(DerivedFactors.AmpFactorArms, ...
    [mu, a_0, a_N, a_F], ArmVars));
MuQ1Func    = matlabFunction(subs(mu*DerivedFactors.q_1Arms, ...
    [mu, a_0, a_N, a_F], ArmVars));

dialogueWindow  = waitbar(0, ...
    "Monte Carlo Simulation - Amplification Factor and mu q_1 (with covariance)");

% For the manually entered configuration
local_config.mu         = mu_manual;
local_config.a_0        = a_0_bar;
local_config.a_N        = a_N_bar;
local_config.a_F        = a_F_bar;

for i = 1:nc
    waitbar((i / nc), dialogueWindow, ...
    "Monte Carlo Simulation - Amplification Factor and mu q_1 (with covariance)")
    % For the amount of simulations per geometry
    % Noise with the specified variance is added to all four factors
    noisyConfig.mu      = Nnoise(local_config.mu, sigma_mu);
    noisyConfig.a_0     = Nnoise(local_config.a_0, sigma_a_0);
    noisyConfig.a_N     = Nnoise(local_config.a_N, sigma_a_N);
    noisyConfig.a_F     = Nnoise(local_config.a_F, sigma_a_F);

    %These 'noisy' values are used to calculate the 'noisy'
    %amplification factor and value for mu * q_1. 
    
    % WARNING - for speed purposes matlab functions are used. This means
    % that MuQ1Func only requires 3 inputs (all but a_0)
    local_config.results(:,i) = [AmpFacFunc(noisyConfig.mu, noisyConfig.a_0, ...
        noisyConfig.a_N, noisyConfig.a_F); MuQ1Func(noisyConfig.mu, ...
        noisyConfig.a_N, noisyConfig.a_F)];
end

% Close waitbar
close(dialogueWindow)

% Calculate covariance matrix manually. First calculate average
average     = mean(local_config.results, 2)

% Calculate the sample covariance, by definition
configs.covariance_matrix(1,1) = 1 / (nc - 1) * ...
    dot(local_config.results(1,:) - average(1,1), ...
    transpose(local_config.results(1,:) - average(1,1)));
configs.covariance_matrix(1,2) = 1 / (nc - 1) * ...
    dot(local_config.results(1,:) - average(1,1), ...
    transpose(local_config.results(2,:) - average(2,1)));
configs.covariance_matrix(2,1) = 1 / (nc - 1) * ...
    dot(local_config.results(2,:) - average(2,1), ...
    transpose(local_config.results(1,:) - average(1,1)));
configs.covariance_matrix(2,2) = 1 / (nc - 1) * ...
    dot(local_config.results(2,:) - average(2,1), ...
    transpose(local_config.results(2,:) - average(2,1)));

disp(join(["Covariance Matrix after", string(nc), "simulations"]))
configs.covariance_matrix

if debug_plot
    scatter(local_config.results(2,:), local_config.results(1,:), 10, 'black', '.')    
    title(join(["Monte Carlo Distribution (N = ", string(nc), ...
    ")"],""))
    grid('on')
    xlabel("$\mu \, q_1$", 'Interpreter', 'Latex')
    ylabel("$\xi$", 'Interpreter', 'Latex')
end

%% Noise function
function [delta] = Nnoise(mean, stddev)
delta = mean + stddev*randn; 
end