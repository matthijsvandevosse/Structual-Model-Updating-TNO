%% Load Sensitivity matrix

modeIndex = [8; 16; 18; 21; 28]; % Indexes of these measured modes
modeIndex = [8; 16; 18; 28]; % Indexes of these measured modes
n_modes = length(modeIndex); % Number of measured modes
%% Assemble structure matrices

structModel.M0 = sparse(Minit);
structModel.K0 = sparse(Kinit);
structModel.M_j = Mdiff;
structModel.M_j = [];
structModel.K_j = Kdiff;


%% Actual values for stiffness updating variables, each alpha represents 
% relative change from nominal stiffness parameter.
alpha_act = zeros(1,length(structModel.K_j));
n_alpha = length(alpha_act) + length(structModel.M_j) ;


%% Optimization structure parameter;
optimzOpts.tolFun = 1e-5;
optimzOpts.tolX = 1e-6;
optimzOpts.toolBox  = 'lsqnonlin';
optimzOpts.optAlgorithm = 'trust-region-reflective';
% optimzOpts.optAlgorithm = 'Levenberg-Marquardt';
optimzOpts.gradSel = 'on';
optimzOpts.maxIter = 50e3;
optimzOpts.maxFunEvals = 12e5;


%% Simulate "experimental data" with old experiments


load("G_ss_Modal_fitted_4.mat")

lambdaExp = diag(-G_ss_Modal.A(size(G_ss_Modal.A)/2+1:end,1:size(G_ss_Modal.A)/2));
freqExp =sqrt(lambdaExp)/2/pi;

%conversion to rad/sec
lambdaExp = (freqExp*2*pi*2*pi).^2

L = G_ss_Modal.C(:,1:size(G_ss_Modal.A)/2);
R = G_ss_Modal.B(size(G_ss_Modal.B)/2+1:end,:);


%%

psi_m = [R(1:3,:); R(5,:)]';

lambdaExp = [lambdaExp(1:3); lambdaExp(5)];
%%
expModes.lambdaExp = lambdaExp(1:n_modes);
expModes.psiExp = psi_m(:,1:n_modes);
expModes.measDOFs = measDOFs;

unmeasDOFs = setdiff(1 : N, measDOFs(:,1));
num_measDOFs = length(measDOFs);
num_unmeasDOFs = length(unmeasDOFs);



expModes.lambdaWeights = 10*[1 1 1 1];
expModes.lambdaWeights = expModes.lambdaWeights(1:n_modes);

expModes.psiWeights = ones(n_modes,1);

expModes.psiWeights = [1 1 1 1 ];
expModes.psiWeights = expModes.psiWeights(1:n_modes);

expModes.resWeights = ones(n_modes,1);

expModes.Actuator_pos = Actuator_pos(:,1:2);


%% Model updating parameter
updatingOpts.formID = 1.4;       % 1: Modal property diff (MAC) ;
                                % 2: Modal property diff (V_mDiff);
                                 % 3: Modal dynamic residual;
                                 % 1.3: Zero Matching
                                 % 1.4: Delta Psi Matching

updatingOpts.modeMatch = 2;      % 1: Without forced matching;
                                 % 2: With forced matching;
updatingOpts.simModesForExpMatch = modeIndex;

if(updatingOpts.formID < 3)
    updatingOpts.x_lb = -5*ones(n_alpha,1);
    updatingOpts.x_ub =  5*ones(n_alpha,1);
    
else
    updatingOpts.x_lb = [-2*ones(n_alpha,1); -2* ones(num_unmeasDOFs * n_modes,1)];
    updatingOpts.x_ub =  [5*ones(n_alpha,1); 5 * ones(num_unmeasDOFs * n_modes,1)];
end

%% Start optimization
numRuns = 1;
randSeed = round(100*rand(1));
filename = ['TEST' num2str(updatingOpts.formID) '_JAC' optimzOpts.gradSel '_' optimzOpts.optAlgorithm '.mat'];


n_x = length(updatingOpts.x_lb);

% zero initialization
optimzOpts.x0 = zeros(n_x, 1);

% random initialization
% optimzOpts.x0 = 0.05*rand(n_x, 1);


updtResults = StructModelUpdating(structModel, expModes, updatingOpts, optimzOpts);
x = updtResults.xOpt;
fval = updtResults.fvalOpt;
optmzSolvOutput = updtResults.output;
    
%%
x_k = (x(1:n_alpha))
% x_m = x(n_alpha+1,end)
%%
structModel.K = structModel.K0;
for i = 1 : n_alpha
    structModel.K = structModel.K + x(i) * structModel.K_j{i};
end

structModel.M = structModel.M0;
% for i = 1 : n_beta
%     structModel.M = structModel.M + x(i + n_alpha) * structModel.M_j{i};
% end


[psi_solved, lambda_solved] = eigs(structModel.K,  structModel.M, max(modeIndex),  (160*2*pi)^2, 'IsSymmetricDefinite', 1);
[psi_o, lambda_o] = eigs(structModel.K0,  structModel.M0, max(modeIndex),  (160*2*pi)^2, 'IsSymmetricDefinite', 1);

