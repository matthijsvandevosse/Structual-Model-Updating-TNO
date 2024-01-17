%% Load Sensitivity matrix

modeIndex = [1 2 3 4 5 6]; % Indexes of these measured modes
n_modes = length(modeIndex); % Number of measured modes
%% Assemble structure matrices

structModel.M0 = sparse(Minit);
structModel.K0 = sparse(Kinit);
structModel.M_j = [];
structModel.K_j = Kdiff;


%% Actual values for stiffness updating variables, each alpha represents 
% relative change from nominal stiffness parameter.
alpha_act = zeros(1,length(structModel.K_j));
n_alpha = length(alpha_act);


%% Optimization structure parameter;
optimzOpts.tolFun = 1e-3;
optimzOpts.tolX = 1e-5;
optimzOpts.toolBox  = 'lsqnonlin';
optimzOpts.optAlgorithm = 'trust-region-reflective';
% optimzOpts.optAlgorithm = 'Levenberg-Marquardt';
optimzOpts.gradSel = 'on';
optimzOpts.maxIter = 50e3;
optimzOpts.maxFunEvals = 12e5;


%% Simulate "experimental data" with old experiments
% 
% 
% load("eigenFrequencies.mat")
% freqExp = eigenFrequencies;
% freqExp(4) = eigenFrequencies(5);
% freqExp(5) = 176;
% freqExp(6) = 176;
% lambdaExp = (2*pi*(freqExp)).^2;
% load("Modeshapes_V1.mat")
% L(:,4) = L(:,5);
% L(:,6) = L(:,5);
% for i = 1:modeIndex(end)
%     [m,index] = max(abs(L([1 2 3 5],i)));
%     L(:,i) = L(:,i) / m;
% end
% 
% 
% psiExpAll = zeros(length(measDOFs),1);
% n = 1;
% for i = modeIndex
% psiExpAll([1], n) = L(1,i);
% psiExpAll([2], n) = L(2,i);
% psiExpAll([3], n) = L(3,i);
% psiExpAll([4], n) = L(5,i);
% n = n+1;
% end
% 
% psi_m = psiExpAll;
% Simulate "experimental data" with new experiment data
addpath("/Users/matthijsvandevosse/FRF_Balk_Processing/FRFs")

load("G_ss_modal_fitted_2.mat")

lambdaExp = diag(-G_ss_Modal.A(9:end,1:8));
freqExp =sqrt(lambdaExp)/2/pi;

L = G_ss_Modal.C(:,1:8);





for i = 1:modeIndex(end)
    [m,index] = max(abs(L([1 2 3 5],i)));
    L(:,i) = L(:,i) / m;
end

% L = flip(L,1);

psiExpAll = zeros(length(measDOFs),1);
n = 1;
for i = modeIndex
psiExpAll([1], n) = L(1,i);
psiExpAll([2], n) = L(2,i);
psiExpAll([3], n) = L(3,i);
psiExpAll([4], n) = L(5,i);
n = n+1;
end

% psi_m = [ psiExpAll(1,:) + psiExpAll(4,:); psiExpAll(2,:) - psiExpAll(3,:)];


psi_m = [psiExpAll(1,:); (psiExpAll(1,:) - psiExpAll(4,:)); (psiExpAll(2,:) - psiExpAll(3,:))];
%%
expModes.lambdaExp = lambdaExp(1:n_modes);
expModes.psiExp = psi_m(:,1:n_modes);
expModes.measDOFs = measDOFs;

unmeasDOFs = setdiff(1 : N, measDOFs(:,1));
num_measDOFs = length(measDOFs);
num_unmeasDOFs = length(unmeasDOFs);

% expModes.lambdaWeights = [1 1 1];
expModes.lambdaWeights = [30 30 20 20  15];
expModes.lambdaWeights = [50 50 30 30 0 10];
expModes.lambdaWeights = expModes.lambdaWeights(1:n_modes);
% expModes.lambdaWeights = 2*[10 10 10 10];
% expModes.psiWeights = [1 1 1];
expModes.psiWeights = ones(n_modes,1);

expModes.psiWeights = [1000 200 400 200 0.1 10];
expModes.psiWeights = [2000 400 400 200 0 100 ];
expModes.psiWeights = expModes.psiWeights(1:n_modes);
% expModes.psiWeights = [1000 200 400 200];
expModes.resWeights = ones(n_modes,1);



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
    updatingOpts.x_lb = -2*ones(n_alpha,1);
    updatingOpts.x_ub =  2*ones(n_alpha,1);
    
else
    updatingOpts.x_lb = [-2*ones(n_alpha,1); -2* ones(num_unmeasDOFs * n_modes,1)];
    updatingOpts.x_ub =  [5*ones(n_alpha,1); 5 * ones(num_unmeasDOFs * n_modes,1)];
end

%% Start optimization
numRuns = 1;
randSeed = round(100*rand(1));
filename = ['TEST' num2str(updatingOpts.formID) '_JAC' optimzOpts.gradSel '_' optimzOpts.optAlgorithm '.mat'];


n_x = length(updatingOpts.x_lb);
optimzOpts.x0 = zeros(n_x, 1);


updtResults = StructModelUpdating(structModel, expModes, updatingOpts, optimzOpts);
x = updtResults.xOpt;
fval = updtResults.fvalOpt;
optmzSolvOutput = updtResults.output;
    
%%
x_2 = x(1:n_alpha)
x_2(2) = abs(x_2(2))

Ksolved = (structModel.K0);
for n = 1:n_alpha
Ksolved = Ksolved + x_2(n)*(structModel.K_j{n});
end

structModel.K = Ksolved;

[Psi_solved, lambdasolved] = eigs(Ksolved, structModel.M0, 10, 1e3);

[Psi_o, lambda_o] = eigs(structModel.K0, structModel.M0, 10, 1e3);