%% Load Sensitivity matrix

modeIndex = [1 2 3 4 5]; % Indexes of these measured modes
n_modes = length(modeIndex); % Number of measured modes
%% Assemble structure matrices

structModel.M0 = sparse(Minit);
structModel.K0 = sparse(Kinit);

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


%% Simulate "experimental data"


load("eigenFrequencies.mat")
freqExp = eigenFrequencies(modeIndex);
freqExp(4) = eigenFrequencies(5);
freqExp(5) = 176;
lambdaExp = (2*pi*(freqExp)).^2;
load("Modeshapes_V1.mat")
L(:,4) = L(:,5);
for i = 1:modeIndex(end)
    [m,index] = max(abs(L([1 2 3 5],i)));
    L(:,i) = L(:,i) / m;
end


psiExpAll = zeros(length(measDOFs),1);
n = 1;
for i = modeIndex
psiExpAll([1], n) = L(1,i);
psiExpAll([2], n) = L(2,i);
psiExpAll([3], n) = L(3,i);
psiExpAll([4], n) = L(5,i);
% psiExpAll([5], n) = 0.9999*L(1,i);
% psiExpAll([6], n) = 0.9999*L(2,i);
% psiExpAll([7], n) = 0.9999*L(3,i);
% psiExpAll([8], n) = 0.9999*L(5,i);
n = n+1;
end

psi_m = psiExpAll;

%%
expModes.lambdaExp = lambdaExp;
expModes.psiExp = psi_m;
expModes.measDOFs = measDOFs;
unmeasDOFs = setdiff(1 : N, measDOFs);
num_measDOFs = length(measDOFs);
num_unmeasDOFs = length(unmeasDOFs);
% expModes.lambdaWeights = [1 1 1];
expModes.lambdaWeights = [20 20 20 20 15];
% expModes.lambdaWeights = 2*[10 10 10 10];
% expModes.psiWeights = [1 1 1];
expModes.psiWeights = ones(n_modes,1);

expModes.psiWeights = [1000 200 400 200 0.1];
% expModes.psiWeights = [1000 200 400 200];
expModes.resWeights = ones(n_modes,1);



%% Model updating parameter
updatingOpts.formID = 1;       % 1: Modal property diff (MAC) ;
                                % 2: Modal property diff (V_mDiff);
                                 % 3: Modal dynamic residual;
updatingOpts.modeMatch = 2;      % 1: Without forced matching;
                                 % 2: With forced matching;
updatingOpts.simModesForExpMatch = modeIndex;
if(updatingOpts.formID < 3)
    updatingOpts.x_lb = -1.5*ones(n_alpha,1);
    updatingOpts.x_ub =  1.5*ones(n_alpha,1);
    
else
    updatingOpts.x_lb = [-2*ones(n_alpha,1); -2* ones(num_unmeasDOFs * n_modes,1)];
    updatingOpts.x_ub =  [5*ones(n_alpha,1); 5 * ones(num_unmeasDOFs * n_modes,1)];
end

%% MultiStart optimization
numRuns = 1;
randSeed = round(100*rand(1));
filename = ['TEST' num2str(updatingOpts.formID) '_JAC' optimzOpts.gradSel '_' optimzOpts.optAlgorithm '.mat'];

MultiRunModelUpdating

x_2 = x(1:n_alpha)

Ksolved = full(structModel.K0);
for n = 1:n_alpha
Ksolved = Ksolved + x_2(n)*full(structModel.K_j{n});
end

structModel.K = Ksolved;

[Psi_solvedxy, lambdasolved] = eig(full(structModel.M0)\Ksolved, "nobalance");
 
[Psi_o, lambda_o] = eig(full(structModel.M0)\structModel.K0);