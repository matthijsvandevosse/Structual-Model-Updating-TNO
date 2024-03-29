clear structModel
%% Load Sensitivity matrix

modeIndex = [1:8]; % Indexes of these measured modes
n_modes = length(modeIndex); % Number of measured modes
%% Assemble structure matrices

structModel.M0 = sparse(Minit);
structModel.K0 = sparse(Kinit);

structModel.M_j{1} = 0.1*structModel.M0;

% structModel.M_j = [];
structModel.K_j = Kdiff;


%% Actual values for stiffness updating variables, each alpha represents 
% relative change from nominal stiffness parameter.
alpha_act = zeros(1,length(structModel.K_j));
n_alpha = length(alpha_act);
n_beta = length(structModel.M_j);

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
eigenFrequencies(6) = 290;
eigenFrequencies(7) = 438;
eigenFrequencies(8) = 643;
freqExp = eigenFrequencies(modeIndex);
freqExp(4) = eigenFrequencies(5);
freqExp(5) = 172;
lambdaExp = (2*pi*(freqExp)).^2;
load("Modeshapes_V1.mat")
L(:,4) = L(:,5);
L(:,5) =   [ -0.0482   -0.5369    1.0000   -0.6948   -0.2233];
L(:,6) =   [   -0.5244    0.5093    0.1579   -0.5494    1.0000];
L(:,7) =   [   -0.5918    1.0000   -0.7506    0.9519   -0.8093];
L(:,8) =    [   1.0000   -0.9733    0.3060    0.6149   -0.8124];

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

%% Zeros
clear modeIndexZeros
clear lambdaExpZeros
clear zerosDofs

lambdaExpZeros{1} = [130].^2;
modeIndexZeros{1} = [ 2];
zerosDofs{1} = [measDOFs(3), measDOFs(3)];

% 
lambdaExpZeros{2} = [181 ].^2;
modeIndexZeros{2} = [ 2];
zerosDofs{2} = [measDOFs(1), measDOFs(1)];


lambdaExpZeros{3} = [187].^2;
modeIndexZeros{3} = [2];
zerosDofs{3} = [measDOFs(4), measDOFs(4)];
% % 





%%
expModes.lambdaExpZeros = lambdaExpZeros;
expModes.zerosDofs = zerosDofs;
expModes.modeIndexZeros = modeIndexZeros;
expModes.lambdaWeightsZeros{1} = [1];
expModes.lambdaWeightsZeros{2} = [1];
expModes.lambdaWeightsZeros{3} = [1];
% expModes.lambdaWeightsZeros{1} = [0];
% expModes.lambdaWeightsZeros{2} = [0];
% expModes.lambdaWeightsZeros{3} = [0];

expModes.lambdaExp = lambdaExp;
expModes.psiExp = psi_m;
expModes.measDOFs = measDOFs;
unmeasDOFs = setdiff(1 : N, measDOFs);
num_measDOFs = length(measDOFs);
num_unmeasDOFs = length(unmeasDOFs);
% expModes.lambdaWeights = [1 1 1];
expModes.lambdaWeights = [20 20 20 20 10 10 10 10];
% expModes.lambdaWeights = 2*[10 10 10 10];
% expModes.psiWeights = [1 1 1];
expModes.psiWeights = ones(n_modes,1);

expModes.psiWeights = [1000 200 400 200 50 50 50 50];

% expModes.psiWeights = [1000 200 400 200];
expModes.resWeights = ones(n_modes,1);



%% Model updating parameter
updatingOpts.formID = 1.9;       % 1: Modal property diff (MAC) ;
                                % 2: Modal property diff (V_mDiff);
                                 % 3: Modal dynamic residual;
                                 % 1.9: Wiht zero matching
updatingOpts.modeMatch = 2;      % 1: Without forced matching;
                                 % 2: With forced matching;
updatingOpts.simModesForExpMatch = modeIndex;
if(updatingOpts.formID < 3)
    updatingOpts.x_lb = -1*ones(n_alpha+n_beta,1);
    updatingOpts.x_ub =  2*ones(n_alpha+n_beta,1);
else
    updatingOpts.x_lb = [-2*ones(n_alpha,1); -2* ones(num_unmeasDOFs * n_modes,1)];
    updatingOpts.x_ub =  [5*ones(n_alpha,1); 5 * ones(num_unmeasDOFs * n_modes,1)];
end


%% MultiStart optimization
numRuns = 1;
randSeed = round(100*rand(1));
filename = ['TEST' num2str(updatingOpts.formID) '_JAC' optimzOpts.gradSel '_' optimzOpts.optAlgorithm '.mat'];

MultiRunModelUpdating

x_1_K = x(1:n_alpha)
x_1_M = x(n_alpha+1:n_alpha+n_beta)
%%
Ksolved = full(structModel.K0);
Msolved = full(structModel.M0);
for n = 1:n_alpha
Ksolved = Ksolved + x_1_K(n)*full(structModel.K_j{n});
end
for n = 1:n_beta
Msolved = Msolved + x_1_M(n)*full(structModel.M_j{n});
end

structModel.K = Ksolved;
structModel.M = Msolved;

[Psi_solvedxy, lambdasolved] = eig(full(structModel.M)\Ksolved, "nobalance");
 
[Psi_o, lambda_o] = eig(full(structModel.M0)\structModel.K0);