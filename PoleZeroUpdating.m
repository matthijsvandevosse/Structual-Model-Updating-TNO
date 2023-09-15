%% Load Sensitivity matrix

modeIndex = [1 2 3 4]; % Indexes of these measured modes
n_modes = length(modeIndex); % Number of measured modes
%% Assemble structure matrices

structModel.M0 = sparse(Minit);
structModel.K0 = sparse(Kinit);

structModel.K_j = Kdiff;
structModel.M_j = 0.05*structModel.M_0;

%% Actual values for stiffness updating variables, each alpha represents 
% relative change from nominal stiffness parameter.
alpha_act = zeros(1,length(structModel.K_j)+length(structModel.M_j));
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
% freqExp(5) = 170;
lambdaExp = (2*pi*(freqExp)).^2;
load("Modeshapes_V1.mat")
L(:,4) = L(:,5);
for i = 1:n_modes
    [m,index] = max(abs(L([1 2 3 5],i)));
    L(:,i) = L(:,i) / m;
end


psiExpAll = zeros(length(measDOFs),1);
n = 1;
for i = 1:n_modes
psiExpAll([1], n) = L(1,i);
psiExpAll([2], n) = L(2,i);
psiExpAll([3], n) = L(3,i);
psiExpAll([4], n) = L(5,i);
n = n+1;
end

psi_m = psiExpAll;
%% Zeros
clear modeIndexZeros
clear lambdaExpZeros
clear zerosDofs

lambdaExpZeros{1} = [125].^2;
modeIndexZeros{1} = [ 2];
zerosDofs{1} = [measDOFs(3), measDOFs(3)];

lambdaExpZeros{2} = [181 ].^2;
modeIndexZeros{2} = [ 2];
zerosDofs{2} = [measDOFs(1), measDOFs(1)];

lambdaExpZeros{3} = [187].^2;
modeIndexZeros{3} = [ 2];

% lambdaExpZeros{1} = [32.3 ].^2;
% modeIndexZeros{1} = [1 ];
% zerosDofs{1} = [measDOFs(3), measDOFs(3)];
% 
% lambdaExpZeros{2} = [ ].^2;
% modeIndexZeros{2} = [1 ];
% zerosDofs{2} = [measDOFs(1), measDOFs(1)];
% 
% lambdaExpZeros{3} = [21.2 ].^2;
% modeIndexZeros{3} = [1 ];


zerosDofs{3} = [measDOFs(4), measDOFs(4)];
%%
expModes.lambdaExp = lambdaExp;
expModes.lambdaExpZeros = lambdaExpZeros;
expModes.zerosDofs = zerosDofs;
expModes.modeIndexZeros = modeIndexZeros;
expModes.lambdaWeightsZeros{1} = [5];
expModes.lambdaWeightsZeros{2} = [25];
expModes.lambdaWeightsZeros{3} = [5];
expModes.psiExp = psi_m;
expModes.measDOFs = measDOFs;
unmeasDOFs = setdiff(1 : N, measDOFs);
num_measDOFs = length(measDOFs);
num_unmeasDOFs = length(unmeasDOFs);

expModes.lambdaWeights = [20 20 30 15];
% expModes.lambdaWeights = zeros(1,4);
expModes.psiWeights = ones(n_modes,1);

expModes.psiWeights = [500 200 400 200];
% expModes.psiWeights = zeros(1,4);


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

x_1 = x(1:n_alpha)
%%
Ksolved = full(structModel.K0);
for n = 1:n_alpha
Ksolved = Ksolved + x_1(n)*full(structModel.K_j{n});
end

structModel.K = Ksolved;

[Psi_solvedxy, lambdasolved] = eig(full(structModel.M0)\Ksolved, "nobalance");
 
[Psi_o, lambda_o] = eig(full(structModel.M0)\structModel.K0);