%% Load Sensitivity matrix

modeIndex = [2]; % Indexes of these measured modes
n_modes = length(modeIndex); % Number of measured modes
%% Assemble structure matrices
N = length(Minit);

measDOFs = N/2+[6; 8; 10; 14];

Mzeroinit = Minit;

Mzeroinit(actDOFs(1),:) = [];
Mzeroinit(:,measDOFs(1)) = [];

Kzero = Ksolved;
% Kzero = Kinit;
Kzero(actDOFs(1),:) = [];
Kzero(:,measDOFs(1)) = [];

structModel.M0 = sparse(Mzeroinit);
structModel.K0 = sparse(Kzero);

Kdiffzero = Kdiff;
for n = 1:length(Kdiffzero)
    Kdiffzero{n}(actDOFs(1),:) = [];
    Kdiffzero{n}(:,measDOFs(1)) = [];
end
structModel.K_j = Kdiffzero;


%% Optimization structure parameter;
optimzOpts.tolFun = 1e-4;
optimzOpts.tolX = 1e-4;
optimzOpts.toolBox  = 'lsqnonlin';
optimzOpts.optAlgorithm = 'trust-region-reflective';
% optimzOpts.optAlgorithm = 'Levenberg-Marquardt';
optimzOpts.gradSel = 'on';
optimzOpts.maxIter = 50e3;
optimzOpts.maxFunEvals = 12e5;


%% Actual values for stiffness updating variables, each alpha represents 
% relative change from nominal stiffness parameter.
alpha_act = zeros(1,length(structModel.K_j));
n_alpha = length(alpha_act);


%% Simulate "experimental data"


lambdaExp = [130].^2;
psi_m = zeros(length(measDOFs),n_modes);
psi_m(2,:) = 1;

%%
expModes.lambdaExp = lambdaExp;
expModes.psiExp = psi_m;
expModes.measDOFs = measDOFs;
unmeasDOFs = setdiff(1 : N, measDOFs);
num_measDOFs = length(measDOFs);
num_unmeasDOFs = length(unmeasDOFs);

expModes.lambdaWeights = 100*ones(n_modes,1);

expModes.psiWeights = zeros(n_modes,1);

expModes.resWeights = zeros(n_modes,1);



%% Model updating parameter
updatingOpts.formID = 1;       % 1: Modal property diff (MAC) ;
                                % 2: Modal property diff (V_mDiff);
                                 % 3: Modal dynamic residual;
updatingOpts.modeMatch = 2;      % 1: Without forced matching;
                                 % 2: With forced matching;


updatingOpts.simModesForExpMatch = modeIndex;
if(updatingOpts.formID < 3)
    updatingOpts.x_lb = -2*ones(n_alpha,1);
    updatingOpts.x_ub =  10*ones(n_alpha,1);
else
    updatingOpts.x_lb = [-2*ones(n_alpha,1); -2* ones(num_unmeasDOFs * n_modes,1)];
    updatingOpts.x_ub =  [5*ones(n_alpha,1); 5 * ones(num_unmeasDOFs * n_modes,1)];
end

%% MultiStart optimization
numRuns = 1;
filename = ['TEST' num2str(updatingOpts.formID) '_JAC' optimzOpts.gradSel '_' optimzOpts.optAlgorithm '.mat'];

MultiRunModelUpdating

x_2 = x(1:n_alpha)

%%
%%
Ksolvedzero = full(structModel.K0);
for n = 1:n_alpha
Ksolvedzero = Ksolvedzero + x_2(n)*full(structModel.K_j{n});
end

structModel.K = Ksolvedzero;

%%

[psi, lambda] = eig(full(structModel.M0)\full(structModel.K0), "nobalance");

[lambda_init,dummyInd] = (sort(diag(lambda), 'ascend'));
zeros_init = sqrt(lambda_init);
zeros_init(1:2)
%%


    [psi, lambda] = eig(full(structModel.M0)\full(structModel.K), "nobalance");
    [lambda, ind] = sort( diag(lambda), 'ascend' );

zeros_solved = sqrt(lambda);
zeros_solved(1:2)