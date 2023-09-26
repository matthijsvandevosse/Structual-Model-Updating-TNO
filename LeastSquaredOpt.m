%% Add documentation
%
%
%
%
%
sensorloc.Model = [measDOFs(1) measDOFs(3) measDOFs(4)];
sensorloc.Data = [1 3 5];
freq = 10:6:28;
% freq = [15 18 20];
structModel.D0 = 4.8e-2;
n = 0;
clear structModel.alpha{n} structModel.Z{n}
for index = freq
    n = n+1;
    omega = data_opt.G_ref.Frequency(index)*2*pi;
    structModel.alpha{n} =  inv((-omega^2*structModel.M + structModel.K)/30000);
    structModel.Z{n} = (-omega^2*structModel.M + structModel.K)/30000;
end
%% 0.0045 - 2373.07
structModel.deltaS = [];
test = zeros(N,N);
loc = [0 0.125 0.25 0.375 0.5];
loc = [0 0.5];
for i = 1:length(loc)-1
nodes = [structuralmodelinit.Mesh.Nodes(1,:)  <= loc(i+1) & structuralmodelinit.Mesh.Nodes(1,:)  > loc(i)];
 structModel.deltaS{i} = diag([nodes nodes]);
 % structModel.deltaS{i-1+length(loc)} = diag([ nodes zeros(1,N/2)]);
end


%%
% Default optimization options if not provided
optimzOpts = struct ('maxIter', 400, 'maxFunEvals', 100, ...
    'tolFun', 1e-4, 'tolX', 1e-4, 'gradSel', 'on',...
    'optAlgorithm', 'Levenberg-Marquardt', 'x0', zeros(1,length(structModel.deltaS)));

% updatingOpts.x_lb = [0 + 0i];
% updatingOpts.x_ub = [1 + 1i];


optimzOpts.toolBox = 'lsqnonlin';

maxIterName = 'MaxIterations';
maxFunName = 'MaxFunctionEvaluations';

% If the bounds are not empty, MATLAB forces to use
% trust-region-reflective algorithm.
if strcmp(optimzOpts.optAlgorithm, 'Levenberg-Marquardt')
    updatingOpts.x_lb = [];
    updatingOpts.x_ub = [];
end


fun = @(x) OptmzObjJac_matthijs(x, data_opt, structModel, freq, sensorloc, optimzOpts.toolBox);

if (strcmp(optimzOpts.toolBox, 'lsqnonlin'))
    options = optimoptions( 'lsqnonlin', 'tolFun', optimzOpts.tolFun, 'tolX', optimzOpts.tolX,...
        'Algorithm', optimzOpts.optAlgorithm, 'Disp', 'iter', 'Jacobian', optimzOpts.gradSel,...
        maxIterName, optimzOpts.maxIter, maxFunName, optimzOpts.maxFunEvals );
    if strcmp(optimzOpts.optAlgorithm, 'Levenberg-Marquardt')
        updatingOpts.x_lb = [];
        updatingOpts.x_ub = [];
    end

    % If the bounds are not empty, MATLAB forces to use
    % trust-region-reflective algorithm.
    [updtRslts.xOpt, updtRslts.fvalOpt, updtRslts.residual, updtRslts.exitFlag, ...
        updtRslts.output, ~, S_temp] = lsqnonlin(fun, optimzOpts.x0, ...
        updatingOpts.x_lb, updatingOpts.x_ub, options );
    updtRslts.gradient = 2 * full(S_temp)' * updtRslts.residual;

    x = updtRslts.xOpt;


elseif (strcmp(optimzOpts.toolBox, 'fmincon'))
    nvar = length(optimzOpts.x0);
    if ( strcmp(optimzOpts.optAlgorithm,'interior-point') && nvar >= 1000 )
        prompt = ['Warning: The number of optimization variables is too large for\n' ...
            'interior point method in the fmincon toolbox. It is recommended to abort\n'...
            'and change to the trust-region-reflective algorithm in lsqnonlin toolbox.\n'...
            '\nWould you like to continue? Y/N [N]:'];
        userInput = input(prompt, 's');
        if isempty (userInput)
            userInput = 'N';
        end
        if( upper(userInput) ~= 'Y')
            error('Program terminates.');
        end
    end

    if(strcmp(optimzOpts.gradSel, 'on'))
        optimzOpts.gradSel = true;
    else
        optimzOpts.gradSel = false;
    end
    options = optimoptions( 'fmincon', 'tolFun', optimzOpts.tolFun, 'tolX', optimzOpts.tolX,...
        'Algorithm', optimzOpts.optAlgorithm, 'Disp', 'iter', 'SpecifyObjectiveGradient',optimzOpts.gradSel,...
        maxIterName, optimzOpts.maxIter, maxFunName, optimzOpts.maxFunEvals );
    [updtRslts.xOpt, updtRslts.fvalOpt, updtRslts.exitFlag, updtRslts.output,~,updtRslts.gradient] = fmincon( fun, optimzOpts.x0,[],[],[],[], ...
        updatingOpts.x_lb, updatingOpts.x_ub,[], options );
end

for alpha = 1:length(x)
structModel.D = structModel.D0*eye(N) + x(alpha).*structModel.deltaS{alpha};
end