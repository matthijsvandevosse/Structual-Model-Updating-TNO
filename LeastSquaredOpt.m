%% Add documentation
%
%
%
%
%
nModes = 50;
data_opt.G_ref = G_bla;

sensorloc.Model = [measDOFs(1) measDOFs(3) measDOFs(4)];
sensorloc.Data = [1 3 5];
freq = [10 15 20 25 30 ];% 40 50 100 200 300];
% freq = 12;
n = 0;
structModel.alpha = [];
structModel.Z = [];

[Psi_solved, lambdasolved] = eigs(structModel.K, structModel.M0,nModes,1e3);

[lambdasolved,dummyInd] = sort(diag(lambdasolved), 'ascend');

Psi_solved = Psi_solved(:,dummyInd);

for i = 1:size(Psi_solved,2)
    m = (Psi_solved(:,i))'*full(structModel.M0)*(Psi_solved(:,i));
    Psi_solved(:,i) = Psi_solved(:,i)./sqrt(m);
    Dm(i,i) = Psi_solved(:,i)'*structModel.D0*Psi_solved(:,i);
end


clear R G


for ii = freq
    n = n+1;
    w = data_opt.G_ref.Frequency(ii)*2*pi;
    for m = 1:nModes%length(lambda)
        R = Psi_solved(:,m)*Psi_solved(:,m)';
        if m == 1
             structModel.alpha{n} = gain* R./( -w^2 + 1i*Dm(m,m)*w + lambdasolved(m) );
        else
             structModel.alpha{n} =  structModel.alpha{n} + gain* R./( -w^2 + 1i*Dm(m,m)*w + lambdasolved(m) );  
        end
    end
        structModel.Z{n} = (-w^2*structModel.M + structModel.K)./gain;
end



%% 
strut_damping1 = 4e4;
strut_damping2 = 20e3;

strut_nodes1 = nodes(coord(:,1) > .25);
strut_nodes2 = nodes(coord(:,2)<0);

% structModel.D0(strut_nodes1,strut_nodes1) = strut_damping1*structModel.D0(strut_nodes1,strut_nodes1);
% structModel.D0(strut_nodes2,strut_nodes2) = strut_damping2*structModel.D0(strut_nodes2,strut_nodes2);

% structModel.deltaS{1} = sparse(0.*structModel.D0);
% structModel.deltaS{1}(strut_nodes1,strut_nodes1) = structModel.D0(strut_nodes1,strut_nodes1);
% structModel.deltaS{2} = sparse(0.*structModel.D0);
% structModel.deltaS{2}(strut_nodes2,strut_nodes2) = structModel.D0(strut_nodes2,strut_nodes2);
structModel.deltaS{1} = sparse(.1.*structModel.D0);



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