function [b, S] = OptmzObjJac_matthijs(x, data, structModel, freq, sensorloc, optToolBox)
% Documentation!!


% Output:
%	r: the objective residual vector r(x)
%   jac: the Jacobian matrix of the residual vector d_r / d_x


[S,b] =  ModelUpdatingObjective_matthijs(x, data, structModel,freq, sensorloc);
b = sparse(b);
S = sparse(S);
x
% fmincon use scalar as objective function output
if(strcmp(optToolBox,'fmincon'))
    if nargout > 1
    	S = 2 * S' * b;
    end
    b = norm(b)^2;
end

end

