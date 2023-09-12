function [S, b] =  ModelUpdatingObjective_matthijs(x, data, structModel,freq, sensorloc)
b = zeros(18*length(freq),1);
S = zeros(18*length(freq),1);
n = 0;
    N = size(structModel.M0);

for alpha = 1:length(x)
structModel.D = structModel.D0*eye(N) + x(alpha).*structModel.deltaS{alpha};
end

for index = freq
    n = n+1;
    %% Model response


    omega = data.G_ref.Frequency(index)*2*pi;
    % alpha_damped = inv((-omega^2*structModel.M0 + omega*(structModel.D0+x)*eye(N)*1i + structModel.K)/30000);
    Zdamped = (-omega^2*structModel.M0 + omega*(structModel.D)*1i + structModel.K)/30000;
    % alpha_undamped_real =  real(alpha_damped_real) + imag(alpha_damped_real)*1./(real(alpha_damped_real))*imag(alpha_damped_real);

    %% Measured response
    alpha_damped_real = data.G_ref.ResponseData(sensorloc.Data,:,index);
    alpha_undamped_real = structModel.alpha{n}(sensorloc.Model,sensorloc.Model);

    b1 = 1i*(alpha_damped_real - alpha_undamped_real);
    b2 = structModel.Z{n}\(omega/30000*(structModel.D))/Zdamped;
    % b2 = structModel.alpha{n}*(omega/30000*(structModel.D0+x)*eye(N))*alpha_damped;

    b_temp = b1 - b2(sensorloc.Model,sensorloc.Model);

    for j = 1:length(x)
    S_full = -structModel.Z{n}\(omega/30000*structModel.deltaS{j})/Zdamped;
    S_temp = S_full(sensorloc.Model,sensorloc.Model);
    Sr = reshape(real(S_temp), [] ,1);
    Si = reshape(imag(S_temp), [] ,1); 
    S(1+(n-1)*18:n*18,j) = [Sr; Si];
    end
    %%
    br = reshape(real(b_temp), [] ,1);
    bi = reshape(imag(b_temp), [] ,1);
    b(1+(n-1)*18:n*18,1) = [br; bi];


    %% complex damping
    % b_temp = reshape((b_temp), [] ,1);
    % S_temp = reshape((S_temp), [] ,1);
    % 
    % b(1+(n-1)*9:n*9,1) = [b_temp];
    % S(1+(n-1)*9:n*9,1) = [S_temp];
end
end