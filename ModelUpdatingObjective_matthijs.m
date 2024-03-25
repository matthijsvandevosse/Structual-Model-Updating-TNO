function [S, b] =  ModelUpdatingObjective_matthijs(x, data, structModel,freq, sensorloc)
b = zeros(18*length(freq),1);
S = zeros(18*length(freq),1);
n = 0;
N = size(structModel.M);

for alpha = 1:length(x)
structModel.D = structModel.D0*eye(N) + x(alpha).*structModel.deltaS{alpha};
end

x

for index = freq
    n = n+1;
    %% Model response


    omega = data.G_ref.Frequency(index)*2*pi;
    Zdamped = (-omega^2*structModel.M + omega*(structModel.D)*1i + structModel.K)/800;

    %% Measured response
    alpha_damped_real = data.G_ref.ResponseData(sensorloc.Data,:,index);
    alpha_undamped_real = structModel.alpha{n}(sensorloc.Model,sensorloc.Model);

    b1 = 1i*(alpha_damped_real - alpha_undamped_real);
    b2 = structModel.Z{n}\(omega/800*(structModel.D))/Zdamped;
 
    b_temp = b1 - b2(sensorloc.Model,sensorloc.Model);

    for j = 1:length(x)
    S_full = -structModel.Z{n}\(omega/800*structModel.deltaS{j})/Zdamped;
    S_temp = S_full(sensorloc.Model,sensorloc.Model);
    Sr = reshape(real(S_temp), [] ,1);
    Si = reshape(imag(S_temp), [] ,1); 
    S(1+(n-1)*18:n*18,j) = [Sr; Si];
    end
    %%
    br = reshape(real(b_temp), [] ,1);
    bi = reshape(imag(b_temp), [] ,1);
    b(1+(n-1)*18:n*18,1) = [br; bi];


end
end