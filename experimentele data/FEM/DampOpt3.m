D_new = 4e-2*eye(N,N);
x = 1;
%%
n = 0;
m = 0;
clear S b D_new_param omega f
%%
stop = false;
while ~stop
    if m >3
        if mean(abs(diff(f(m-3:end))))<0.01
            stop = true;
        end
    end
    m = m+1;
    n = 0;
for index = [14]
    n=n+1;
omega = G_ref.Frequency(index)*2*pi;
alpha = inv((-omega^2*structModel.M0 + Ksolved)/30000);
alpha_damped = inv((-omega^2*structModel.M0 + omega*D_new*1i + Ksolved)/30000);
% Z = (-omega^2*structModel.M0 + Ksolved )/30000;
% Zdamped = (-omega^2*structModel.M0 + omega*D_new*i + Ksolved)/30000;
alpha_damped_real = G_ref.ResponseData([1 3 5],:,index);

% alpha_undamped_real =  real(alpha_damped_real) + imag(alpha_damped_real)*1./(real(alpha_damped_real))*imag(alpha_damped_real);
alpha_undamped_real = alpha([measDOFs(1) measDOFs(3) measDOFs(4)],[measDOFs(1) measDOFs(3) measDOFs(4)]);
% alpha([measDOFs(1) measDOFs(3) measDOFs(4)],[measDOFs(1) measDOFs(3) measDOFs(4)]);
% alpha_damped([measDOFs(1) measDOFs(3) measDOFs(4)],[measDOFs(1) measDOFs(3) measDOFs(4)])


b1 = 1i*(alpha_damped_real - alpha_undamped_real);
b2 = alpha*(omega/30000*D_new)*alpha_damped;
b2([measDOFs(1) measDOFs(3) measDOFs(4)],[measDOFs(1) measDOFs(3) measDOFs(4)]);


b(:,:,n) = b1 - b2([measDOFs(1) measDOFs(3) measDOFs(4)],[measDOFs(1) measDOFs(3) measDOFs(4)]);

S_full = alpha*(omega/30000*eye(N,N))*alpha_damped;
S(:,:,n) = S_full([measDOFs(1) measDOFs(3) measDOFs(4)],[measDOFs(1) measDOFs(3) measDOFs(4)]);
f(n) = abs(b(1,1,n));
xr(n) = (real(b(1,1,n)))./(real(S(1,1,n)));
xi(n) = (imag(b(1,1,n)))./(imag(S(1,1,n)));
x(n) = mean([xr(n) xi(n)]);
end
[fmax,MaxIndex] = max(f)
mean(f)
f = (f/fmax).^2;
x_step = mean(x.*f);

D_new = (0.5*(x_step))*eye(N,N)+D_new;

D_new_param(m) = D_new(end,end)
end

D_new = D_new_param(end)*eye(N);

