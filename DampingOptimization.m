D_new = D;
x = 1;
clear S b D_new_param f
%%
n = 0;
for index = [14:19]
    n=n+1;
    x=1;
stop = false;
m = 0;
while ~stop
    m = m+1;
    if m >3
        mean(diff(f(n,m-3:end)))
        if mean(diff(f(n,m-3:end)))>-0.07
            stop = true;
        end
    end
omega(n) = G_ref.Frequency(index)*2*pi;
alpha = inv((-omega(n)^2*structModel.M0 + Ksolved)/30000);
alpha_damped = inv((-omega(n)^2*structModel.M0 + omega(n)*D_new*1i + Ksolved)/30000);
% Z = (-omega^2*structModel.M0 + Ksolved )/30000;
% Zdamped = (-omega^2*structModel.M0 + omega*D_new*i + Ksolved)/30000;
alpha_damped_real = G_ref.ResponseData([1 3 5],:,index);

alpha_undamped_real =  real(alpha_damped_real) + imag(alpha_damped_real)*1./(real(alpha_damped_real))*imag(alpha_damped_real);
alpha_undamped_real = alpha([measDOFs(1) measDOFs(3) measDOFs(4)],[measDOFs(1) measDOFs(3) measDOFs(4)]);
% alpha([measDOFs(1) measDOFs(3) measDOFs(4)],[measDOFs(1) measDOFs(3) measDOFs(4)]);
% alpha_damped([measDOFs(1) measDOFs(3) measDOFs(4)],[measDOFs(1) measDOFs(3) measDOFs(4)])


b1 = 1i*(alpha_damped_real - alpha_undamped_real);
b2 = alpha*(omega(n)/30000*D_new)*alpha_damped;
b2([measDOFs(1) measDOFs(3) measDOFs(4)],[measDOFs(1) measDOFs(3) measDOFs(4)]);


b(:,:) = b1 - b2([measDOFs(1) measDOFs(3) measDOFs(4)],[measDOFs(1) measDOFs(3) measDOFs(4)]);

S_full = alpha*(omega(n)/30000*eye(N,N))*alpha_damped;
S(:,:) = S_full([measDOFs(1) measDOFs(3) measDOFs(4)],[measDOFs(1) measDOFs(3) measDOFs(4)]);
x = ((b(1,1)))./((S(1,1)))
D_new = (mean(.25*x))*eye(N,N)+D_new;
f(n,m) = abs(b(1,1));
end
error(n) = abs(alpha_damped_real(1,1)) - abs(alpha([measDOFs(1)],[measDOFs(1)]));
D_new_param(n) = D_new(end,end)
end

error(:) = abs(error(:)/max(error(:)));