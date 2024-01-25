    [psi, lambda] = eigs(Kinit,  Minit, 100,  (160*2*pi)^2, 'IsSymmetricDefinite', 1);
    lambda = diag(lambda);
    sqrt(lambda)/2/pi

%%
clear R FRF

freq = logspace(0,4,10000);
Dm(1:length(lambda)) = 0;
for ii = 1:length(freq)
    w = freq(ii);
    for m = 1:100%length(lambda)
        R(:,:,m) = psi(node_mapping.z(input_nodes_top),m)*(psi(node_mapping.z(input_nodes_top),m) - psi(node_mapping.z(input_nodes_back),m)  )';
       
        if m == 1
            FRF(:,:,ii) = 10000000000* R(:,:,m)./( -w^2 + 1i*Dm(m)*w + lambda(m) );
            
        else
            FRF(:,:,ii) = FRF(:,:,ii) +10000000000* R(:,:,m)./( -w^2 + 1i*Dm(m)*w + lambda(m) );
            
        end
    end
end

%% 
figure(10)
semilogx(freq/2/pi,mag2db(abs(squeeze(FRF(5,5,:)))), 'LineWidth', 1)
xscale log
% xlim([100, 3000])
%%
FRF_file = 'G_bla_mech';




save_result  = 'off';                   % 'on'/'off' to save result
version      = 1;                      % version number of result file

FRF_file_info = whos('-file',FRF_file);
names = {FRF_file_info.name};
% FRF_name = input('name (or number in above list) of the FRF data  ');
FRF_name = 'G_bla';

if ~ischar(FRF_name)
    FRF_name = names{FRF_name};
end

G_frf = load(FRF_file,FRF_name);

[G_frf2,omega,Ts] = frdata(G_frf.(FRF_name));
omega = omega/2/pi;