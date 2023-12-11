%% x = [position; velocity]

Minv_o = inv((structModel.M0));
% Minv = inv(structModel.M);
C_0 = [zeros(5,2*N)];
C_0(1,measDOFs(1)) =      35000;
C_0(2,measDOFs(2)) =      31000;
C_0(3,measDOFs(3)) =      31000;
C_0(4,nonmeasDOFs(1)) =   31000;
C_0(5,measDOFs(4)) =      31000;

% C = [zeros(5,2*N)];
% C(1,measDOFs(1)) =      30000+6000*x_1_M(1);
% C(2,measDOFs(2)) =      30000+6000*x_1_M(1);
% C(3,measDOFs(3)) =      30000+6000*x_1_M(1);
% C(4,nonmeasDOFs(1)) =   30000+6000*x_1_M(1);
% C(5,measDOFs(4)) =      30000+6000*x_1_M(1);


%%
structModel.D0 = 8e-2.*eye(N,N);
structModel.D = 8e-2.*eye(N,N);
% structModel.D0 = 0e-2.*eye(N,N);

%%
% A = [zeros(N,N) eye(N); -inv(structModel.M)*structModel.K -inv(structModel.M)*structModel.D];

A_o = [zeros(N,N) eye(N); -Minv_o*structModel.K0 zeros(N,N)];

A_p = [zeros(N,N) eye(N); -Minv_o*Kp -Minv_o*structModel.D0];


% B = [zeros(N,3); Minv(:,actDOFs)];
B_o = [zeros(N,3); Minv_o(:,actDOFs)];
load("Modal_Systeem_V1.mat")
% sys1 = ss(A,B,C,0);
syso = ss(A_o,B_o,C_0,0);
sysp = ss(A_p,B_o,C_0,0);

%%
w = logspace(0,4,800);
% h = freqresp(sys1, w);
% h_o = freqresp(syso, w);
h_p = freqresp(sysp, w);
hm = freqresp(Gmod7, w);

%%
% A_new = [zeros(N) eye(N); -inv(structModel.M)*structModel.K -inv(structModel.M)*(structModel.D)];
%
% sys3 = ss(A_new,B,C,0);
%
% h2 = freqresp(sys3, w);

%%

% h2_frf = frd(h2, w);
% h4_frf = frd(h4, w);
% 
% h_frf = frd(h, w);
hp_frf = frd(h_p, w);
% ho_frf = frd(h_o, w);
hm_frf = frd(hm, w);
%%

% [h2_mag,h2_pha   ,wfrf]   = bode(h2_frf,w);
% [h4_mag,h4_pha   ,wfrf]   = bode(h4_frf,w);
[hp_mag,hp_pha   ,wfrf]   = bode(hp_frf,w);
% [h_mag,h_pha   ,wfrf]   = bode(h_frf,w);
% [ho_mag,ho_pha   ,wfrf]   = bode(ho_frf);
[hm_mag,hm_pha   ,wfrf]   = bode(hm_frf);