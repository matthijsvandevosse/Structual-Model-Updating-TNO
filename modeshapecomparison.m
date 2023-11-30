%%
clear all 
% load('G_ref_cal.mat')
% load('G_ref_wgn.mat')

% load('G_ref_wgn_0.01.mat')

% load('G_ref_tf_estimate.mat')
% load('G_ref_tf_estimate_new_cal_2.mat')
% load('G_ref_tf_estimate_new_cal_paul.mat')
load('/Users/matthijsvandevosse/FRF_Balk_Processing/G_ref_high_pass_0_04.mat')
load('G_ref.mat')
figure
plot(imag(squeeze(G_ref.ResponseData(:,1,:)))')
Gfrf_mag = abs(G_ref.ResponseData);
Gfrf_pha = rad2deg(angle(G_ref.ResponseData));
ffrf = G_ref.Frequency;

% for i = 1:5
%     offset(i) =  mean(Gfrf_mag(i,:,1));
% end
% offset = offset'./offset(1)
% for i = 1:5
% G_ref.ResponseData(i,:,:) = G_ref.ResponseData(i,:,:)./(offset(i));
% end

Gfrf_mag = abs(G_ref.ResponseData);
Gfrf_pha = rad2deg(angle(G_ref.ResponseData));
%%
n=0;
index_mode = [12 22 132 208 373 690 734 1164 1752 2572];
index_mode = [15 27 165 259 467 536 867 921 1446];
%'G_ref_cal.mat'
index_mode = [7 13 82 130 233 434 724];
% %G_ref_wgn.mat
% index_mode = [8 15 82 129 233 434 724];
% % G_ref_tf_estimate.mat
% index_mode = [17 29 165 260 468 869 923];
% % % G_ref_tf_estimate_new_cal.mat
% % index_mode = [8 13 67 105 188 350 576];
% % G_ref_tf_estimate_new_cal_2.mat
% index_mode = [25 45 266 408 743 831 1405 1483];
% % G_ref_tf_estimate_new_cal_paul.mat
index_mode = [12 22 133 204 371 702 740];
% high pass
index_mode = [7 13 83 124 228 257 436 459];

for index = index_mode
    n = n+1;
    for i = 1:3
        column_2norm(i) = norm(G_ref.ResponseData(:,i,index));
    end
    for i = 1:5
        row_2norm(i) = norm(G_ref.ResponseData(i,:,index));
    end

    R(:,:,n) = imag(squeeze(G_ref.ResponseData(:,:,index)));
    R(:,:,n) = -R(:,:,n);
    R(:,:,n) = R(:,:,n)./ sqrt([1 1.2 1.1]);
    Modeshape_sensor_old(:,n) = R(:,:,n)/column_2norm;
    Modeshape_act_old(:,n) = R(:,:,n)'/row_2norm;
    
    [U,S,V] = svd(R(:,:,n));

    Modeshape_sensor(:,n) = U(:,1)*sqrt(abs(S(1,1)));
    Modeshape_act(:,n) = V(:,1)*sqrt(abs(S(1,1)));
    
    Modeshape_sensor(:,n)*Modeshape_act(:,n)'
    R(:,:,n)

end

%% Test bode diagram with modal form
load("Modeshapes_V1.mat")
load("eigenFrequencies.mat")

% load("L_test.mat")
%  L = L*0.01;
ModalFitData.eigenFrequencies = eigenFrequencies; 


nModes = length(index_mode);
Dm = ones(1, nModes);
Dm_i = ones(1, nModes); %old

% Dm_i = ones(1, nModes);
Dm_i(1:2) = 1.5;

lambda = ModalFitData.eigenFrequencies([1 2 3 4 5]).^2;
lambda_i = G_ref.Frequency(index_mode).^2;

L(:,6:nModes) = zeros(5, nModes-5);
lambda(6:nModes) = zeros(nModes-5,1);

clear G
for i = 1:nModes

    Psi_squared = L(:,i)*(L(:,i))';
    Psi_squared(:,[1]) = 0.02*Psi_squared(:,[1])
    Psi_squared(:,[3]) = 0.02*1.2*Psi_squared(:,[3])
    Psi_squared(:,[5]) = 0.02*1.1*Psi_squared(:,[5])
    
    
    Psi_squared_i = Modeshape_sensor(:,i)*Modeshape_act(:,i)'*sqrt(lambda_i(i))*Dm_i(i);

    % Psi_squared_i = Modeshape_sensor(:,i)*Modeshape_sensor([1 3 5],i)'*sqrt(lambda_i(i))*Dm_i(i);

    % Psi_squared_i = Modeshape_sensor(:,i)*Modeshape_sensor([1 3 5],i)'*sqrt(lambda_i(i))*Dm_i(i);
    % Psi_squared_i = R(:,:,i)*sqrt(lambda_i(i))*Dm_i(i);


    % if i == 2
    %    Psi_squared_i = Psi_squared(:,[1 3 5]);
    % end

    if i > 1
        G = G + Psi_squared(:,[1 3 5]) * tf(1, [1 Dm(i) lambda(i)]);
        G_i = G_i +  Psi_squared_i * tf(1, [1 Dm_i(i) lambda_i(i)]);
    else
        G = Psi_squared(:,[1 3 5])  * tf(1, [1 Dm(i) lambda(i)]);
        G_i =  Psi_squared_i * tf(1, [1 Dm_i(i) lambda_i(i)]);
    end
end

w = ffrf;
[G_mag,G_pha, w] = bode(G, w);
[Gi_mag,Gi_pha, w] = bode(G_i, w);

for i = 1:1:5
    for j = 1:3

        if G_pha(i,j) >520
            G_pha(i,j,:) = G_pha(i,j,:) -540;
        end
        if G_pha(i,j) >340
            G_pha(i,j,:) = G_pha(i,j,:) -360;
        end
        if G_pha(i,j) >160
            G_pha(i,j,:) = G_pha(i,j,:) -180;
        end
        if Gi_pha(i,j) >520
            Gi_pha(i,j,:) = Gi_pha(i,j,:) -540;
        end
        if Gi_pha(i,j) >340
            Gi_pha(i,j,:) = Gi_pha(i,j,:) -360;
        end
        if Gi_pha(i,j) >160
            Gi_pha(i,j,:) = Gi_pha(i,j,:) -180;
        end
    end
end


Out = [1 1 2 2 3 3];
In = [1 1 2 2 3 3];
figure(26);
tiledlayout(5,3,'TileSpacing','Compact','Padding','Compact');
for i = 1:1:5
    for j = 1:3
        nexttile
        hold on

        set(gca, 'XScale', 'log')
        % grid on

        set(gca,'ytick',[-50,0,50])


        plot(ffrf,squeeze(mag2db(Gfrf_mag(i,j,:))), 'b',   'LineWidth', 1.5)

        plot(ffrf,mag2db(squeeze(Gi_mag(i,j,:))),'k', 'LineWidth',2 )
         plot(ffrf,mag2db(squeeze(G_mag(i,j,:))) , 'LineWidth',2 )


        % xline(data_opt.G_ref.Frequency(freq)*2*pi)

        if i ==1 & j == 1
            lgd = legend(["FRF",  "New","Model Paul"],Location="southwest");
        end

        lgd.FontSize = 18;
        % legend(["FRF", "Updated Zeros and Poles",],Location="southwest")
        if i==3 && j==2;
            xlabel('f [Hz]')
        elseif i~=3
            xticklabels({})
        end
        if j==1 && i==3
            ylabel('Mag [dB]')
        elseif j~=1
            yticklabels({})
        end
        xlim([1,2000])
        ylim([-75,80])
        set(gca, 'Fontsize', 13)
        set(gca, 'LineWidth', 1.25)
        box on

    end
end

%%
f = figure(27);
tiledlayout(5,3,'TileSpacing','Compact','Padding','Compact');
for i = 1:1:5
    for j = 1:1:3
        nexttile
        hold on

        set(gca, 'XScale', 'log')
        set(gca,'ytick',[-360,-180,0])


        plot(ffrf/2/pi,squeeze((Gfrf_pha(i,j,:)))+5*40/100*ffrf/2/pi/2/pi, 'b',   'LineWidth', 1.5)
        plot(w,squeeze((Gi_pha(i,j,:))),   'LineWidth', 1.5)

        if i ==1 & j == 1
            lgd = legend(["FRF", "Init", "Updated",],Location="southwest");
        end
        lgd.FontSize = 18;  % xline(data_opt.G_ref.Frequency(freq)*2*pi)

        if i==3 && j==2;
            xlabel('f [Hz]')
        elseif i~=3
            %             set(gca,'xtick',[])
            xticklabels({})
        end
        if j==1 && i==3
            ylabel('Phase [^\circ]')
        elseif j~=1
            %             set(ax2,'xticklab',get(ax1,'xticklab'))
            yticklabels({})
        end
        xlim([1,2000])
        ylim([-380,0])
        set(gca, 'Fontsize', 13)
        set(gca, 'LineWidth', 1.25)
        box on
        % fig = get(groot,'CurrentFigure');
        % set(fig, 'PaperUnits', 'centimeters');
        % set(fig, 'Units', 'centimeters');
        % x_width=8.85;          %x_width of the figure, set at template column width
        % y_width=0.75*x_width;   %y_width of the figure, free to choose
        % set(fig, 'PaperPosition', [0 0 x_width y_width]);
        % set(fig, 'PaperSize', [x_width y_width]);
        % set(fig, 'InnerPosition', [0 0 x_width y_width]);
        % set(gca,'FontSize',9)

    end
end

%%
for i = 1:length(index_mode)
    if i < 3
    [m,index] = max(abs(L([1 2 3 5],i)));
    L_norm(:,i) = L(:,i) / L(index,i);
    else
         [m,index] = max(abs(Modeshape_sensor(:,i)));  
    end
    Modeshape_sensor_norm(:,i) = Modeshape_sensor(:,i)./(Modeshape_sensor(index,i));
    [m,index] = max(abs(Modeshape_act(:,i)));
    Modeshape_act_norm(:,i) = Modeshape_act(:,i)./(Modeshape_act(index,i));
end

%%
load('sensor_0.02.mat')

figure();
tiledlayout(4,3,'TileSpacing','Compact','Padding','Compact');
for i = 1:length(index_mode)
    legend(['equilibrium = 0V'],[ 'equilibrium = -0.02V'], "Paul")
        nexttile
    hold on
        plot(Modeshape_sensor_norm(:,i))
    plot(Modeshape_sensor_norm2(:,i))
    if i <6
        plot(L_norm(:,i))
    end
end
%%
figure(100);
tiledlayout(4,3,'TileSpacing','Compact','Padding','Compact');
for i = 1:length(index_mode)
    legend('sensor', 'actuator')
        nexttile
    hold on
    plot([1 3 5], Modeshape_act_norm(:,i))
    
plot(Modeshape_sensor_norm(:,i))

end