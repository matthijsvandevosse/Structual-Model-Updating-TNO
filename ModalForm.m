load("Solved.mat")
%%
structModel.K = Ksolved;

[Psi_solvedxy, lambdasolved] = eig(full(structModel.M)\structModel.K, "nobalance");

[lambdasolved,dummyInd] = sort(diag(lambdasolved), 'ascend');

Psi_solvedxy = Psi_solvedxy(:,dummyInd);

structModel.D0 = 5.6e-2*eye(N);



% Dm = inv(Psi_solvedxy)*(full(structModel.M0)\structModel.D0)*(Psi_solvedxy);

for i = 1:length(Psi_solvedxy)
    m = (Psi_solvedxy(:,i))'*full(structModel.M0)*(Psi_solvedxy(:,i));

    Psi_solvedxy(:,i) = Psi_solvedxy(:,i)./sqrt(m);
    Dm(i,i) = Psi_solvedxy(:,i)'*structModel.D0*Psi_solvedxy(:,i)

    % Km = inv(Psi_solvedxy)*(full(structModel.M0)\structModel.K)*(Psi_solvedxy);
    % Km = transpose(Psi_solvedxy(:,1))*(full(structModel.M0)\structModel.K)*(Psi_solvedxy(:,1))
end
%%
%%
clear G
nModes =5;

% Psi_solvedxy =  Psi_solvedxy * L(1,1)/Psi_solvedxy([measDOFs(1)],1);


for i = 1:nModes
    % R = inv(full(structModel.M0))*Psi_solvedxy(:,i);
    % Psi_squared = inv(full(structModel.M0))*Psi_solvedxy(:,i)*Psi_solvedxy(:,i)';
    Psi_squared = 30000*Psi_solvedxy(:,i)*Psi_solvedxy(:,i)';

    if i > 1
        % G = G+ 0.0012*Psi_squared([measDOFs(1) measDOFs(3) measDOFs(4)],[measDOFs(1) measDOFs(3) measDOFs(4)])* tf(1, [1 Dm(i,i) lambdasolved(i)]);
        G = G+ Psi_squared([measDOFs(1) measDOFs(3) measDOFs(4)],actDOFs)* tf(1, [1 Dm(i,i) lambdasolved(i)]);

    else
        % G =Psi_solvedxy([measDOFs(1) measDOFs(3) measDOFs(4)],i) * Psi_solvedxy([measDOFs(1) measDOFs(3) measDOFs(4)],i)'* tf(1, [1 Dm(i,i) lambdasolved(i)]);
        % G = 0.0012*Psi_squared([measDOFs(1) measDOFs(3) measDOFs(4)],[measDOFs(1) measDOFs(3) measDOFs(4)])* tf(1, [1 Dm(i,i) lambdasolved(i)]);
        G = Psi_squared([measDOFs(1) measDOFs(3) measDOFs(4)],actDOFs)* tf(1, [1 Dm(i,i) lambdasolved(i)]);

    end
end

%%%
[G_mag,G_pha, ffrf] = bode(G, w);

for i = 1:3
    for j = 1:3
        if G_pha(i,j,1) > 600
            G_pha(i,j,:) = G_pha(i,j,:) - 720;
        elseif G_pha(i,j,1) > 200
            G_pha(i,j,:) = G_pha(i,j,:) - 360;
        end
    end
end
Out= [1 3 5];
Inp = [1:3];



f = figure(24);
tiledlayout(4,3,'TileSpacing','Compact','Padding','Compact');
for i = 1:3
    for j = 1:1:3
        nexttile
        hold on

        set(gca, 'XScale', 'log')
        set(gca,'ytick',[-360,-180,0])


        plot(ffrf,squeeze((G_pha( i,j,:))), 'b',   'LineWidth', 1.2)
        plot(w,(squeeze(h_pha(Out(i),Inp(j),:))), 'LineWidth',1.2)
        plot(w,(squeeze(h2_pha(Out(i),Inp(j),:))), 'LineWidth',1.2)


        xline(naturalfrequency*2*pi)

        legend(["FRF", "No damping", "Updated with optimization", "Updated with hand"], Location="southwest")
        if i==3 && j==2;
            xlabel('$f$ [Hz]')
        elseif i~=3
            %             set(gca,'xtick',[])
            xticklabels({})
        end
        if j==1 && i==2
            ylabel('Phase [$^\circ$]')
        elseif j~=1
            %             set(ax2,'xticklab',get(ax1,'xticklab'))
            yticklabels({})
        end
        xlim([1,2000])
        ylim([-400,0])
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


f = figure(23);
tiledlayout(3,3,'TileSpacing','Compact','Padding','Compact');
for i = 1:3
    for j = 1:3
        nexttile
        hold on

        set(gca, 'XScale', 'log')
        % grid on

        set(gca,'ytick',[-50,0,50])


        plot(ffrf,squeeze(mag2db(G_mag( i,j,:))), 'b',   'LineWidth', 1.2)
        plot(w,mag2db(squeeze(h_mag(Out(i),Inp(j),:))), 'LineWidth',1.2)
        plot(w,mag2db(squeeze(h2_mag(Out(i),Inp(j),:))), 'LineWidth',1.2)


        legend(["Decoupled", "No damping", "damping"],Location="southwest")
        if i==3 && j==2;
            xlabel('$f$ [Hz]')
        elseif i~=3
            xticklabels({})
        end
        if j==1 && i==2
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

