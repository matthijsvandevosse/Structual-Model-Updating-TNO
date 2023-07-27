clear all; close all;
%% IFAC World Congress Basis
%% VIDI 2023 Maart
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex'); 
set(groot, 'defaulttextInterpreter','latex');

%loading data
Files = dir(fullfile(pwd, '*.mat'))

Files = Files(1:end);
data = cell(length(Files),1);

for k = 1:length(Files)
   data{k} = load((Files(k).name)).ModalFitData;
end



plotWithOutput = {[1 2 3 4 5],...
                  [1 2 3 4],...
                  [1 2 4 5],...
                  [1 2 4]};              
              
newData = cell(length(plotWithOutput));
for k = 1:length(newData)
    for k1 = 1:length(data)
       if isequal(plotWithOutput{k},data{k1}.outputsToUse)
            newData{k} = data{k1}; 
       end
    end
end
              
data = newData;
                   
nm = 5;
colors = [0      0      0     ;...
          0      0.4470 0.7410;...
          0.8500 0.3250 0.0980;...
          0.9290 0.6940 0.1250;...
          0.4940 0.1840 0.5560;...	
          0.4660 0.6740 0.1880;...
          0.3010 0.7450 0.9330;...
          0.6350 0.0780 0.1840];

% Plotting point comparison  
h1 = figure('Renderer', 'painters', 'Position', [0 0 1000 550]);
legendData = cell(nm,length(data),3);
for k = 1:nm
    subplot(2,3,k);
    for s = 1:length(data)
        points = data{s}.points;
        points{4}(3,2) = -20;
        
        outputsToUse = data{s}.outputsToUse; 
        eigenFrequencies = data{s}.eigenFrequencies;
        
        p = sortrows(points{k},1);
        
        if isequal(outputsToUse, [1 2 3 4 5])           
            legendData{k,s,1} = plot(p(:,1)+3,p(:,2),'linewidth',4,'color',colors(s,:));
        else           
            legendData{k,s,1} = plot(p(:,1)+3,p(:,2),'linewidth',2,'color',colors(s,:));
        end       
        hold on;
        legendData{k,s,2} = plot(points{k}(1:length(outputsToUse),1)+3,points{k}(1:length(outputsToUse),2),'o','markerSize',6,'LineWidth',2,'Color','r','lineStyle','none');
        if k > 2
            legendData{k,s,3} = plot(points{k}(length(outputsToUse)+1:end,1)+3,points{k}(length(outputsToUse)+1:end,2),'o','markerSize',6,'LineWidth',2,'Color','g','lineStyle','none');
        end   
    end 
    title(['$i=' num2str(k) '$']);% ', Freq: ' num2str(eigenFrequencies(k),3) 'Hz']);
    xlim([1 5]);
    xlabel('$F,X[-]$');
    set(gca,'FontSize',15);
end
% Create legend entries that define the outputs used 
subplot(2,3,5);
usableEntries = [];
legendsNames = cell(1,length(data)+2);
for s = 1:length(data)
    usableEntries = [usableEntries legendData{5,s,1}];
    data{s}.outputsToUse
    outs = '$\bar{\mathcal{L}}^i_p,\: p\in$[';
    for g = 1:5
        if any(data{s}.outputsToUse == g)
            outs = append(outs,num2str(g));
        else
            outs = append(outs,'-');
        end
        if g < 5
            outs = append(outs,' ');
        end
    end
    outs = append(outs,']');
    % Append the the identified mode shape frequencies for each line 
    outs = append(outs, string(newline));
    fqs = round(data{s}.eigenFrequencies,2);
    for f = 1:length(fqs)
        outs = append(outs,num2str(fqs(f),'%5.2f'));
        if f < length(fqs)
            outs = append(outs,'$|$');
        end
    end   
    legendsNames{s} = outs;
end


% Add entries for the circle markers 
usableEntries = [usableEntries legendData{5,1,2}];
legendsNames{length(data)+1} = '$\bar{\mathcal{L}}^i_p$'; 
k=1;
while(isempty(legendData{5,k,3}))
    k=k+1;
end
usableEntries = [usableEntries legendData{5,k,3}];
legendsNames{length(data)+2} = '$\bar{\mathcal{R}}_i^q$'; 

h = legend(usableEntries,legendsNames);
rect = get(h,'Position');
rect(1,1) = rect(1,1)+0.31;
rect(1,2) = rect(1,2)+0.02;
set(h, 'Position', rect);

%print -depsc2 indentifiedModeShapesRealBeam.eps
%saveas(h1,'indentifiedModeShapesRealBeam','png');

% return
%% Comparing mode shapes in the input and output matrix together (uses the first entry in data)

h3 = figure('Renderer', 'painters', 'Position', [0 0 1000 550]);
for k = 1:nm
    subplot(2,3,k);
    eigenFrequencies = data{1}.eigenFrequencies;
    
    L = data{1}.G_modal.C(:,1:end/2);
    R = data{1}.G_modal.B(end/2+1:end,:);  
    
    plot(1:size(L,1),L(:,k),'linewidth',2,'color',colors(1,:));
    hold on;
    p2 = plot(data{1}.SensorLocs(:,1)+3,L(:,k),'o','markerSize',6,'LineWidth',2,'Color','r','lineStyle','none');    
    p3 = plot(data{1}.ActuatorLocs(:,1)+3,R(k,:),'o','markerSize',6,'LineWidth',2,'Color','g','lineStyle','none');   
    
    title(['$i=' num2str(k), ',\: f=' num2str(round(eigenFrequencies(k),1)) '$Hz']);% ', Freq: ' num2str(eigenFrequencies(k),3) 'Hz']);
    xlim([1 5]);
    xlabel('$F,X[-]$');
    set(gca,'FontSize',15);
    
    if k <= 5
       ylim([-90 90]);         
    end
       
    if k == 1
       legend([p2 p3],'$\bar{\mathcal{L}}^i_{1,2,3,4,5}$','$\bar{\mathcal{R}}_i^{1,3,5}$','Location','northwest'); 
    end    
end





%% Bode plots 

G_ref = data{1}.G_ref;

[Gfrf_mag,Gfrf_pha   ,wfrf]   = bode(G_ref);


load('./Files_Modal/Modeshapes_V1.mat')
load('./Files_Modal/Modal_Systeem_V1.mat')

[Gmod6_mag,Gmod6_pha]   = bode(Gmod7, wfrf);


ffrf = wfrf/2/pi;
% return

%% Plot
Inp = [1,2,3];
Out = [1,3,5];

f = figure(22);
tiledlayout(3,3,'TileSpacing','Compact','Padding','Compact');
for i = 1:3
    for j = 1:1:3
        nexttile
        hold on
%         plot(ffrf,squeeze(mag2db(Gfrf_mag( Out(i),Inp(j),:))), 'b',   'LineWidth', 1.2)
%         plot(ffrf,squeeze(mag2db(Gmod1_mag(Out(i),Inp(j),:))), 'k--', 'LineWidth', 1.2)
%         plot(ffrf,squeeze(mag2db(Gmod1d_mag(i,j))), 'k', 'LineWidth', 1.2)
%         plot(ffrf,squeeze(mag2db(Gmod1_mag(i,j,:))),'r--','LineWidth',1.2)
%         plot(ffrf,(mag_Cexp_ls) ,'r','LineWidth',1.2)
        % plot(f,mag_W2_rb,'r','LineWidth',1.2)
        set(gca, 'XScale', 'log')
        % grid on
        %            axis([2 400 -180 0])
%         xlim([1,45])
        set(gca,'ytick',[-50,0,50])

        % xlabel('$f$ [Hz]')
        if i==1
            plot(ffrf,squeeze(mag2db(Gfrf_mag( Out(i),Inp(j),:))), 'b',   'LineWidth', 1.2)
        else
            plot(ffrf,squeeze(mag2db(Gfrf_mag( Out(i),Inp(j),:))), 'Color',0.6*[1,1,1],   'LineWidth', 1.2)
        end
        
%         plot(ffrf,squeeze(mag2db(Gmod6_mag(Out(i),Inp(j),:))), 'r--', 'LineWidth', 1.2)
        plot(ffrf,squeeze(mag2db(Gmod6_mag(i,j,:))), 'r--', 'LineWidth', 1.2)

        if i==3 && j==2;
            xlabel('$f$ [Hz]')
        elseif i~=3
%             set(gca,'xtick',[])
            xticklabels({})
        end
        if j==1 && i==2
            ylabel('Mag [dB]')
        elseif j~=1
%             set(ax2,'xticklab',get(ax1,'xticklab'))
            yticklabels({})
        end
        xlim([1,45])
        ylim([-50,50])
        set(gca, 'Fontsize', 13)
        set(gca, 'LineWidth', 1.25)
        box on
        fig = get(groot,'CurrentFigure');
        set(fig, 'PaperUnits', 'centimeters');
        set(fig, 'Units', 'centimeters');
        x_width=8.85;          %x_width of the figure, set at template column width
        y_width=0.75*x_width;   %y_width of the figure, free to choose
        set(fig, 'PaperPosition', [0 0 x_width y_width]); 
        set(fig, 'PaperSize', [x_width y_width]); 
        set(fig, 'InnerPosition', [0 0 x_width y_width]);
        set(gca,'FontSize',9)

    end
end
set(f, 'PaperUnits', 'centimeters');
set(f, 'Units', 'centimeters');
x_width=9.85;          %x_width of the figure, set at template column width
y_width=1.14*0.75*x_width;   %y_width of the figure, free to choose
set(f, 'PaperPosition', [0 0 x_width y_width]); 
set(f, 'PaperSize', [x_width y_width]); 
set(f, 'InnerPosition', [0 0 x_width y_width]);
% set(gca,'FontSize',9)
set(gcf,'renderer','painter')
% saveas(gcf, [filePathImg1 ['Fig_Bode']], 'pdf');





f = figure(23);
tiledlayout(3,3,'TileSpacing','Compact','Padding','Compact');
for i = 1:3
    for j = 1:1:3
        nexttile
        hold on
%         plot(ffrf,squeeze(mag2db(Gfrf_mag( Out(i),Inp(j),:))), 'b',   'LineWidth', 1.2)
%         plot(ffrf,squeeze(mag2db(Gmod1_mag(Out(i),Inp(j),:))), 'k--', 'LineWidth', 1.2)
%         plot(ffrf,squeeze(mag2db(Gmod1d_mag(i,j))), 'k', 'LineWidth', 1.2)
%         plot(ffrf,squeeze(mag2db(Gmod1_mag(i,j,:))),'r--','LineWidth',1.2)
%         plot(ffrf,(mag_Cexp_ls) ,'r','LineWidth',1.2)
        % plot(f,mag_W2_rb,'r','LineWidth',1.2)
        set(gca, 'XScale', 'log')
        % grid on
        %            axis([2 400 -180 0])
%         xlim([1,45])
        set(gca,'ytick',[-50,0,50])

        % xlabel('$f$ [Hz]')
        if i==1
            plot(ffrf,squeeze(mag2db(Gfrf_mag( Out(i),Inp(j),:))), 'b',   'LineWidth', 1.2)
        else
            plot(ffrf,squeeze(mag2db(Gfrf_mag( Out(i),Inp(j),:))), 'Color',0.6*[1,1,1],   'LineWidth', 1.2)
        end
        
%         plot(ffrf,squeeze(mag2db(Gmod6_mag(Out(i),Inp(j),:))), 'r--', 'LineWidth', 1.2)
        plot(ffrf,squeeze(mag2db(Gmod6_mag(i,j,:))), 'r--', 'LineWidth', 1.2)

        if i==3 && j==2;
            xlabel('$f$ [Hz]')
        elseif i~=3
%             set(gca,'xtick',[])
            xticklabels({})
        end
        if j==1 && i==2
            ylabel('Mag [dB]')
        elseif j~=1
%             set(ax2,'xticklab',get(ax1,'xticklab'))
            yticklabels({})
        end
        xlim([1,45])
        ylim([-50,50])
        set(gca, 'Fontsize', 13)
        set(gca, 'LineWidth', 1.25)
        box on
        fig = get(groot,'CurrentFigure');
        set(fig, 'PaperUnits', 'centimeters');
        set(fig, 'Units', 'centimeters');
        x_width=8.85;          %x_width of the figure, set at template column width
        y_width=0.75*x_width;   %y_width of the figure, free to choose
        set(fig, 'PaperPosition', [0 0 x_width y_width]); 
        set(fig, 'PaperSize', [x_width y_width]); 
        set(fig, 'InnerPosition', [0 0 x_width y_width]);
        set(gca,'FontSize',9)

    end
end
set(f, 'PaperUnits', 'centimeters');
set(f, 'Units', 'centimeters');
x_width=9.85;          %x_width of the figure, set at template column width
y_width=1.25*0.75*x_width;   %y_width of the figure, free to choose
set(f, 'PaperPosition', [0 0 x_width y_width]); 
set(f, 'PaperSize', [x_width y_width]); 
set(f, 'InnerPosition', [0 0 x_width y_width]);
% set(gca,'FontSize',9)
set(gcf,'renderer','painter')
% saveas(gcf, [filePathImg2 ['Fig_Bode']], 'pdf');



f = figure(24);
tiledlayout(3,3,'TileSpacing','Compact','Padding','Compact');
for i = 1:3
    for j = 1:1:3
        nexttile
        hold on
%         plot(ffrf,squeeze(mag2db(Gfrf_mag( Out(i),Inp(j),:))), 'b',   'LineWidth', 1.2)
%         plot(ffrf,squeeze(mag2db(Gmod1_mag(Out(i),Inp(j),:))), 'k--', 'LineWidth', 1.2)
%         plot(ffrf,squeeze(mag2db(Gmod1d_mag(i,j))), 'k', 'LineWidth', 1.2)
%         plot(ffrf,squeeze(mag2db(Gmod1_mag(i,j,:))),'r--','LineWidth',1.2)
%         plot(ffrf,(mag_Cexp_ls) ,'r','LineWidth',1.2)
        % plot(f,mag_W2_rb,'r','LineWidth',1.2)
        set(gca, 'XScale', 'log')
        % grid on
        %            axis([2 400 -180 0])
%         xlim([1,45])
        set(gca,'ytick',[-360,-180,0])
        if squeeze(abs(Gmod6_pha( (i),(j),1))) > 30
            Gmod6_pha( (i),(j),:) = Gmod6_pha( (i),(j),:) -360*ones(size(Gmod6_pha( (i),(j),:)));
        end
        if squeeze(abs(Gmod6_pha( (i),(j),1))) > 30
            Gmod6_pha( (i),(j),:) = Gmod6_pha( (i),(j),:) -360*ones(size(Gmod6_pha( (i),(j),:)));
        end
        % xlabel('$f$ [Hz]')
        if i==1
            plot(ffrf,squeeze((Gfrf_pha( Out(i),Inp(j),:)))+40/100*ffrf, 'b',   'LineWidth', 1.2)
        else
            plot(ffrf,squeeze((Gfrf_pha( Out(i),Inp(j),:)))+40/100*ffrf, 'Color',0.6*[1,1,1],   'LineWidth', 1.2)
        end
        
%         plot(ffrf,squeeze(mag2db(Gmod6_mag(Out(i),Inp(j),:))), 'r--', 'LineWidth', 1.2)
        plot(ffrf,squeeze((Gmod6_pha(i,j,:))), 'r--', 'LineWidth', 1.2)

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
        xlim([1,45])
        ylim([-400,0])
        set(gca, 'Fontsize', 13)
        set(gca, 'LineWidth', 1.25)
        box on
        fig = get(groot,'CurrentFigure');
        set(fig, 'PaperUnits', 'centimeters');
        set(fig, 'Units', 'centimeters');
        x_width=8.85;          %x_width of the figure, set at template column width
        y_width=0.75*x_width;   %y_width of the figure, free to choose
        set(fig, 'PaperPosition', [0 0 x_width y_width]); 
        set(fig, 'PaperSize', [x_width y_width]); 
        set(fig, 'InnerPosition', [0 0 x_width y_width]);
        set(gca,'FontSize',9)

    end
end
set(f, 'PaperUnits', 'centimeters');
set(f, 'Units', 'centimeters');
x_width=9.85;          %x_width of the figure, set at template column width
y_width=1.14*0.75*x_width;   %y_width of the figure, free to choose
set(f, 'PaperPosition', [0 0 x_width y_width]); 
set(f, 'PaperSize', [x_width y_width]); 
set(f, 'InnerPosition', [0 0 x_width y_width]);
% set(gca,'FontSize',9)
set(gcf,'renderer','painter')
% saveas(gcf, [filePathImg1 ['Fig_Bode_Phase']], 'pdf');




%%


f = figure(25);
tiledlayout(3,1,'TileSpacing','Compact','Padding','Compact');
for i = 1:3
        nexttile
        hold on
        scatter(Lx(3:2:end),L(3:2:end,i),50,0.6*[1 1 1],'LineWidth',2)
%         scatter(Rx,R(i,:),'r')
        scatter(Lx(1),L6(1,i),50,'b','LineWidth',2)
        scatter(Rx,R6(i,:),50,'r+','LineWidth',2)
%         plot(ffrf,squeeze(mag2db(Gfrf_mag (Out(i),Inp(j),:))), 'b',   'LineWidth', 1.2)
%         plot(ffrf,squeeze(mag2db(Gmod2_mag(Out(i),Inp(j),:))), 'r--', 'LineWidth', 1.2)
%         plot(ffrf,squeeze(mag2db(Gmod1_mag(i,j,:))),'r--','LineWidth',1.2)
%         plot(ffrf,(mag_Cexp_ls) ,'r','LineWidth',1.2)
        % plot(f,mag_W2_rb,'r','LineWidth',1.2)
%         set(gca, 'XScale', 'log')
        % grid on
        %            axis([2 400 -180 0])
%         xlim([0.2,70])
        % xlabel('$f$ [Hz]')
        if i==1
            ylabel('$\phi_{\mathrm{a},1}(\rho)$')
%             set(gca,'xtick',[0,0.25,0.5,0.75,1])
            set(gca,'xtick',[0,0.5,1])
            xticklabels({})
        elseif i==2
            ylabel('$\phi_{\mathrm{a},2}(\rho)$')
%             set(gca,'xtick',[0,0.25,0.5,0.75,1])
            set(gca,'xtick',[0,0.5,1])
            xticklabels({})
        elseif i==3
            ylabel('$\phi_{\mathrm{a},3}(\rho)$')
            xlabel('$\rho$')
%             set(gca,'xtick',[0,0.25,0.5,0.75,1])
            set(gca,'xtick',[0,0.5,1])
        end
        ylim([-80,90])
        
        set(gca, 'Fontsize', 9)
        set(gca, 'LineWidth', 1.25)
        box on
%         fig = get(groot,'CurrentFigure');
%         set(fig, 'PaperUnits', 'centimeters');
%         set(fig, 'Units', 'centimeters');
%         x_width=8.85;          %x_width of the figure, set at template column width
%         y_width=0.75*x_width;   %y_width of the figure, free to choose
%         set(fig, 'PaperPosition', [0 0 x_width y_width]); 
%         set(fig, 'PaperSize', [x_width y_width]); 
%         set(fig, 'InnerPosition', [0 0 x_width y_width]);
%         set(gca,'FontSize',9)
end
set(f, 'PaperUnits', 'centimeters');
set(f, 'Units', 'centimeters');
x_width=8.85;          %x_width of the figure, set at template column width
y_width=1.05*0.75*x_width;   %y_width of the figure, free to choose
set(f, 'PaperPosition', [0 0 x_width y_width]); 
set(f, 'PaperSize', [x_width y_width]); 
set(f, 'InnerPosition', [0 0 x_width y_width]);
% set(gca,'FontSize',9)
set(gcf,'renderer','painter')
% saveas(gcf, [filePathImg1 ['Fig_Modeshapes']], 'pdf');















return
bodeplot(G_ref,'b',options); % 
set(findall(gcf,'Color',[0 0 1]),'Color',colors(5,:));
xlim([0.5,70])

hold on; 
bodeplot(data{1}.G_modal([3,5],:),'b:',options);
set(findall(gcf,'Color',[0 0 1]),'Color',colors(1,:));
xlim([0.5,70])

bodeplot(data{3}.G_modal([3,5],:),'b:',options);
set(findall(gcf,'Color',[0 0 1]),'Color',colors(3,:));
xlim([0.5,70])

bodeplot(data{4}.G_modal([3,5],:),'b:',options);
set(findall(gcf,'Color',[0 0 1]),'Color',colors(4,:));

xlim([0.5,70])
set(findall(gcf,'type','line'),'linewidth',1.5);

%print -depsc2 BodeComparisonRealBeam.eps
%saveas(h2,'BodeComparisonRealBeam','png');



