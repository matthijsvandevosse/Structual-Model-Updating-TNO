load('BeamModal_O12345I123.mat')

frf = permute(G_ref.ResponseData(:,:,5:300),[3 1 2]);

f = G_ref.Frequency(5:300);
frf_obj = frd(frf, f);
w = logspace(0,4,900)/2/pi;
w2 = f;
fs = 1/G_ref.Ts ;

%%
G = idfrd(permute(frf,[2 3 1]),f,0,'FrequencyUnit','Hz');

np = 42;
np = 16;
nz = 12;

sys = tfest(G,np,nz)

figure(2)
bode(G, sys)

xlim([1 500])


%%

Modes = 15
[fn,dr,ms,ofrf] = modalfit(sys, f, Modes);
fn
ofrf_p = permute(ofrf,([2 3 1]));

hm2_frf = frd(ofrf_p, f);

[hm2_mag,hm2_pha   ,wfrf]   = bode(hm2_frf, f);


%%
[frf]  = freqresp(sys,f, 'Hz');
hm2_frf = frd(frf, f);
frf = permute(frf,[3 1 2]);

[hm2_mag,hm2_pha]   = bode(hm2_frf, f);

%%
eigenFrequencies = [ModalFitData.eigenFrequencies([1 2 3])]
[fn,dr,ms,ofrf] = modalfit(frf, f, fs, 12, 'FitMethod', 'lsrf', 'FreqRange', [1 50], 'PhysFreq', eigenFrequencies);
fn

ofrf_p = permute(ofrf,([2 3 1]));

hm2_frf = frd(ofrf_p, f);
[hm2_mag,hm2_pha   ,wfrf]   = bode(hm2_frf, f);
%%
hm = freqresp(Gmod7, w);
hm_frf = frd(hm, w);
[hm_mag,hm_pha, w]   = bode(hm_frf, w);
w= w/2/pi;

%%
[Gfrf_mag,Gfrf_pha ,ffrf]   = bode(G_ref);
ffrf = ffrf/2/pi;
%%

Out = [1 1 2 2 3 3];
In = [1 1 2 2 3 3];
figure(23);
tiledlayout(5,3,'TileSpacing','Compact','Padding','Compact');
for i = 1:1:5
    for j = 1:3
        nexttile
        hold on

        set(gca, 'XScale', 'log')
        % grid on

        set(gca,'ytick',[-50,0,50])


        plot(ffrf,squeeze(mag2db(Gfrf_mag(i,j,:))), 'b',   'LineWidth', 1.5)
        % if i == 1 || i == 3 || i == 5
        % plot(w,mag2db(squeeze(hm_mag((In(i)),j,:))),'r', 'LineWidth',2 )
        % end
          plot(wfrf,mag2db(squeeze(hm2_mag(i,j,:))),'k', 'LineWidth',2 )
        
        xline(eigenFrequencies     )

        % xline(data_opt.G_ref.Frequency(freq)*2*pi)

        if i ==1 & j == 1
            lgd = legend(["FRF", "Model Paul", "New",],Location="southwest");
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

figure(24);
tiledlayout(5,3,'TileSpacing','Compact','Padding','Compact');
for i = 1:1:5
    for j = 1:3
        nexttile
        hold on

        set(gca, 'XScale', 'log')
        % grid on

        set(gca,'ytick',[-360,0,90])


        plot(ffrf,squeeze((Gfrf_pha(i,j,:))), 'b',   'LineWidth', 1.5)
        % if i == 1 || i == 3 || i == 5
        % plot(w,(squeeze(hm_pha((In(i)),j,:))),'r', 'LineWidth',2 )
        % end
          plot(wfrf,(squeeze(hm2_pha(i,j,:))),'k', 'LineWidth',2 )
        
        xline(eigenFrequencies     )

        % xline(data_opt.G_ref.Frequency(freq)*2*pi)

        if i ==1 & j == 1
            lgd = legend(["FRF", "Model Paul", "New",],Location="southwest");
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
        ylim([-360,90])
        set(gca, 'Fontsize', 13)
        set(gca, 'LineWidth', 1.25)
        box on

    end
end

%%

figure
plot(imag(squeeze(G_ref.ResponseData(:,1,:)))')

%%
n=0;
index_mode = [12 21 131 208 373 690 734 1164 1752 2572];
sign(1,:) = [1; 1; 1];
sign(2,:) = [1; -1; -1];
sign(3,:) = [1; -1; 1];
sign(4,:) = [1; -1; 1];
sign(5,:) = [1; -1; 1];
sign(6,:) = [1; -1; 1];
sign(7,:) = [1; -1; 1];
sign(8,:) = [1; -1; -1];
sign(9,:) = [1; 1; 1];
sign(10,:) = [1; 1; -1];
for index = index_mode
    n = n+1;
    for i = 1:3
        column_2norm(i) = norm(G_ref.ResponseData(:,i,index));
    end
    for i = 1:5
        row_2norm(i) = norm(G_ref.ResponseData(i,:,index));
    end

    Modeshape = imag(squeeze(G_ref.ResponseData(:,:,index)));
    Modeshape = Modeshape.*sign(n,:)
    Modeshape_sensor(:,n) = Modeshape/column_2norm
    Modeshape_act(:,n) = Modeshape'/row_2norm;
end

%%
ms_sensor = imag(ms);
%%
load("Modeshapes_V1.mat")
sign2 = [-1 -1 1 1 -1];
for i = 1:5
    [m,index] = max(abs(L([1 2 3 5],i)));
    L(:,i) = sign2(i)*L(:,i) / m;
end

for i = 1:length(index_mode)
    if i < 6
        [~, index] = max(abs(ms_sensor(:,i)));

        ms_sensor(:,i) = ms_sensor(:,i)./ms_sensor(index,i);
    end
    
    [~, index] = max(abs(Modeshape_sensor(:,i)));
    Modeshape_sensor(:,i) = Modeshape_sensor(:,i)./Modeshape_sensor(index,i);
   
end

figure(99);
tiledlayout(4,3,'TileSpacing','Compact','Padding','Compact');
for i = 1:length(index_mode)
    legend('imag','fitted', 'paul')
        nexttile
    hold on
    plot(Modeshape_sensor(:,i))
    if i < 6
    plot((ms_sensor(:,i)))
    end
    if i <6
    plot(L(:,i))
    end


end


