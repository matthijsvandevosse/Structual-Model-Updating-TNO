load('BeamModal_O12345I123.mat')

frf = permute(G_ref.ResponseData,[3 1 2]);
f = G_ref.Frequency;
w = f;
fs = 1/G_ref.Ts ;
%%
[fn,dr,ms,ofrf] = modalfit(frf, f, fs ,15, 'FitMethod', 'lsrf', 'FreqRange', [1 200], 'PhysFreq', [eigenFrequencies; 1080/2/pi], DriveIndex=[3 3]);
fn
%%
ofrf_p = permute(ofrf,([2 3 1]));
%%
hm_frf = frd(ofrf_p, f);
[hm_mag,hm_pha   ,wfrf]   = bode(hm_frf, f);
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

        plot(w,mag2db(squeeze(hm_mag(i,j,:))), 'LineWidth',2 )

        xline(fn ...
            )

        % xline(data_opt.G_ref.Frequency(freq)*2*pi)

        if i ==1 & j == 1
            lgd = legend(["FRF", "Init", "Updated",],Location="southwest");
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

figure
plot(imag(squeeze(G_ref.ResponseData(:,1,:)))')

%%
n=0;
index_mode = [12 21 131 208 373 690 734 1164];
sign(1,:) = [1; 1; 1];
sign(2,:) = [1; -1; -1];
sign(3,:) = [1; -1; 1];
sign(4,:) = [1; -1; 1];
sign(5,:) = [1; -1; 1];
sign(6,:) = [1; -1; 1];
sign(7,:) = [1; -1; 1];
sign(8,:) = [1; 1; -1];
sign(9,:) = [1; 1; -1];

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
for i = 1:length(index_mode)
    if i < 7
        [~, index] = max(abs(ms_sensor(:,i)));

        ms_sensor(:,i) = ms_sensor(:,i)./ms_sensor(index,i);
    end

    [~, index] = max(abs(Modeshape_sensor(:,i)));
    Modeshape_sensor(:,i) = Modeshape_sensor(:,i)./Modeshape_sensor(index,i);
end

figure();
tiledlayout(3,3,'TileSpacing','Compact','Padding','Compact');
for i = 1:length(index_mode)
    nexttile
    if i < 7
    plot((ms_sensor(:,i)))
    end
    hold on
    plot(Modeshape_sensor(:,i))
end
