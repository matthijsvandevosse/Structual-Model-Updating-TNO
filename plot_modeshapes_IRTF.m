

clf(figure(6))
clf(figure(5))
figure(6)



displaymode = 2;

x = expModes.Actuator_pos(:,1);
y = expModes.Actuator_pos(:,2);
T = delaunay(x,y);
T = unique(T,'rows');
trisurf(T, x, y, -psi_m_solved(:, displaymode))
hold on
xlabel("X")
ylabel("Y")
title(['Simulated mode at ' num2str(freq_solved(displaymode)) 'hz' ])
figure(5)
trisurf(T, x, y, -expModes.psiExp(:, displaymode))
hold on
xlabel("X")
ylabel("Y")
title(['Measured mode at ' num2str(freq_m(displaymode)) 'hz' ])
drawnow


%%
figure(7)
clf(figure(7))

displaymode = 18;
x = expModes.Actuator_pos(:,1);
y = expModes.Actuator_pos(:,2);
T = delaunay(x,y);
T = unique(T,'rows');
trisurf(T, x, y, -psi_solved(expModes.measDOFs(:,1), displaymode) + psi_solved(expModes.measDOFs(:,2), displaymode))
hold on
xlabel("X")
ylabel("Y")
title(['Simulated mode at ' num2str(sqrt(lambda_solved(displaymode))/2/pi) 'hz' ])

figure(8)
clf(figure(8))


x = expModes.Actuator_pos(:,1);
y = expModes.Actuator_pos(:,2);
T = delaunay(x,y);
T = unique(T,'rows');
trisurf(T, x, y, -psi_o(expModes.measDOFs(:,1), displaymode) + psi_o(expModes.measDOFs(:,2), displaymode))
hold on
xlabel("X")
ylabel("Y")
title(['Simulated mode at ' num2str(sqrt(lambda_solved(displaymode))/2/pi) 'hz' ])
