% save('scene16.mat', 'Uncertainty', 'TRN_error', 'uav_data', 'LLA_log')

orgin_data = load("scene10.mat");
lci_data = load("scene11.mat");

utm2lat = 110961.7060516;
utm2lon = 89476.51124;

%% LL
figure
plot(orgin_data.uav_data.time(:), orgin_data.LLA_log.signals(1).values(:, 1));
hold on;
plot(orgin_data.uav_data.time(:), orgin_data.LLA_log.signals(1).values(:, 2));
hold on;
plot(orgin_data.uav_data.time(:), orgin_data.LLA_log.signals(1).values(:, 3));
hold on;
plot(orgin_data.uav_data.time(:), lci_data.LLA_log.signals(1).values(:, 2));
hold on;
plot(orgin_data.uav_data.time(:), lci_data.LLA_log.signals(1).values(:, 3));
hold off; grid on; legend('True', 'INS', 'TRN only', 'LCI','LCI + TRN');
xlabel('time'); ylabel('latitude [deg]');
set(gca,'FontSize',18)

figure
plot(orgin_data.uav_data.time(:), orgin_data.LLA_log.signals(2).values(:, 1));
hold on;
plot(orgin_data.uav_data.time(:), orgin_data.LLA_log.signals(2).values(:, 2));
hold on;
plot(orgin_data.uav_data.time(:), orgin_data.LLA_log.signals(2).values(:, 3));
hold on;
plot(orgin_data.uav_data.time(:), lci_data.LLA_log.signals(2).values(:, 2));
hold on;
plot(orgin_data.uav_data.time(:), lci_data.LLA_log.signals(2).values(:, 3));
hold off; grid on; legend('True', 'INS', 'TRN only', 'LCI','LCI + TRN');
xlabel('time'); ylabel('longitude [deg]');
set(gca,'FontSize',18)

%% LL Error
figure
plot(orgin_data.uav_data.time(:), abs(orgin_data.TRN_error.signals(1).values(:)));
hold on;
plot(orgin_data.uav_data.time(:), abs(lci_data.TRN_error.signals(1).values(:)));
hold on;
plot(orgin_data.uav_data.time(:),  abs(orgin_data.LLA_log.signals(1).values(:, 1) - orgin_data.LLA_log.signals(1).values(:, 2))*utm2lat );
hold off; grid on; legend('TRN only', 'LCI + TRN', 'INS');
xlabel('time'); ylabel('North [m]');
set(gca,'FontSize',18)

figure
plot(orgin_data.uav_data.time(:), abs(orgin_data.TRN_error.signals(2).values(:)));
hold on;
plot(orgin_data.uav_data.time(:), abs(lci_data.TRN_error.signals(2).values(:)));
hold on;
plot(orgin_data.uav_data.time(:),  abs(orgin_data.LLA_log.signals(2).values(:, 1) - orgin_data.LLA_log.signals(2).values(:, 2))*utm2lat );
hold off; grid on; legend('TRN only', 'LCI + TRN', 'INS');
xlabel('time'); ylabel('East [m]');
set(gca,'FontSize',18)