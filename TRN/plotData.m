clc; close all

mode       = 1;    % 0: No flight path, 1: Yes flight path

t          = Data_TRNpos.time;
Lat        = Data_TRNpos.signals(1).values;
Lon        = Data_TRNpos.signals(2).values;
Alt        = Data_Alt.signals.values;
TRN_output = Data_TRNoutput.signals.values;

% Pos(Sensor+TRN)
Pos_INS    = [Lat(:,1), Lon(:,1)];
Pos_GPS    = [Lat(:,2), Lon(:,2)];
Pos_TRN    = [Lat(:,3), Lon(:,3)];
Pos_EGI    = [Lat(:,4), Lon(:,4)];

% Alt(Sensor)
Alt_RALT   = Alt(:,1);
Alt_INS    = Alt(:,2);
Alt_GPS    = Alt(:,3);
Alt_EGI    = Alt(:,4);

% TRN output
TRNEstpos  = [TRN_output(:,1) TRN_output(:,2)];
TRNEstAlt  = TRN_output(:,3);
TRNEstAGL  = TRN_output(:,4);
TRNHorErr  = TRN_output(:,5);   % TRN HorUncertainty
TRNVerErr  = TRN_output(:,6);   % TRN VerUncertainty

%% plot
figure(1)
subplot(2,1,1)
plot(t, Pos_INS(:,1));
hold on
plot(t, Pos_GPS(:,1));
plot(t, Pos_TRN(:,1));
plot(t, Pos_EGI(:,1));
hold off
legend('INS','GPS','TRN','EGI')
grid on
title('Latitude[deg]')
xlabel('time(sec)'); ylabel('[deg]')
subplot(2,1,2)
plot(t, Pos_INS(:,2));
hold on
plot(t, Pos_GPS(:,2));
plot(t, Pos_TRN(:,2));
plot(t, Pos_EGI(:,2));
hold off
legend('INS','GPS','TRN','EGI')
grid on
title('Longitude[deg]')
xlabel('time(sec)'); ylabel('[deg]')

if (mode == 1)
    figure
    plot(Pos_INS(:,2),Pos_INS(:,1));
    hold on
    plot(Pos_GPS(:,2), Pos_GPS(:,1));
    plot(Pos_TRN(:,2), Pos_TRN(:,1));
    plot(Pos_EGI(:,2), Pos_EGI(:,1));
    hold off
    legend('INS','GPS','TRN','EGI')
    grid on
    title('Flight path')
    xlabel('Longitude(deg)'); ylabel('Latitude(deg)')
end

figure
plot(t, Alt_RALT);
hold on
plot(t, Alt_INS);
plot(t, Alt_GPS);
plot(t, TRNEstAlt);
plot(t, Alt_EGI);
hold off
legend('RALT','INS','GPS','TRN','EGI')
grid on
title('Altitude[m]')
xlabel('time(sec)'); ylabel('[m]')

figure
subplot(2,1,1)
plot(t, TRNEstpos(:,1))
grid on
title('TRN EstLatitude[deg]')
xlabel('time(sec)'); ylabel('[deg]')
subplot(2,1,2)
plot(t, TRNEstpos(:,2))
grid on
title('TRN EstLongitude[deg]')
xlabel('time(sec)'); ylabel('[deg]')

figure
subplot(2,1,1)
plot(t, TRNEstAlt)
grid on
title('TRN EstAltitude[m]')
xlabel('time(sec)'); ylabel('[m]')
subplot(2,1,2)
plot(t, TRNEstAGL)
grid on
title('TRN EstHeightAGL[m]')
xlabel('time(sec)'); ylabel('[m]')

figure
subplot(2,1,1)
plot(t, TRNHorErr)
grid on
title('TRN HorUncertainty[m]')
xlabel('time(sec)'); ylabel('[m]')
subplot(2,1,2)
plot(t, TRNVerErr)
grid on
title('TRN VerUncertainty[m]')
xlabel('time(sec)'); ylabel('[m]')

