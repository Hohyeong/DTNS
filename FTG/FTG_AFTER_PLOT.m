
idxTrueState    = 1;
idxPureINS      = 2;

labels = {'Latitude (deg)', 'Longitude (deg)', 'Altitude (m)', ...
    'Roll (deg)', 'Pitch (deg)', 'Heading (deg)', ...
    'Velocity X (m/s)', 'Velocity Y (m/s)', 'Velocity Z (m/s)', ...
    'Acceleration X (m/s^2)', 'Acceleration Y (m/s^2)', 'Acceleration Z (m/s^2)'};

figure;
for ff = 1:12
    subplot(4, 3, ff);
    plot(out.logsout{idxPureINS}.Values.Data(:,ff) - out.logsout{idxTrueState}.Values.Data(:,ff), '-b');
    xlabel('Time (sec)');
    ylabel(labels{ff});
    grid minor;
end
sgtitle('Pure INS Errors (1 Hour)');


figure;
for ff = 1:12
    subplot(4, 3, ff);
    hold on;
    plot(out.logsout{idxPureINS}.Values.Data(:,ff), '-b');
    plot(out.logsout{idxTrueState}.Values.Data(:,ff), '-r');
    hold off;
    xlabel('Time (sec)');
    ylabel(labels{ff});
    grid minor;
end
sgtitle('Pure INS vs True State (1 Hour)');



figure;
for ff = 1:12
    subplot(4, 3, ff);
    hold on;
    plot(out.logsout{idxTrueState}.Values.Data(:,ff), '-r');
    hold off;
    xlabel('Time (sec)');
    ylabel(labels{ff});
    grid minor;
end
sgtitle('True State');
