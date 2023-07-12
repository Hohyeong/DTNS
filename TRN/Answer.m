%% Q2 answer

error_0 = load("error_0.mat");
error_1 = load("error_1.mat");

err_0_1(1, :) = (error_1.TRN_error.signals(1).values - error_0.TRN_error.signals(1).values);
err_0_1(2, :) = (error_1.TRN_error.signals(2).values - error_0.TRN_error.signals(2).values);

err_0_1_sort = err_0_1(1:2, 1:5000);


figure
subplot(2, 1, 1)
plot(err_0_1_sort(1, :))
grid on; ylabel('lat error [m]'); xlabel('time');
title('R0 / M,N Compare');

subplot(2, 1, 2)
plot(err_0_1_sort(2, :))
grid on; ylabel('lon error [m]'); xlabel('time');

%% 