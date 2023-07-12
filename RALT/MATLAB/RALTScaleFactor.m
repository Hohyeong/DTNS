
altitude = 1:12000;
noise = zeros(1, 12000);

for ii = 1:12000
    
    if altitude(ii) <= 100
        noise(ii) = 2;
    elseif (altitude(ii) > 100) && (altitude(ii) <= 5000)
        noise(ii) = altitude(ii) * (2/100);
    elseif (altitude(ii) > 5000) && (altitude(ii) <= 10000)
        noise(ii) = 100;
    elseif (altitude(ii) > 10000)
        noise(ii) = altitude(ii) * (1/100);
    else
        noise(ii) = 0.6096 * (1/3);
    end
    
    
end


figure;
plot(altitude, noise, '-k');
xlabel('Altitude [ft]'); ylabel('Variance of Error(3\sigma) [ft]');
grid on;