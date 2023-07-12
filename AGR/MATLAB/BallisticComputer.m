

x0 = 1000;
z0 = 1000;
V0 = 250;

figure;
for ii = -30:5:30
gamma = ii;
T = findFinalTime(z0, V0, gamma);
t = 0:1:ceil(T);

[x, z] = bombKinematics(x0, z0, V0, gamma, t);

hold on;
if ii == -30
    plot(x, z, '--ro', 'MarkerSize', 3);
elseif ii == 0
    plot(x, z, '--ro', 'MarkerSize', 3);
elseif ii == 30
    plot(x, z, '--ro', 'MarkerSize', 3);
else
    plot(x, z, '--ko', 'MarkerSize', 3);
end
end
grid minor;
xlabel('Horizontal Range (m)');
ylabel('Vertical Range (m)');
axis equal
axis([0 10000, 0, 2000]);


function [x, z] = bombKinematics(x0, z0, V0, gamma, t)

g = 9.8;
x = x0 + V0 * cosd(gamma) .* t;
z = z0 + V0 * sind(gamma) .* t - 0.5 * g * t.^2;
end

function T = findFinalTime(z0, V0, gamma)

g = 9.8;

T = (V0*sind(gamma) + sqrt(V0^2*(sind(gamma))^2 - 2*g*(0-z0) ) ) / g;

end