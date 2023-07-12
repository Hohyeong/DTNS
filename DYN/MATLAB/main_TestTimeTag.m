time = 0:0.04:1000;


timeTag = uint16(mod(time/(6.4e-5), 65535));

I = find(timeTag == 0);

figure;
plot(time(1:1000), timeTag(1:1000), '.b', 'MarkerSize', 2);
grid minor
xlabel('Time (sec)'); ylabel('Time Tag (uint16)')