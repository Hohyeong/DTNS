% Profile Hulling Activation Flag
IS_HULLED = true; % Hulling is deactivated in order not to interfere
% MATLAB Script Debugging Flag
ON_DEBUG = true;

if ON_DEBUG
    wtp = ones(41, 1, 'single') * 300;
    for ii = 2:41
        if abs(randn) < 0.2 
            wtp(ii) = wtp(ii-1);
        else
            wtp(ii) = single( wtp(ii-1) + ceil(randn * 10) );
        end
    end
    swh = 30.48;
    params = [1.5, 120, 4.0, 4.0, 75];
    % (1) Reaction Time (sec)
    % (2) Max Pull Up G (g)
    % (3) Roll Rate (deg/s)
    % (4) G Onset Rate (g/s)
    % (5) Max Flight Path Angle (deg)
    VerticalUncertainty = 20;
end

wtpLength = size(wtp, 1);
wtpTime = 0:0.5:20;

% Initialize
hulled = zeros(wtpLength, 1, 'single');
hullIdx = zeros(wtpLength, 1, 'uint8');
hulltemp = wtp;

if ON_DEBUG
    a = 4.0 * 9.8;
    v = 50.0;
    theta = 75 * pi/180;

    testTime = -10:0.5:10;
    y = zeros(size(testTime));

    for ii = 1:length(testTime)
        t = testTime(ii);

        y(ii) = GetProbeHeight(a, theta, v, t);

    end
%     figure;
%     plot(time, y, '.-r');
end


% If Hulling Activated
if IS_HULLED
    
    % Sort profile in descending order
    [~, hullIdx] = sort(wtp, 'descend');

    for kk = 2:length(hullIdx)
        changed = false;
        for ii = 2:length(hullIdx)
            % 제일 높은 1번째 인덱스는 건들 필요 없음

            % I번째 검사지점
            curI = hullIdx(ii);
            dh = single(0.0);

            % I번째보다 높은 (I-1)부터 1까지 비교하였을 때
            for jj = (ii-1):-1:1
                curJ = hullIdx(jj);
                dt = (curJ - curI) * 0.5;

                if hulltemp(curJ) > ( hulltemp(curI) + GetProbeHeight(a, theta, v, dt))

                    dhTemp = ( hulltemp(curJ) - ( hulltemp(curI) + GetProbeHeight(a, theta, v, dt)) );
                    if dhTemp > dh
                        dh =dhTemp;
                    else
                    end
                else
                end

            end
            if dh > 0
                hulltemp(curI) = hulltemp(curI) + dh;
                changed = true;
            else

            end

        end

        % Terimination Condition
        if changed
            [~, hullIdx] = sort(hulltemp, 'descend');
            continue
        else
            break
        end

    end

end


hulled = hulltemp + swh + 2 * VerticalUncertainty;


if ON_DEBUG
    figure('WindowStyle', 'docked')
    hold on;
    plot(wtpTime, wtp, '.-b', 'DisplayName', 'Worst-case Terrain Profile');
    plot(wtpTime, hulltemp, '-r', 'DisplayName', 'Hulled Profile');

    time = -1.5:0.5:1.5;

    for ii = 1:length(hulltemp)
        y = zeros(size(time));
        t = time + (ii-1)*0.5;
        for jj = 1:length(time)
            y(jj) = hulltemp(ii) + GetProbeHeight(a, theta, v, time(jj));
        end
        if ii == 1
            plot(t, y, ':', 'Color', [0.4 0.4 0.4], 'DisplayName', 'Probe profile');
        else
            plot(t, y, ':', 'Color', [0.4 0.4 0.4], 'HandleVisibility', 'off');
        end

    end
    % plot(wtpTime, hulled, '-m', 'DisplayName', 'Hulled + SWH + Uncertainty Profile');
    legend('Location', 'best');
    hold off;
    grid minor
    xlabel('Time (sec)');
    ylabel('Height (m)');
end

%% subfunctions
function y = GetProbeHeight(a, theta, v, t)

R = v^2 / a;
w = a / v;
T = theta / w;

if (abs(t) > T)
    y = single( R - R*cos(theta) + v*sin(theta)*( abs(t) -abs(T) ) );
else
    y = single( R - R*cos(w*t) );
end

end
