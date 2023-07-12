clc;

hDEM = out.hDEM.signals.values;
time = out.hDEM.time;
a=area(time,hDEM);
a(1).FaceColor = [0.4660 0.6740 0.1880];
a(1).EdgeColor = 'k';
title('Terrain Height')
xlabel('time[sec]'); ylabel('[m]')
set(gca, 'color', 'none')
grid on
xlim([0,inf])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GPS Notvalid -> (Flat)50m,3m / (Smooth,Moderate)30m,5m / (Rough)30m,10m        %
%                                                                                %
% GPS Valid -> (Flat)50m,3m / (Smooth)30m,3m / (Moderate)30m,5m / (Rough)30m,10m %
%                                                                                %
% GPS Pcode -> (Flat,Smooth)20m,3m / (Moderate)20m,5m / (Rough)20m,10m           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
GPSDataNotValid = out.GPS_state.signals(1).values;
Ccode = out.GPS_state.signals(2).values;
Pcode = out.GPS_state.signals(3).values;
roughness = mean(out.TRNError.signals(3).values);
TRNHorErr = abs(median(sort(out.TRNError.signals(1).values)));
TRNVerErr = abs(median(sort(out.TRNError.signals(2).values)));

if GPSDataNotValid == 1                              % GPS Notvalid
    if roughness < 2                                 % Flat
        if TRNHorErr < 50
            fprintf("TRNHorTest[CEP] result\n standard 50m, test %fm[Pass]\n",TRNHorErr);
        else
            fprintf("TRNHorTest[CEP] result\n standard 50m, test %fm[Fail]\n",TRNHorErr);
        end
        if TRNVerErr < 3
            fprintf("TRNVerTest[LEP] result\n standard 3m, test %fm[Pass]\n",TRNVerErr);
        else
            fprintf("TRNVerTest[LEP] result\n standard 3m, test %fm[Fail]\n",TRNVerErr);
        end
    elseif (roughness >= 2) && (roughness < 10)       % Smooth/Moderate
        if TRNHorErr < 30
            fprintf("TRNHorTest[CEP] result\n standard 30m, test %fm[Pass]\n",TRNHorErr);
        else
            fprintf("TRNHorTest[CEP] result\n standard 30m, test %fm[Fail]\n",TRNHorErr);
        end
        if TRNVerErr < 5
            fprintf("TRNVerTest[LEP] result\n standard 5m, test %fm[Pass]\n",TRNVerErr);
        else
            fprintf("TRNVerTest[LEP] result\n standard 5m, test %fm[Fail]\n",TRNVerErr);
        end
    elseif (roughness >= 10) && (roughness < 20)      % Rough
        if TRNHorErr < 30
            fprintf("TRNHorTest[CEP] result\n standard 30m, test %fm[Pass]\n",TRNHorErr);
        else
            fprintf("TRNHorTest[CEP] result\n standard 30m, test %fm[Fail]\n",TRNHorErr);
        end
        if TRNVerErr < 10
            fprintf("TRNVerTest[LEP] result\n standard 10m, test %fm[Pass]\n",TRNVerErr);
        else
            fprintf("TRNVerTest[LEP] result\n standard 10m, test %fm[Fail]\n",TRNVerErr);
        end
    end
end
if GPSDataNotValid == 0                               % GPS Valid
    if roughness < 2                                  % Flat
        if TRNHorErr < 50
            fprintf("TRNHorTest[CEP] result\n standard 50m, test %fm[Pass]\n",TRNHorErr);
        else
            fprintf("TRNHorTest[CEP] result\n standard 50m, test %fm[Fail]\n",TRNHorErr);
        end
        if TRNVerErr < 3
            fprintf("TRNVerTest[LEP] result\n standard 3m, test %fm[Pass]\n",TRNVerErr);
        else
            fprintf("TRNVerTest[LEP] result\n standard 3m, test %fm[Fail]\n",TRNVerErr);
        end
    elseif (roughness >= 2) && (roughness < 5)        % Smooth
        if TRNHorErr < 30
            fprintf("TRNHorTest[CEP] result\n standard 30m, test %fm[Pass]\n",TRNHorErr);
        else
            fprintf("TRNHorTest[CEP] result\n standard 30m, test %fm[Fail]\n",TRNHorErr);
        end
        if TRNVerErr < 3
            fprintf("TRNVerTest[LEP] result\n standard 3m, test %fm[Pass]\n",TRNVerErr);
        else
            fprintf("TRNVerTest[LEP] result\n standard 3m, test %fm[Fail]\n",TRNVerErr);
        end
    elseif (roughness >= 5) && (roughness < 10)        % Moderate
        if TRNHorErr < 30
            fprintf("TRNHorTest[CEP] result\n standard 30m, test %fm[Pass]\n",TRNHorErr);
        else
            fprintf("TRNHorTest[CEP] result\n standard 30m, test %fm[Fail]\n",TRNHorErr);
        end
        if TRNVerErr < 5
            fprintf("TRNVerTest[LEP] result\n standard 5m, test %fm[Pass]\n",TRNVerErr);
        else
            fprintf("TRNVerTest[LEP] result\n standard 5m, test %fm[Fail]\n",TRNVerErr);
        end
    elseif (roughness >= 10) && (roughness < 20)       % Rough
        if TRNHorErr < 30
            fprintf("TRNHorTest[CEP] result\n standard 30m, test %fm[Pass]\n",TRNHorErr);
        else
            fprintf("TRNHorTest[CEP] result\n standard 30m, test %fm[Fail]\n",TRNHorErr);
        end
        if TRNVerErr < 10
            fprintf("TRNVerTest[LEP] result\n standard 10m, test %fm[Pass]\n",TRNVerErr);
        else
            fprintf("TRNVerTest[LEP] result\n standard 10m, test %fm[Fail]\n",TRNVerErr);
        end
    end
end
if (any(GPSDataNotValid) == 0) && (any(Pcode) == 1)   % GPS Pcode
    if roughness < 5                                  % Flat/Smooth
        if TRNHorErr < 20
            fprintf("TRNHorTest[CEP] result\n standard 20m, test %fm[Pass]\n",TRNHorErr);
        else
            fprintf("TRNHorTest[CEP] result\n standard 20m, test %fm[Fail]\n",TRNHorErr);
        end
        if TRNVerErr < 3
            fprintf("TRNVerTest[LEP] result\n standard 3m, test %fm[Pass]\n",TRNVerErr);
        else
            fprintf("TRNVerTest[LEP] result\n standard 3m, test %fm[Fail]\n",TRNVerErr);
        end
    elseif (roughness >= 5) && (roughness < 10)        % Moderate
        if TRNHorErr < 20
            fprintf("TRNHorTest[CEP] result\n standard 20m, test %fm[Pass]\n",TRNHorErr);
        else
            fprintf("TRNHorTest[CEP] result\n standard 20m, test %fm[Fail]\n",TRNHorErr);
        end
        if TRNVerErr < 5
            fprintf("TRNVerTest[LEP] result\n standard 5m, test %fm[Pass]\n",TRNVerErr);
        else
            fprintf("TRNVerTest[LEP] result\n standard 5m, test %fm[Fail]\n",TRNVerErr);
        end
    elseif (roughness >= 10) && (roughness < 20)       % Rough
        if TRNHorErr < 20
            fprintf("TRNHorTest[CEP] result\n standard 20m, test %fm[Pass]\n",TRNHorErr);
        else
            fprintf("TRNHorTest[CEP] result\n standard 20m, test %fm[Fail]\n",TRNHorErr);
        end
        if TRNVerErr < 10
            fprintf("TRNVerTest[LEP] result\n standard 10m, test %fm[Pass]\n",TRNVerErr);
        else
            fprintf("TRNVerTest[LEP] result\n standard 10m, test %fm[Fail]\n",TRNVerErr);
        end
    end
end


