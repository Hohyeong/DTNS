close all
clc; % No clear(Simulink data delete)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Run 'Roughness2_path.slx' %%%%%%%%   %% Roughness_MAP_path gen
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Pos_map_path = [Roughness2_path.signals(4).values Roughness2_path.signals(3).values];
% roughness_map_path = [Roughness2_path.signals(1).values];
% hDEM_map_path = [Roughness2_path.signals(2).values];
% map_path_X=Pos_map_path(:,1);
% map_path_Y=Pos_map_path(:,2);
% map_path_Z=roughness_map_path;
% save('Roughness_map_path.mat','map_path_X','map_path_Y','map_path_Z')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Run 'Roughness2_range.slx' %%%%%%%%   %% Roughness_MAP_range gen
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Pos_map_range = [Roughness2_range.signals(4).values Roughness2_range.signals(3).values];
% roughness_map_range = [Roughness2_range.signals(1).values];
% hDEM_map_range = [Roughness2_range.signals(2).values];
% map_range_X=Pos_map_range(:,1);
% map_range_Y=Pos_map_range(:,2);
% map_range_Z=roughness_map_range;
% save('Roughness_map_range.mat','map_range_X','map_range_Y','map_range_Z')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Run 'TRN_test.slx' %%%%%%%%%   %% Roughness_UAV flight path
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Pos_uav = [out.uav_data.signals(4).values out.uav_data.signals(3).values];
roughness_uav = [out.uav_data.signals(1).values];
hDEM_uav = [out.uav_data.signals(2).values];
uav_X=Pos_uav(:,1);
uav_Y=Pos_uav(:,2);
uav_Z=roughness_uav;
save('Roughness_uav.mat','uav_X','uav_Y','uav_Z')



%% plot 2D
% 
load('Roughness_map_range.mat')
figure(1)
geoscatter(map_range_Y(1:length(map_range_Y)-2),map_range_X(1:length(map_range_X)-2),'filled')     % MAP area plot 
alpha(.05)
hold on
geoplot(uav_Y,uav_X,'r','LineWidth',2)                                     % UAV path plot 
geolimits([33,44],[126,132])
geobasemap('colorterrain')
title('MAP');


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% map range %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot 3D map_path data
clc;

load('Roughness_map_range.mat')

k = 0;
x = 89;
y = 126;

map_range_Z_new = zeros(y, x);
for i = 1 : 1 : y
    for j = 1 : 1 : x
        k = k + 1;
        map_range_Z_new(i, j) = map_range_Z(k);
    end
%     map_range_Y_new(i) = map_range_Y(i*y, 1);
end

%% pos plot map_path and uav data about roughness
load('Roughness_uav.mat')
f=figure(2);
f.Position(3:4) = [560 840];
hold on
% pcolor(map_range_X(1:x, 1), map_range_Y(1:92:(length(map_range_Y)) - 83, 1), map_range_Z_new)
pcolor(map_range_X(1:x, 1), map_range_Y(1:x: x*y), map_range_Z_new)
% pcolor(map_range_Y(1:y, 1), map_range_X(1:x, 1), map_range_Z_new)
colorbar
caxis([0 15]);
plot(uav_X,uav_Y,'r','LineWidth',2)
axis([124.2,131.6,33.2,43.7])
title('Roughness'); xlabel('Longitude [deg]'); ylabel('Latitude [deg]');


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% map path %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot 3D map_path data
% load('Roughness_map_path.mat');
% map_path_Z_new = zeros(747, 747);
% k = 0;
% 
% for i = 1 : 1 : 747
%     for j = 1 : 1 : 747
%         k = k + 1;
%         map_path_Z_new(i, j) = map_path_Z(k);
%     end
%     map_path_Y_new(i) = map_path_Y((i*747)-746,1);
% end
% 
% %% pos plot map_path and uav data about roughness
% figure(1)
% hold on
% pcolor(map_path_X(1:747, 1), map_path_Y_new(1,1:747), map_path_Z_new)
% colorbar
% plot(uav_X,uav_Y,'r','LineWidth',2)
% axis([128.325,128.575,34.825,35.2417])   % moderate
% % axis([127.8,128.8,35,35.4])
% title('Roughness'); xlabel('Longitude [deg]'); ylabel('Latitude [deg]');

%%
% figure(2)
% surf(map_X(1:19, 1), map_Y_new, map_Z_new)
% colorbar
% title('Roughness'); xlabel('lon [deg]'); ylabel('lat [deg]'); zlabel('roughness')
% 
% %% pos plot_map or uav data
% figure(3)
% subplot(2,1,1)
% plot(uav_Y)
% findpeaks(uav_Y)          % max
% title('Latitude');
% subplot(2,1,2)
% plot(uav_X)
% findpeaks(uav_X)          % max
% title('Longitude');

%% pos plot map_path and uav data about roughness
% figure(4)
% hold on
% pcolor(map_path_X(1:747, 1), map_path_Y_new(1,1:747), map_path_Z_new)
% colorbar
% plot(uav_X,uav_Y,'r','LineWidth',2)
% axis([128.325,128.575,34.825,35.2417])   % moderate
% % axis([127.8,128.8,35,35.4])
% title('Roughness'); xlabel('Longitude [deg]'); ylabel('Latitude [deg]');

mean_roughness = mean(out.uav_data.signals(1).values)
