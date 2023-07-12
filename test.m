
freq_uint8 = uint8(25);

y = zeros(100, 1);

for ii = 1:100

    y(ii) = GenerateUpdateFlag(freq_uint8);

end

figure;
plot(y, 'ro');



function updateFlag = GenerateUpdateFlag(freq_uint8)
%GenerateUpdateFlag RALT 데이터 갱신 플래그 생성
% Input
% - freq_uint8  1 = update data at every simulation step
%               25 = update data at every 25 simulation steps
%               Accepts  (8bit unsigned integer=uint8) only
% Output
% - updateFlag  True/False
persistent cnt

assert (isa(freq_uint8, 'uint8'))

if isempty(cnt)
    cnt = uint8(0);
end

% Initialize
updateFlag = false;

if mod(cnt, freq_uint8) == uint8(0)
    updateFlag = true;
    cnt = uint8(1);
else
    updateFlag = false;
    cnt = cnt + uint8(1);
end

end