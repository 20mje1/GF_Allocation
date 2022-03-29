clear;
close;
cfgVoIP = networkTrafficVoIP;
cfgVoIP.ExponentialMean=1250;
cfgVoIP.HasJitter=1;
% Time interval between two consecutive packet transfers in milliseconds
dt = zeros(100, 1);

% Packet sizes in bytes
packetSize = zeros(100, 1);

for packetCount = 1:100
    [dt(packetCount),packetSize(packetCount)] = generate(cfgVoIP);
end
stem(packetSize); 
title('Packet Size Versus Packet Number');
xlabel('Packet Number');
ylabel('Packet Size in Bytes');