function [y] = multibank(x, gain1, gain2, gain3, gain4)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
n = 99;
fp = 0.48;
[h0,h1,g0,g1] = firpr2chfb(n,fp);

% band 1 -> 0, pi/8
x0 = conv(x, h0); 
x0 = downsample(x0, 2);
x_band1 = gain1.*downsample(conv(downsample(conv(x0, h0), 2), h0),2);
x_band1 = x_band1(50:length(x_band1)-1);

x_band2 = gain2.*downsample(conv(downsample(conv(x0, h0), 2), h1),2);
x_band2 = x_band2(50:length(x_band2)-1);

x_band3 = gain3.*downsample(conv(x0, h1),2);
x_band3 = x_band3(50:length(x_band3)-1);

x_band4 = gain4.*downsample(conv(x, h1), 2);
x_band4 = x_band4(50:length(x_band4)-1);

y_band4 = conv(upsample(x_band4, 2), g1);
y_band1 = conv(upsample(x_band1, 2), g0);
y_band2 = conv(upsample(x_band2, 2), g1);
y_band12 = conv(upsample(y_band1 + y_band2,2), g0);

y_band3 = conv(upsample(x_band3, 2), g1);
y_band123 = y_band12 + [y_band3 zeros(1, length(y_band12)-length(y_band3))];

y_band123_final = conv(upsample(y_band123, 2),g0);
y = conv(upsample(y_band123,2), g0) +[y_band4 zeros(1, length(y_band123_final)-length(y_band4))];
% y = y(110:110+length(x)-1);
end

