%% 
x = 0:999;
x = sin(x);
n = 99;
fp = 0.48;
[h0,h1,g0,g1] = firpr2chfb(n,fp);
[h00,h01,g00,g01] = firpr2chfb(n,fp);
[h000,h001,g000,g001] = firpr2chfb(n,fp);
fvtool(h0,1,h1,1,g0,1,g1,1);

x0 = conv(x, h0);
