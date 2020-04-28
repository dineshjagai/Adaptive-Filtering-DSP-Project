%% 
f = [0 0.125 0.125 1];
m = [1 1 0 0];

b1 = fir2(30,f,m);
[h1,w] = freqz(b1,1);

plot(f,m,w/pi,abs(h1))
xlabel('\omega / \pi')
lgs = {'Ideal','fir2 default'};
legend(lgs)

fvtool(b1)