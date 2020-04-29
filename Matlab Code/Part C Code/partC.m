%% 
n = 0:99;
w = 0.2*pi;
w1 = 0.4*pi;
w2 = 0.1*pi;
x = sin(w.*n) + sin(w1.*pi);

out = multibank(x, 1, 1, 1, 1);


figure
stem(x)
hold on
stem(out)
legend('Input, x[n]', 'Output, y[n]')
xlabel('n')
title('Output vs Input of Multirate Filter Bank')

%% 

speech_out = multibank(nspeech2', 1, 1, 1, 1);

figure 
plot(linspace(-pi, pi, length(nspeech2'))./pi, abs(fft(nspeech2')))
hold on
plot(linspace(-pi, pi, length(speech_out))./pi, abs(fft(speech_out)))
legend('Input Speech', 'Output Speech')
xlabel('w (rad/s) (Normalized to pi)')
title('Effect of Filter Bank on Speech')
