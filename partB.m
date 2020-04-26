%% Setup
amplitudes = [-1 1];
n_in = 0:999;
train = zeros(1, 1000);
for i=1:1000
    val = randi(2);
    train(i) = amplitudes(val);
end

h = [0.3 1 0.7 0.3 0.2];
x = conv(train, h);
x = awgn(x, 37);
n_out = 1:length(x);

figure
stem(n_in, train)
xlabel('n')
ylabel('s[n]')
title('Input Signal')

figure
stem(n_out, x);
xlabel('n')
ylabel('x[n]')
title('Output Signal w/ Interference')

%% 

learning_rate = 0.01;
h_filter = zeros(1, length(h));
delay = zeros


