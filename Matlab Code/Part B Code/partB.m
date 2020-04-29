%% Generate training signal
n_in = 0:999;
train = randi([0 1], 1, 1000);
train(train == 0) = -1;

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

%%  Train Adaptive Filter
SNR = 20;
filter_length = 12;
h_filter = ones(1, filter_length);
out = conv(train, h);
out = awgn(out, SNR);
max_mu = 2 / (norm(out)^2);
mu = 0.05;

desired = [zeros(1, 2) train];

for i=length(h_filter):length(train)
   outn = out(i-length(h_filter)+1:i);
   s_hat = outn*h_filter';
   e = desired(i) - s_hat;
   h_filter = h_filter + mu*e*outn;
end

%% Test Adaptive Filter

H_channel=freqz(h);       
w=(0:length(H_channel)-1)./(length(H_channel));         
H_channel_inv_mag=-10*log10(real(H_channel).^2 + imag(H_channel).^2);  
H_channel_inv_phase=-imag(H_channel)./real(H_channel);                    
H_adaptive=freqz(h_filter);       
H_adaptive_mag=10*log10(real(H_adaptive).^2 + imag(H_adaptive).^2);    
H_adaptive_phase=imag(H_adaptive)./real(H_adaptive);   

figure 
plot(w,H_channel_inv_mag)
hold on
plot(w,H_adaptive_mag)
legend('Desired Spectrum','Adaptive Spectrum');
xlabel('w (rad/s)');
ylabel('Magnitude (dB)');
title('Magnitude response');

%% Plots of channel and filter

fvtool(h);
fvtool(h_filter);


%% Overall LTI system

h_overall = conv(h, h_filter);
fvtool(h_overall);

%% 

test = randi([0 1], 1, 3000);
test(test == 0) = -1;

channel_out = conv(test, h);
acc_channel = 0;
for i=1:length(test)
    if (channel_out(i) > 0)
        channel_out(i) = 1;
    else
        channel_out(i) = -1;
    end
    
    if (channel_out(i) == test(i)) 
        acc_channel = acc_channel + 1;
    end
end
acc_channel = acc_channel / length(test);
equalized = conv(channel_out, h_filter);

acc_eq = 0;
for i=1:length(test)
    if (equalized(i) > 0)
        equalized(i) = 1;
    else
        equalized(i) = -1;
    end
    
    if (equalized(i) == test(i)) 
        acc_eq = acc_eq + 1;
    end
end
acc_eq = acc_eq / length(test);
