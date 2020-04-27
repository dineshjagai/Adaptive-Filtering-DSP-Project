%Main file, used to run our functions
close all
clear variables
clc 

%% -----------------------------------------------------------------------%% 
%Part A
% constants
num_iter = 700000; 
r_range = [0.85, 0.86, 0.87, 0.88, 0.89, 0.90, 0.91, 0.92, 0.93,... 
           0.94, 0.95, 0.96, 0.97, 0.98]; 
       
% w_noise = .2*pi 
% w_desire = pi
f_noise = 0.4; % frequency of sinusoid signal we want to remove(High inference)
f_desire = 0.2; % frequency of sinusoid signal we want
mu = 0.000009; %small mu reauired to converge

% sequences
e = zeros(1,num_iter); % intermediate Sequence
y = zeros(1, num_iter); % output sequence 
a = zeros(1, num_iter); % Coefficients of a 
range = 1:num_iter;
x_noise = 20 * sin(2*pi*f_noise*range); % sinusoidal noise with Strong interference 
x_desired = sin(2*pi*f_desire*range); % signal desired 
x = x_noise + x_desired; 
w = linspace(-pi, pi, num_iter);
% Plot of Signal Generated 
figure

plot(w/pi, abs(fftshift(fft(x)))) % Note the Frequences at w = 0.8pi 
                                  % i.e. f = 0.4Hz and w = 0.6 pi i.e. f =
                                  % 0.3 hZ


% Simple Algorithm to minimize E(y[n]^2) 
r = 0.95; 
for n = 3:num_iter - 1
    e(n) = x(n) + a(n)*x(n - 1) + x(n - 2);
    y(n) = e(n) - r*a(n)*y(n-1) - r*r*y(n - 2); 
    if ((a(n) >= -2) && (a(n) <= 2))
        a(n + 1) = a(n) - mu*y(n)*x(n-1);
    else 
        a(n) = 0;
    end
    
end

a_fin = a(num_iter);

z = exp(1i*w); 
H_adapt = ((1+a_fin*z.^(-1)+z.^(-2))./(1+r*a_fin*z.^(-1)+r^2*z.^(-2))); 
% For visualizing the filter
% fvtool(H_adapt)

figure()
plot(w/pi, 20*log10(abs(H_adapt)));

% Extract the filtered region from y by taking the last 2k samples
% Extract the sample samples  from x for comparison sake
y_filtered_extract = y(num_iter - 2000 :num_iter); 
x_extract = x(num_iter - 2000 : num_iter); 
% Take the fft of x and y and center them
Y_filter = fftshift(fft(y_filtered_extract, 2001));
X_extract = fftshift(fft(x_extract, 2001));
figure
w_filt = linspace(-pi, pi, 2001)/pi;  
plot(w_filt, abs(X_extract))
hold on
plot(w_filt,abs(Y_filter));
hold off
legend('Input Signal', 'Output signal');
figure
plot(1:num_iter, a);

%% -----------------------------------------------------------------------%% 
%Part B
% constants
    num_iter = 70000; 
% keep w in range [-pi to pi] 
    for f = 0:0.03:0.45    
    f_noise =  f; % frequency of sinusoid signal we want to remove(High inference)
    f_desire = 0.2; % frequency of sinusoid signal we want
    mu = 0.000009; %small mu reauired to converge

    % sequences
    e = zeros(1,num_iter); % intermediate Sequence
    y = zeros(1, num_iter); % output sequence 
    a = zeros(1, num_iter); % Coefficients of a 
    range = 1:num_iter;
    x_noise = 20 * sin(2*pi*f_noise*range); % sinusoidal noise with Strong interference 
    x_desired = sin(2*pi*f_desire*range); % signal desired 
    x = x_noise + x_desired; 
    w = linspace(-pi, pi, num_iter);
    % Plot of Signal Generated 
    
    % Note the Frequences at w = 0.8pi i.e. f = 0.4Hz and w = 0.6 pi i.e. f =
    % 0.3 hZ
%     plot(w/pi, abs(fftshift(fft(x))))


    % Simple Algorithm to minimize E(y[n]^2) 
    r = 0.95; 
    for n = 3:num_iter - 1
        e(n) = x(n) + a(n)*x(n - 1) + x(n - 2);
        y(n) = e(n) - r*a(n)*y(n-1) - r*r*y(n - 2); 
        if ((a(n) >= -2) && (a(n) <= 2))
            a(n + 1) = a(n) - mu*y(n)*x(n-1);
        else 
            a(n) = 0;
        end

    end

    a_fin = a(num_iter);

    z = exp(1i*w); 
    H_adapt = ((1+a_fin*z.^(-1)+z.^(-2))./(1+r*a_fin*z.^(-1)+r^2*z.^(-2))); 
    % For visualizing the filter
    % fvtool(H_adapt)

%     figure()
%     plot(w/pi, 20*log10(abs(H_adapt)));

    % Extract the filtered region from y by taking the last 2k samples
    % Extract the sample samples  from x for comparison sake
    y_filtered_extract = y(num_iter - 2000 :num_iter); 
    x_extract = x(num_iter - 2000 : num_iter); 
    % Take the fft of x and y and center them
    Y_filter = fftshift(fft(y_filtered_extract, 2001));
    X_extract = fftshift(fft(x_extract, 2001));
    figure
    w_filt = linspace(-pi, pi, 2001)/pi;  
    plot(w_filt, abs(X_extract))
    hold on
    plot(w_filt,abs(Y_filter));
    hold off
    legend('Input Signal', 'Output signal');
    figure
    plot(1:num_iter, a);
    end
%% -----------------------------------------------------------------------%% 
%Part C
% constants
r_one = 0.95; 
r_two = 0.94;
mu_one = 0.000009; %small mu reauired to converge
mu_two = 0.000009; %small mu reauired to converge
num_iter = 700000; 
f_noise_one = 0.4; % frequency of sinusoid signal we want to remove(High inference)
f_noise_two = 0.1; % frequency of second sinusoid signal we want to remove(High inference)
f_desire = 0.2; % frequency of sinusoid signal we want

% sequences
e_one = zeros(1,num_iter); % intermediate Sequence one
e_two = zeros(1,num_iter); % intermediate Sequence two
y_one = zeros(1, num_iter); % output sequence one
y_two = zeros(1, num_iter); % output sequence two
a_one = zeros(1, num_iter); % Coefficients of a one
a_two = zeros(1, num_iter); % Coefficients of a two
range = 1:num_iter;
x_noise_one = 20 * sin(2*pi*f_noise_one*range); % sinusoidal one noise with Strong interference 
x_noise_two = 20 * sin(2*pi*f_noise_two*range); % sinusoidal two noise with Strong interference 
x_desired = sin(2*pi*f_desire*range); % signal desired 
x = x_noise_one + x_desired;
x_full = x_noise_one + x_noise_two + x_desired;
w = linspace(-pi, pi, num_iter);

% Plot of Signal Generated 
figure
plot(w/pi, abs(fftshift(fft(x))))

% Simple Algorithm to minimize E(y[n]^2) 

for n = 3:num_iter - 1
    e_one(n) = x(n) + a_one(n)*x(n - 1) + x(n - 2);
    y_one(n) = e_one(n) - r_one*a_one(n)*y_one(n-1) - r_one*r_one*y_one(n - 2); 
    if ((a_one(n) >= -2) && (a_one(n) <= 2))
        a_one(n + 1) = a_one(n) - mu_one*y_one(n)*x(n-1);
    else 
        a_one(n) = 0;
    end    
end
z = exp(1i*w); 
a_fin = a_one(num_iter);
H_adapt_one = ((1+a_fin*z.^(-1)+z.^(-2))./(1+r_one*a_fin*z.^(-1)+r_one^2*z.^(-2))); 
x_two = x_noise_two + y_one;

for n = 3:num_iter - 1
    e_two(n) = x_two(n) + a_two(n)*x_two(n - 1) + x_two(n - 2);
    y_two(n) = e_two(n) - r_two*a_two(n)*y_two(n-1) - r_two*r_two*y_two(n - 2); 
    if ((a_two(n) >= -2) && (a_two(n) <= 2))
        a_two(n + 1) = a_two(n) - mu_two*y_two(n)*x_two(n-1);
    else 
        a_two(n) = 0;
    end
end
    
a_fin = a_two(num_iter);
H_adapt_two = ((1+a_fin*z.^(-1)+z.^(-2))./(1+r_two*a_fin*z.^(-1)+r_two^2*z.^(-2))); 
% For visualizing the filter
% fvtool(H_adapt)
H_adapt = H_adapt_one .* H_adapt_two;
figure()
plot(w/pi, 20*log10(abs(H_adapt)));

% Extract the filtered region from y by taking the last 2k samples
% Extract the sample samples  from x for comparison sake
y_filtered_extract = y_two(num_iter - 2000 :num_iter); 
x_extract = x_full(num_iter - 2000 : num_iter); 
% Take the fft of x and y and center them
Y_filter = fftshift(fft(y_filtered_extract, 2001));
X_extract = fftshift(fft(x_extract, 2001));
figure
w_filt = linspace(-pi, pi, 2001)/pi;  
plot(w_filt, abs(X_extract))
hold on
plot(w_filt,abs(Y_filter));
hold off
legend('Input Signal', 'Output signal');
figure
plot(1:num_iter, a_one);
figure
plot(1:num_iter, a_two);












