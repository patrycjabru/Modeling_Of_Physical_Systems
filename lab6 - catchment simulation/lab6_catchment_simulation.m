clc
close all
clear 

data = importdata('L06\opady.prn');
realData = importdata('L06\dunaj.prn');

C_in = data(:,2);

t_half = 12*12+3.6;
lambda = log(2)/t_half;

tt = 7;
T = 629;
dt = 1;

Pe = 1.8;

for i = 161:T
    C(i) = exponentional_integral(C_in,i,dt,tt,lambda);
    C2(i) = dispersion_integral(C_in,i,dt,tt,lambda,Pe);
end

figure
%plot(C());
title('Concentration of tritum in months')
ylabel('Tritum concentration')
xlabel('Month')
hold on
plot(C2());
plot(realData(:,2));
legend('dispersion model','real data')

function total = exponentional_integral(c_in, i, dt, tt, lambda)
    sum = 0;
    t = i * dt;
    for j = 1:i-1
        tp = j*dt;
        sum = sum + c_in(j) * ...
                    tt^(-1) * ...
                    exp(-1 * (t - tp) / tt) * ... 
                    exp(-1 * lambda * (t - tp));
    end
    total = sum * dt;
end

function total = dispersion_integral(c_in, i, dt, tt, lambda, Pe)
    sum = 0;
    t = i * dt;
    for j = 1:i-1
        tp = j*dt;
        a = (4*pi*Pe*(t - tp)/tt)^(-1/2);
        b = exp((-(1-(t-tp)/tt)^2)/(4*Pe*(t-tp)/tt));
        sum = sum + c_in(j) * ...
                    tt^(-1) * ...
                    a * 1/((t-tp)*b) * ... 
                    exp(-1 * lambda * (t - tp));
    end
    total = sum * dt;
end