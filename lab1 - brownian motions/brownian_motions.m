clear;
clc;
close all;

Time = 400;
N=5000;
dt = 1;

xN=zeros(N,Time);
yN=zeros(N,Time);
a2N=zeros(N,Time);
for n=1:N
    x=zeros(1,Time);
    y=zeros(1,Time);
    a2=zeros(1,Time);
    for t=2:dt:Time
       x(t) = x(t-1)+randn(1);
       y(t) = y(t-1)+randn(1);
       a2(t) = (x(t)^2+y(t)^2);
    end
    xN(n,:) = x;
    yN(n,:) = y;
    a2N(n,:) = a2;
end
figure(1)
plot(xN(1:6,:)',yN(1:6,:)')
title('Position of particles');
xlabel('x');
ylabel('y');

%autocorelation of coordinatex x and y of the last particle
figure(2)
autocorr(x)
title('Autocorelation of x coordiate of one particle');
figure(3)
autocorr(y)
title('Autocorelation of y coordiate of one particle');

%autocorelation of random generated nnumbers 
r = randn(1000,1);
figure
autocorr(r)
title('Autocorelation of random generated numbers');

meanVector = zeros(1,Time);
for i=1:dt:Time
    column = a2N(:,i);
    meanVector(i) = mean(column);
end
figure
plot(meanVector);
title('Relationship between mean square of displacement and time');
xlabel('Time');
ylabel('Mean square of displacement');
