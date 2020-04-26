close all
clear all
clc
dt = 0.1;
dx = 0.1;

length = 100; %m
width = 5; %m
depth = 1; %m
x_size = length / dx;
T = 100; %units

D = 0.01;
U = 0.1;

Ca = (U*dt)/dx;
Cd = (D*dt)/dx^2;

c = zeros(x_size,T);
tracer_mass = zeros(T);

injection_point = 10; %m
injection_point_x = injection_point/dx;
injected_tracer = 10; %kg

c(injection_point_x,1) = injected_tracer/(width*depth*dx);

F1 = Cd*(1-Ca) - Ca/6*(Ca^2-3*Ca+2);
F2 = Cd*(2-3*Ca) - Ca/2*(Ca^2-2*Ca-1);
F3 = Cd*(1-3*Ca) - Ca/2*(Ca^2-Ca-2);
F4 = Cd*Ca + Ca/6*(Ca^2-1);
for n=2:T
    for j=1:x_size-1
        if j==1 || j==2
            c(j,n)=0;
        else
            c(j,n) = c(j,n-1) + F1*c(j+1,n-1) - F2*c(j,n-1) + F3*c(j-1,n-1) + F4*c(j-2,n-1);
        end
        tracer_mass(n) = tracer_mass(n) + c(j,n);
    end 
end

figure
plot(c(:,T));
xlabel("Node coordinate")
ylabel("Tracer mass")
title("Tracer distribution in the river")

figure
plot(tracer_mass)
title("Tracer mass over time")
xlabel("Time steps")
ylabel("Mass [kg]")
