clc
clear all
close all

A = 0.3;
S = 1366;
%Earth_area = 510072000*1000*1000;
Stefan_Boltzmann_constant = 5.67 * 10^(-8);

Temperature_without_atmosphere = (S*(1-A)/(4*Stefan_Boltzmann_constant))^(1/4);
Temperature_without_atmosphere_celsius = Temperature_without_atmosphere - 273.15;

disp("Mean temperature of the Earth without atmosphere (Celsius degrees)");
disp(Temperature_without_atmosphere_celsius);

x0 = [100,100];
i=1;
for s=0.8*S:1.2*S
    fun = @(x) root2d(x,s);
    x = fsolve(fun,x0);
    output(i,1) = x(1);
    output(i,2) = x(2);
    i = i+1;
end
figure
plot(0.8*S:1.2*S,output(:,1))
xlabel("surface temperature (Kelvins)")
ylabel("solar constant")
title("Relationshio between surface temperature and solar constant")
%disp("Mean temperature of the Earth surface and atmosphere (celsius degrees)");
%disp(x(1)-273.15);
%disp(x(2)-273.15);

function F = root2d(x,S)
a_s = 0.19; %a_s
t_a = 0.53; %t_a
a_a = 0.3; %a_a

t_a_p = 0.06; %t_a'
a_a_p = 0.31; %a_a'

c = 2.7;
Stefan_Boltzmann_constant = 5.67 * 10^(-8);

F(1) = (-t_a)*(1-a_s)*S/4 + c*(x(1) - x(2)) + Stefan_Boltzmann_constant*x(1)^4*(1-a_a_p)-Stefan_Boltzmann_constant*x(2)^4;
F(2) = -(1-a_a-t_a+a_s*t_a)*S/4 - c*(x(1) - x(2)) - Stefan_Boltzmann_constant*x(1)^4*(1-t_a_p-a_a_p)+2*Stefan_Boltzmann_constant*x(2)^4;
end 

