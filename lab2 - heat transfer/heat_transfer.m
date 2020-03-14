clc
clear all
close all

a = 30;
b = 20;
c = 5;
d = 10;
thickness = 2; 

plate_init_temp = 10;
heater_init_temp = 80;

shape_array = zeros(a,a);
shape_array(1,:) = 2;
shape_array(1:d,1) = 2;
shape_array(1:d,a) = 2;
shape_array(d:a,c+1) = 2;
shape_array(d:a,a-c) = 2;
shape_array(d,1:c+1) = 2;
shape_array(d,a-c:a) = 2;
shape_array(a,c:a-c) = 2;

shape_array(d+1:a,1:c) = -1;
shape_array(d+1:a,a-c+1:a) = -1;

shape_array(d/2,(a-d)/2+1:(a+d)/2) = 1;
shape_array(3*d/2,(a-d)/2+1:(a+d)/2) = 1;
shape_array(d/2:3*d/2,(a-d)/2+1) = 1;
shape_array(d/2:3*d/2,(a+d)/2) = 1;

shape_array(round(d/2)+1:round(3*d/2)-1,(a-d)/2+2:(a+d)/2-1) = 1;

T = 1000;
K = 237;
c = 900;
p = 2700;
dt = 1;
dx = 0.05;
dy = 0.05;

factor = K/(c*p);

temp = zeros(a,a,T);
temp(d/2:3*d/2,(a-d)/2+1:(a+d)/2,:) = heater_init_temp;

for i=1:a
    for j=1:a
        if shape_array(i,j)==0
            temp(i,j,1)=plate_init_temp;
        elseif shape_array(i,j)==2
            temp(i,j,1)=plate_init_temp;
        end
    end
end

for t=2:T
    for i=1:a
        for j=1:a
            if shape_array(i,j)==2
                temp(i,j,t)=plate_init_temp;
            elseif shape_array(i,j)==0
                temp(i,j,t)=temp(i,j,t-1)+(factor*dt)/dx^2 * (temp(i+1,j,t-1) - 2*temp(i,j,t-1) + temp(i-1,j,t-1)) + (factor*dt)/dy^2 * (temp(i,j+1,t-1) - 2*temp(i,j,t-1) + temp(i,j-1,t-1));
            end   
        end
    end
end

surf(temp(:,:,T))