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
shape_array(1,:) = 20;
shape_array(1:d,1) = 21;
shape_array(1:d,a) = 22;
shape_array(d:a,c+1) = 21;
shape_array(d:a,a-c) = 22;
shape_array(d,1:c+1) = 23;
shape_array(d,a-c:a) = 23;
shape_array(a,c:a-c) = 23;

shape_array(d+1:a,1:c) = -1;
shape_array(d+1:a,a-c+1:a) = -1;

shape_array(d/2,(a-d)/2+1:(a+d)/2) = 1;
shape_array(3*d/2,(a-d)/2+1:(a+d)/2) = 1;
shape_array(d/2:3*d/2,(a-d)/2+1) = 1;
shape_array(d/2:3*d/2,(a+d)/2) = 1;

shape_array(round(d/2)+1:round(3*d/2)-1,(a-d)/2+2:(a+d)/2-1) = 1;

shape_array(d/2-1,(a-d)/2+1:(a+d)/2) = 3;
shape_array(3*d/2+1,(a-d)/2+1:(a+d)/2) = 4;
shape_array(d/2:3*d/2,(a-d)/2+1-1) = 5;
shape_array(d/2:3*d/2,(a+d)/2+1) = 6;

T = 1000;
%alumina
%K = 237;
%cw = 900;
%p = 2700;

%cooper
%K = 401;
%cw = 380;
%p = 8920;

%stainless steal
K = 58;
cw = 450;
p = 7860;

dx = 0.01;
dy = 0.01;
%dt = 0.27;
dt = 0.25*min(dx,dy)^2*(cw*p)/K;

factor = K/(cw*p);

temp = zeros(a,a,T);
temp(d/2:3*d/2,(a-d)/2+1:(a+d)/2,:) = heater_init_temp;

%Boundary condition type 1
for i=1:a
    for j=1:a
        if shape_array(i,j)==0 || shape_array(i,j)==3 || shape_array(i,j)==4 || shape_array(i,j)==5 || shape_array(i,j)==6
            temp(i,j,1)=plate_init_temp;
        elseif shape_array(i,j)==20 || shape_array(i,j)==21 || shape_array(i,j)==22 || shape_array(i,j)==23 || shape_array(i,j)==24
            temp(i,j,1)=plate_init_temp;
        end
    end
end

didChange = false; 
minimumChange = 0.001;
for t=2:T
    for i=1:a
        for j=1:a
            if shape_array(i,j)==20 || shape_array(i,j)==21 || shape_array(i,j)==22 || shape_array(i,j)==23 || shape_array(i,j)==24
                temp(i,j,t)=plate_init_temp;
            elseif shape_array(i,j)==0 || shape_array(i,j)==3 || shape_array(i,j)==4 || shape_array(i,j)==5 || shape_array(i,j)==6
                temp(i,j,t)=temp(i,j,t-1)+(factor*dt)/dx^2 * (temp(i+1,j,t-1) - 2*temp(i,j,t-1) + temp(i-1,j,t-1)) + (factor*dt)/dy^2 * (temp(i,j+1,t-1) - 2*temp(i,j,t-1) + temp(i,j-1,t-1));
            end
            if abs(temp(i,j,t)-temp(i,j,t-1)) > minimumChange
                didChange = true;
            end
        end
    end
    if ~didChange
        break;
    end
    didChange = false;
end

figure
surf(temp(:,:,t))
title("Temperature distribution in a plate");
xlabel("X");
ylabel("Y");
zlabel("Temperature");
disp("Simulation 1 time: ")
disp(t*dt)

%figure
%for i=1:t
%    surf(temp(:,:,i))
%    zlim([0 100]);
%    pause(0.01)
%end

%Boundary condition type 2

shape_array(1,1) = -2;
shape_array(1,a) = -2;
shape_array(d,1) = -2;
shape_array(d,a) = -2;
shape_array(a,c+1) = -2;
shape_array(a,a-c) = -2;


temp = zeros(a,a,T);
plate_init_temp = 20;

for i=1:a
    for j=1:a
        if shape_array(i,j)==0 || shape_array(i,j)==3 || shape_array(i,j)==4 || shape_array(i,j)==5 || shape_array(i,j)==6
            temp(i,j,1)=plate_init_temp;
        elseif shape_array(i,j)==20 || shape_array(i,j)==21 || shape_array(i,j)==22 || shape_array(i,j)==23 || shape_array(i,j)==24
            temp(i,j,1)=plate_init_temp;
        elseif shape_array(i,j)==-2
            temp(i,j,1)=plate_init_temp;
        end
    end
end

P = 100;
heatingTime = 10/dt; 
minimumChange = 0.001;
for t=2:T
    for i=1:a
        for j=1:a
            %inner part of plate
            if shape_array(i,j)==0
                temp(i,j,t) = calculateTemp(i,j,t,temp,factor,dt,dx,dy);
            %heater border
            elseif shape_array(i,j)==5 
                if t<heatingTime
                    temp(i,j,t) = calculateTempFromPower(i,j,t,temp,P,dt,cw,p,d,thickness);
                else
                    temp(i,j,t) = calculateTempWithoutRightNode(i,j,t,temp,factor,dt,dx,dy);
                end
            elseif shape_array(i,j)==6 
                if t<heatingTime
                    temp(i,j,t) = calculateTempFromPower(i,j,t,temp,P,dt,cw,p,d,thickness);
                else
                    temp(i,j,t) = calculateTempWithoutLeftNode(i,j,t,temp,factor,dt,dx,dy);
                end
            elseif shape_array(i,j)==3
                if t<heatingTime
                    temp(i,j,t) = calculateTempFromPower(i,j,t,temp,P,dt,cw,p,d,thickness);
                else
                    temp(i,j,t) = calculateTempWithoutLowerNode(i,j,t,temp,factor,dt,dx,dy);
                end
            elseif shape_array(i,j)==4
                if t<heatingTime
                    temp(i,j,t) = calculateTempFromPower(i,j,t,temp,P,dt,cw,p,d,thickness);
                else
                    temp(i,j,t) = calculateTempWithoutUpperNode(i,j,t,temp,factor,dt,dx,dy);
                end
            %plate border
            elseif shape_array(i,j)==22
                temp(i,j,t) = calculateTempWithoutRightNode(i,j,t,temp,factor,dt,dx,dy);
            elseif shape_array(i,j)==21
                temp(i,j,t) = calculateTempWithoutLeftNode(i,j,t,temp,factor,dt,dx,dy);
            elseif shape_array(i,j)==23
                temp(i,j,t) = calculateTempWithoutLowerNode(i,j,t,temp,factor,dt,dx,dy);
            elseif shape_array(i,j)==20
                temp(i,j,t) = calculateTempWithoutUpperNode(i,j,t,temp,factor,dt,dx,dy);
            %plate corners
            elseif shape_array(i,j)==-2
                temp(i,j,t) = copyTemp(i,j,t,temp,a);
            end
            if abs(temp(i,j,t)-temp(i,j,t-1)) > minimumChange
                didChange = true;
            end
        end
    end
    if ~didChange
        break;
    end
    didChange = false;
end

figure
surf(temp(:,:,t))
title("Temperature distribution in a plate");
xlabel("X");
ylabel("Y");
disp("Simulation 2 time: ")
disp(t*dt)
disp("delta temp")
disp(calculateMeanTemp(temp,t,a)-20)

%figure
%for i=1:t
%    surf(temp(:,:,i))
%    zlim([0 40])
%    pause(0.01)
%end

function value = calculateTemp(i,j,t,temp,factor,dt,dx,dy)
    value = temp(i,j,t-1)+(factor*dt)/dx^2 * (temp(i+1,j,t-1) - 2*temp(i,j,t-1) + temp(i-1,j,t-1)) + (factor*dt)/dy^2 * (temp(i,j+1,t-1) - 2*temp(i,j,t-1) + temp(i,j-1,t-1));
end

function value = calculateTempWithoutRightNode(i,j,t,temp,factor,dt,dx,dy)
    value = temp(i,j,t-1)+(factor*dt)/dx^2 * (temp(i+1,j,t-1) - 2*temp(i,j,t-1) + temp(i-1,j,t-1)) + (factor*dt)/dy^2 * ( - temp(i,j,t-1) + temp(i,j-1,t-1));
end

function value = calculateTempWithoutLeftNode(i,j,t,temp,factor,dt,dx,dy)
    value = temp(i,j,t-1)+(factor*dt)/dx^2 * (temp(i+1,j,t-1) - 2*temp(i,j,t-1) + temp(i-1,j,t-1)) + (factor*dt)/dy^2 * (temp(i,j+1,t-1) - temp(i,j,t-1));
end

function value = calculateTempWithoutLowerNode(i,j,t,temp,factor,dt,dx,dy)
    value = temp(i,j,t-1)+(factor*dt)/dx^2 * ( - temp(i,j,t-1) + temp(i-1,j,t-1)) + (factor*dt)/dy^2 * (temp(i,j+1,t-1) - 2*temp(i,j,t-1) + temp(i,j-1,t-1));
end

function value = calculateTempWithoutUpperNode(i,j,t,temp,factor,dt,dx,dy)
    value = temp(i,j,t-1)+(factor*dt)/dx^2 * (temp(i+1,j,t-1) - temp(i,j,t-1)) + (factor*dt)/dy^2 * (temp(i,j+1,t-1) - 2*temp(i,j,t-1) + temp(i,j-1,t-1));
end

function value = calculateTempFromPower(i,j,t,temp,P,dt,cw,p,d,thickness)
    value = temp(i,j,t-1) + (P * dt) / (cw * (d/100)^2 * (thickness/1000) * p);
end

function value = copyTemp(i,j,t,temp,a)
    if i-1>0 && temp(i-1,j,t-1)~=0
        value = temp(i-1,j,t-1);
    elseif i+1<a+1 && temp(i+1,j,t-1)~=0
        value = temp(i+1,j,t-1);
    elseif j-1>0 && temp(i,j-1,t-1)~=0
        value = temp(i,j-1,t-1);
    elseif j+1<a+1 && temp(i,j+1,t-1)~=0
        value = temp(i,j+1,t-1);
    end
end

function value = calculateMeanTemp(temp,t,a) 
    numberOfNodes = 0;
    sumTemp = 0;
    for i=1:a
       for j=1:a
           if temp(i,j,t)>0
              numberOfNodes = numberOfNodes + 1;
               sumTemp = sumTemp + temp(i,j,t);
           end    
       end
    end
    value = sumTemp/numberOfNodes;
end