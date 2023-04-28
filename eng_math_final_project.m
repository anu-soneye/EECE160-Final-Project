close all 


%% One Story Test

%Parameters
story_1 = [40 3 10]; %Mass (kg), Spring Constant (N/m), Damping (kg/s)
%End of parameters

%Solution in the form x(t) = (v1)e^(st) + (v2)e^(rt)
A = [0 1 ; (1/story_1(1)) * -story_1(2), (1/story_1(1)) * -story_1(3)];

%Eigenvalue in terms of a + Bi where a is the decay rate and B is the
%angular frequency of the mass spring system
eigen_val = eig(A);
nat_freq_1 = imag(eigen_val(1));

tspan = [0 300];
x0 = [20; 0;];
[t,x] = ode45(@(t,x) odefcn(t,x,story_1, 20 * cos(nat_freq_1 * t)), tspan, x0); 
[t2,x2] = ode45(@(t2,x2) odefcn(t2,x2,story_1, 20 * cos((nat_freq_1*2) * t2)), tspan, x0); 
% in this case, the natural angular frequency is 1.4 = 2 * pi * f (Hz = s^-1)

plot(t,x(:,1),'-o', t2,x2(:,1),'r-o')
title ("One Story Building (M=40 kg, K=3N/m, C=10kg/s)");
legend("frequency = natural", "frequency = 2 * natural");
function dxdt = odefcn(t,x,story_1, F1)
  m1 = story_1(1);
  k1 = story_1(2);
  c1 = story_1(3);
  
  dxdt = zeros(2,1);
  dxdt(1) = x(2);
  dxdt(2) = (1/m1) * (-c1*x(2) - k1 * x(1) + F1);
end


%% Two Story Test
%{
%Parameters
story_1 = [10000 100 0]; %Mass (kg), Spring Constant (N/m), Damping (kg/s)
story_2 = [10000 100 0];
%End of parameters

 m1 = story_1(1);
 k1 = story_1(2);
 c1 = story_1(3);
 
 m2 = story_2(1);
 k2 = story_2(2);
 c2 = story_2(3);
  
%Solution in the form x(t) = (v1)e^(st) + (v2)e^(rt)
A = [0 0 1 0 ; 0 0 0 1; -(k1 + k2)/m1 k2/m1 -(c1+c2)/m1 c2/m1; k2/m2 -k2/m2 c2/m2 -c2/m2];

%Eigenvalue in terms of a + Bi where a is the decay rate and B is the
%angular frequency of the mass spring system
eigen_val = eig(A);
nat_freq_1 = imag(eigen_val(1));
nat_freq_2 = imag(eigen_val(3));

tspan = [0 300];
x0 = [0; 30; 0; 0];
[t,x] = ode45(@(t,x) odefcn(t,x,story_1, story_2, 50 * cos(nat_freq_1 * t)), tspan, x0);
[t2,x2] = ode45(@(t2,x2) odefcn(t2,x2,story_1, story_2, 50 * cos(nat_freq_2 * t2)), tspan, x0);

plot(t,x(:,1),'-r',t,x(:,2))

function dxdt = odefcn(t,x,story_1, story_2, F1)
  m1 = story_1(1);
  k1 = story_1(2);
  c1 = story_1(3);
  m2 = story_2(1);
  k2 = story_2(2);
  c2 = story_2(3);
  
  dxdt = zeros(4,1);
  dxdt(1) = x(3);
  dxdt(2) = x(4);
  dxdt(3) = (1/m1) * (-k1 * x(1) + k2 * (x(2) - x(1)) - c1 * x(3) + c2 * (x(4) - x(3)) + F1);
  dxdt(4) = (1/m2) * (-k2 * (x(2) - x(1)) - c2 * (x(4) - x(3)));

end

%}
%{
%parameters
A=1;
B=1;
C=1;
%time interval and initial conditions
t_interval = [0,10];
init_cond = [0,0,0,0]';
%solution
[t,y] = ode45(@(t,Y) odefcn(t,Y,A,B,C) , t_interval , init_cond);
%plot
plot(t,y(:,1),'b',t,y(:,2),'r');

function dYdt = odefcn(t,Y,A,B,C)
dYdt = [ Y(3);
         Y(4);
         A*(Y(3)-Y(4)) + B*Y(1)^3;
         C*(Y(4)-Y(3)) + B*Y(2)^3];
end
%}