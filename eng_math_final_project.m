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

tspan = [0 150];
x0 = [20; 0;];
[t,x] = ode45(@(t,x) odefcn1(t,x,story_1, 20 * cos(nat_freq_1 * t)), tspan, x0); 
[t2,x2] = ode45(@(t2,x2) odefcn1(t2,x2,story_1, 30 * cos((nat_freq_1) * t2)), tspan, x0); 
[t3,x3] = ode45(@(t3,x3) odefcn1(t3,x3,story_1, 20 * cos((nat_freq_1*2) * t3)), tspan, x0); 
% in this case, the natural angular frequency is 1.4 = 2 * pi * f (Hz = s^-1)

plot(t,x(:,1),'-o', t2,x2(:,1),'g-o', t3,x3(:,1),'r-o')
title('One Story Building (M=40 kg, K=3 N/m, C=10 kg/s)');
legend('f = natural = ' + string(nat_freq_1) + ", A = 20",'f = natural = ' + string(nat_freq_1)+ ", A = 30", 'f = 2 * natural = ' + string(2 * nat_freq_1)+ ", A = 20");


%% Two Story Test

%{
%Parameters
story_1 = [10000 100 10]; %Mass (kg), Spring Constant (N/m), Damping (kg/s)
story_2 = [10000 100 10];
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

tspan = [0 1000];
x0 = [0; 30; 0; 0];
[t,x] = ode45(@(t,x) odefcn2(t,x,story_1, story_2, 200 * cos(nat_freq_1 * t)), tspan, x0);
[t2,x2] = ode45(@(t2,x2) odefcn2(t2,x2,story_1, story_2, 300 * cos(nat_freq_1 * t2)), tspan, x0);
[t3,x3] = ode45(@(t3,x3) odefcn2(t3,x3,story_1, story_2, 200 * cos(nat_freq_2 * t3)), tspan, x0);
[t4,x4] = ode45(@(t4,x4) odefcn2(t4,x4,story_1, story_2, 200 * cos((nat_freq_1 + nat_freq_2) * t4)), tspan, x0);
[t5,x5] = ode45(@(t5,x5) odefcn2(t5,x5,story_1, story_2, 200 * cos((1) * t5)), tspan, x0);

plot(t,x(:,1),t2,x2(:,1), t3,x3(:,1), t4,x4(:,1), t5,x5(:,1));
title('Floor 1 of Two Story Building (M1, M2 =10000 kg | K1, K2=100 N/m | C1, C2=10 kg/s)');
legend('floor1: f = nat_1 = ' + string(nat_freq_1) + ", A = 200", 'floor1: f = nat_1 = ' + string(nat_freq_1)+ ", A = 300",'floor1: f = nat_2 = ' + string(nat_freq_2)+ ", A = 200", 'floor1: f = nat_1 + nat_2 = ' + string(nat_freq_1 + nat_freq_2)+ ", A = 200", 'floor1: f = 1' + ", A = 200");

figure;
plot(t,x(:,2),t2,x2(:,2),t3,x3(:,2),t4,x4(:,2), t5,x5(:,2));
title('Floor 2 of Two Story Building (M1, M2 =10000 kg | K1, K2=100 N/m | C1, C2=10 kg/s)');
legend('floor2: f = nat_1 = ' + string(nat_freq_1) + ", A = 200", 'floor2: f = nat_1 = ' + string(nat_freq_1)+ ", A = 300",'floor2: f = nat_2 = ' + string(nat_freq_2)+ ", A = 200", 'floor2: f = nat_1 + nat_2 = ' + string(nat_freq_1 + nat_freq_2)+ ", A = 200", 'floor2: f = 1' + ", A = 200");

max_1 = max(x)
max_2 = max(x2)
max_3 = max(x3)
max_4 = max(x4)
max_5 = max(x5)


%}


%% Three Story Test

%{
%Parameters
story_1 = [10000 100 10]; %Mass (kg), Spring Constant (N/m), Damping (kg/s)
story_2 = [10000 100 10];
story_3 = [10000 100 10];
%End of parameters

 m1 = story_1(1);
 k1 = story_1(2);
 c1 = story_1(3);
 
 m2 = story_2(1);
 k2 = story_2(2);
 c2 = story_2(3);
 
 m3 = story_3(1);
 k3 = story_3(2);
 c3 = story_3(3);
  
%Solution in the form x(t) = (v1)e^(st) + (v2)e^(rt)
A = [0 0 0 1 0 0; 0 0 0 0 1 0; 0 0 0 0 0 1; -(k1 + k2)/m1, k2/m1, 0, -(c1+c2)/m1, c2/m1, 0; k2/m2, -(k2 + k3)/m2, k3/m2, c2/m2, -(c2 + c3)/m2, c3/m2; 0, k3/m3, -k3/m3, 0, c3/m3, -c3/m3];

%Eigenvalue in terms of a + Bi where a is the decay rate and B is the
%angular frequency of the mass spring system
eigen_val = eig(A);
nat_freq_1 = imag(eigen_val(1));
nat_freq_2 = imag(eigen_val(3));
nat_freq_3 = imag(eigen_val(5));

tspan = [0 2000];
x0 = [0; 30; 60; 0; 0; 0];
[t,x] = ode45(@(t,x) odefcn3(t,x,story_1, story_2, story_3, 200 * cos(nat_freq_1 * t)), tspan, x0);
[t2,x2] = ode45(@(t2,x2) odefcn3(t2,x2,story_1, story_2, story_3, 300 * cos(nat_freq_1 * t2)), tspan, x0);
[t3,x3] = ode45(@(t3,x3) odefcn3(t3,x3,story_1, story_2, story_3, 200 * cos(nat_freq_2 * t3)), tspan, x0);
[t4,x4] = ode45(@(t4,x4) odefcn3(t4,x4,story_1, story_2, story_3, 200 * cos(nat_freq_3 * t4)), tspan, x0);
[t5,x5] = ode45(@(t5,x5) odefcn3(t5,x5,story_1, story_2, story_3, 200 * cos((nat_freq_1 + nat_freq_2 + nat_freq_3) * t5)), tspan, x0);

plot(t,x(:,1),t2,x2(:,1), t3,x3(:,1), t4,x4(:,1), t5,x5(:,1));
title('Floor 1 of Three Story Building (M1, M2, M3 =10000 kg | K1, K2, K3=100 N/m | C1, C2, C3=10 kg/s)');
legend('floor1: f = nat_1 = ' + string(nat_freq_1) + ", A = 200", 'floor1: f = nat_1 = ' + string(nat_freq_1)+ ", A = 300",'floor1: f = nat_2 = ' + string(nat_freq_2)+ ", A = 200", 'floor1: f = nat_3 = ' + string(nat_freq_3)+ ", A = 200", 'floor1: f = nat_1 + nat_2 + nat_3 = ' + string(nat_freq_1 + nat_freq_2 + nat_freq_3) + ", A = 200");

figure;
plot(t,x(:,2),t2,x2(:,2),t3,x3(:,2),t4,x4(:,2), t5,x5(:,2));
title('Floor 2 of Three Story Building (M1, M2 =10000 kg | K1, K2=100 N/m | C1, C2=10 kg/s)');
legend('floor2: f = nat_1 = ' + string(nat_freq_1) + ", A = 200", 'floor2: f = nat_1 = ' + string(nat_freq_1)+ ", A = 300",'floor2: f = nat_2 = ' + string(nat_freq_2)+ ", A = 200", 'floor2: f = nat_3 = ' + string(nat_freq_3)+ ", A = 200", 'floor2: f = nat_1 + nat_2 + nat_3 = ' + string(nat_freq_1 + nat_freq_2 + nat_freq_3) + ", A = 200");

figure;
plot(t,x(:,3),t2,x2(:,3),t3,x3(:,3),t4,x4(:,3), t5,x5(:,3));
title('Floor 3 of Three Story Building (M1, M2 =10000 kg | K1, K2=100 N/m | C1, C2=10 kg/s)');
legend('floor3: f = nat_1 = ' + string(nat_freq_1) + ", A = 200", 'floor3: f = nat_1 = ' + string(nat_freq_1)+ ", A = 300",'floor3: f = nat_2 = ' + string(nat_freq_2)+ ", A = 200", 'floor3: f = nat_3 = ' + string(nat_freq_3)+ ", A = 200", 'floor3: f = nat_1 + nat_2 + nat_3 = ' + string(nat_freq_1 + nat_freq_2 + nat_freq_3) + ", A = 200");

max_1 = max(x)
max_2 = max(x2)
max_3 = max(x3)
max_4 = max(x4)
max_5 = max(x5)
%}

%% Low Frequency Earthquake
%{
A = 20;
w = 0.5;

%Parameters
story1_1 = [10000 200 10]; %Mass (kg), Spring Constant (N/m), Damping (kg/s)
story1_2= [10000 200 10];
story1_3 = [10000 200 10];

story2_1 = [10000 200 30]; %Mass (kg), Spring Constant (N/m), Damping (kg/s)
story2_2= [10000 200 30];
story2_3 = [10000 200 30];

story3_1 = [1000000 200 10]; %Mass (kg), Spring Constant (N/m), Damping (kg/s)
story3_2= [1000000 200 10];
story3_3 = [1000000 200 10];
%End of parameters


tspan = [0 2000];
x1_0 = [0; 30; 0; 0];
[t,x] = ode45(@(t,x) odefcn2(t,x,story1_1, story1_2, A * cos(w * t)), tspan, x1_0);
[t4,x4] = ode45(@(t4,x4) odefcn2(t4,x4,story2_1, story2_2, A * cos(w * t4)), tspan, x1_0); 
[t5,x5] = ode45(@(t5,x5) odefcn2(t5,x5,story3_1, story3_2,A * cos(w * t5)), tspan, x1_0); 

x2_0 = [0; 30; 60; 0; 0; 0];
[t2,x2] = ode45(@(t2,x2) odefcn3(t2,x2,story1_1, story1_2, story1_3, A * cos(w * t2)), tspan, x2_0);
[t3,x3] = ode45(@(t3,x3) odefcn3(t3,x3,story2_1, story2_2, story2_3, A * cos(w * t3)), tspan, x2_0);
[t6,x6] = ode45(@(t6,x6) odefcn3(t6,x6,story3_1, story3_2, story3_3, A * cos(w * t6)), tspan, x2_0);

plot(t,x(:,1),"-",t,x(:,2),"-", t4,x4(:,1), 'x' ,t4,x4(:,2), 'x',  t5,x5(:,1), '*' ,t5,x5(:,2), '*');
title('Two Story Building: F = 20cos(0.5t)');
legend("Floor1: (10000 200 10)", "Floor2: (10000 200 10)", "Floor1: (10000 200 30)", "Floor2: (10000 200 30)", "Floor1: (1000000 200 10)", "Floor2: (1000000 200 10)");

figure;
plot(t2,x2(:,1),"-",t2,x2(:,2),"-", t2,x2(:,3),"-",t3,x3(:,1), 'x' ,t3,x3(:,2), 'x', t3,x3(:,3),'x' ,t6,x6(:,1), '*' ,t6,x6(:,2), '*', t6,x6(:,3), '*');
title('Three Story Building: F = 20cos(0.5t)');
legend("Floor1: (10000 200 10)", "Floor2: (10000 200 10)", "Floor3: (10000 200 10)", "Floor1: (10000 200 30)", "Floor2: (10000 200 30)", "Floor3: (10000 200 30)", "Floor1: (1000000 200 10)", "Floor2: (1000000 200 10)", "Floor3: (1000000 200 10)")
%}

%% High Frequency Earthquake
%{
A = 20;
w = 50;

%Parameters
story1_1 = [10000 200 10]; %Mass (kg), Spring Constant (N/m), Damping (kg/s)
story1_2= [10000 200 10];
story1_3 = [10000 200 10];

story2_1 = [10000 200 30]; %Mass (kg), Spring Constant (N/m), Damping (kg/s)
story2_2= [10000 200 30];
story2_3 = [10000 200 30];

story3_1 = [1000000 200 10]; %Mass (kg), Spring Constant (N/m), Damping (kg/s)
story3_2= [1000000 200 10];
story3_3 = [1000000 200 10];
%End of parameters

tspan = [0 2000];
x1_0 = [0; 30; 0; 0];
[t,x] = ode45(@(t,x) odefcn2(t,x,story1_1, story1_2, A * cos(w * t)), tspan, x1_0);
[t4,x4] = ode45(@(t4,x4) odefcn2(t4,x4,story2_1, story2_2, A * cos(w * t4)), tspan, x1_0); 
[t5,x5] = ode45(@(t5,x5) odefcn2(t5,x5,story3_1, story3_2,A * cos(w * t5)), tspan, x1_0); 

x2_0 = [0; 30; 60; 0; 0; 0];
[t2,x2] = ode45(@(t2,x2) odefcn3(t2,x2,story1_1, story1_2, story1_3, A * cos(w * t2)), tspan, x2_0);
[t3,x3] = ode45(@(t3,x3) odefcn3(t3,x3,story2_1, story2_2, story2_3, A * cos(w * t3)), tspan, x2_0);
[t6,x6] = ode45(@(t6,x6) odefcn3(t6,x6,story3_1, story3_2, story3_3, A * cos(w * t6)), tspan, x2_0);

plot(t,x(:,1),"-",t,x(:,2),"-", t4,x4(:,1), 'x' ,t4,x4(:,2), 'x',  t5,x5(:,1), '*' ,t5,x5(:,2), '*');
title('Two Story Building: F = 20cos(0.5t)');
legend("Floor1: (10000 200 10)", "Floor2: (10000 200 10)", "Floor1: (10000 200 30)", "Floor2: (10000 200 30)", "Floor1: (1000000 200 10)", "Floor2: (1000000 200 10)");

figure;
plot(t2,x2(:,1),"-",t2,x2(:,2),"-", t2,x2(:,3),"-",t3,x3(:,1), 'x' ,t3,x3(:,2), 'x', t3,x3(:,3),'x' ,t6,x6(:,1), '*' ,t6,x6(:,2), '*', t6,x6(:,3), '*');
title('Three Story Building: F = 20cos(0.5t)');
legend("Floor1: (10000 200 10)", "Floor2: (10000 200 10)", "Floor3: (10000 200 10)", "Floor1: (10000 200 30)", "Floor2: (10000 200 30)", "Floor3: (10000 200 30)", "Floor1: (1000000 200 10)", "Floor2: (1000000 200 10)", "Floor3: (1000000 200 10)")

%}

function dxdt = odefcn1(t,x,story_1, F1)
  m1 = story_1(1);
  k1 = story_1(2);
  c1 = story_1(3);
  
  dxdt = zeros(2,1);
  dxdt(1) = x(2);
  dxdt(2) = (1/m1) * (-c1*x(2) - k1 * x(1) + F1);
end

function dxdt = odefcn2(t,x,story_1, story_2, F1)

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

function dxdt = odefcn3(t,x,story_1, story_2, story_3, F1)
  m1 = story_1(1);
  k1 = story_1(2);
  c1 = story_1(3);
  
  m2 = story_2(1);
  k2 = story_2(2);
  c2 = story_2(3);
  
  m3 = story_3(1);
  k3 = story_3(2);
  c3 = story_3(3);
 
  dxdt = zeros(6,1);
  dxdt(1) = x(4);
  dxdt(2) = x(5);
  dxdt(3) = x(6);
  dxdt(4) = (1/m1) * (-k1 * x(1) - c1 * x(4) + k2 * (x(2) - x(1)) + c2 * (x(5) - x(4)) + F1);
  dxdt(5) = (1/m2) * (-k2 * (x(2) - x(1)) - c2 * (x(5) - x(4)) + k3 *(x(3) - x(2)) + c3 * (x(6) - x(5)));
  dxdt(6) = (1/m3) * (-k3 * (x(3) - x(2)) - c3 * (x(6) - x(5)));

end


