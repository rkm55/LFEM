clear;clc;

% Homework question 53-55

% Solve with everything in terms of xsi and across -1 to 1
syms x y
degree = 3; %user defined
F_Leg = zeros(1,degree+1);
M_Leg = zeros(degree+1);
F_Ber = zeros(1,degree+1);
M_Ber = zeros(degree+1);
F_Lag = zeros(1,degree+1);
M_Lag = zeros(degree+1);

x = (1/2)*(y + 1);
f_x = x^3 - (8*x^2)/5 + (3*x)/5;
nodes = linspace(-1,1,degree+1);

N_Leg(1) = x^0;
N_Leg(2) = 2*x - 1;
N_Leg(3) = (1/2)*(3*(2*x - 1)^2 - 1);
N_Leg(4) = (1/6)*(15*(2*x-1)^3 - 9*(2*x-1));
for basis_idx = 0:degree
    N_Ber(basis_idx+1) = nchoosek(degree,basis_idx) * (x)^basis_idx * (1-(x))^(degree-basis_idx);
    val = x^0;
    for j = 0:degree
        if basis_idx == j
            N_Lag(basis_idx+1) = val;
        else
            val = val*((2*x-1) - nodes(j+1)) / (nodes(basis_idx+1) - nodes(j+1));
            N_Lag(basis_idx+1) = val;
        end
    end
end
for i = 1:(degree+1)
    for j = 1:(degree+1)
        M_Leg(i,j) = double(int(N_Leg(i)*N_Leg(j),-1,1)/2);
        M_Ber(i,j) = double(int(N_Ber(i)*N_Ber(j),-1,1)/2);
        M_Lag(i,j) = double(int(N_Lag(i)*N_Lag(j),-1,1)/2);
    end
    F_Leg(i) = double(int(f_x*N_Leg(i),-1,1)/2);
    F_Ber(i) = double(int(f_x*N_Ber(i),-1,1)/2);
    F_Lag(i) = double(int(f_x*N_Lag(i),-1,1)/2);
end
d_Leg = inv(M_Leg)*F_Leg';
d_Ber = inv(M_Ber)*F_Ber';
d_Lag = inv(M_Lag)*F_Lag';

%% part 2 q55
x = linspace(-1,1,4);
f_x = x.^3 - (8.*x.^2)/5 + (3.*x)./5;

xx = linspace(-1,1,100);
f_xx = xx.^3 - (8.*xx.^2)/5 + (3.*xx)./5;

% 3 points
f_x3 = (-8/5).*xx.^2 + (8/5).*xx;
f_x4 = (26930000000/26934006633).*xx.^3 - (142225000/88891111).*xx.^2 + (80822053064/134670033165).*xx - (3888/444455555);

plot(xx,f_xx,'--b',LineWidth=2)
hold on 
plot(xx,f_x3,LineWidth=1.2)
plot(xx,f_x4,LineWidth=1.2)
legend('f(x)','3 points','4 points','Location','southeast')
title('Lagrange Interpolatory Polynomials Plotted')





