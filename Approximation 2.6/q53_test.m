clear;clc;

% Homework question 54

% Solve with everything in terms of x and across 0 to 1
syms x y
degree = 2;
F2 = zeros(1,degree+1);
M2 = zeros(degree+1);
f_x = x^3 - (8*x^2)/5 + (3*x)/5;
N(2) = 2*x - 1;
N(3) = (1/2)*(3*(2*x - 1)^2 - 1);
N(1) = 1;
for i = 1:length(N)
    F2(i) = double(int(f_x*N(i),0,1));
end
for i = 1:3
    for j = 1:3
        M2(i,j) = double(int(N(i)*N(j),0,1));
    end
end
%%
% Solve with everything in terms of x and across 0 to 1
syms x y
degree = 3;
F = zeros(1,degree+1);
M = zeros(degree+1);
f_x = x^3 - (8*x^2)/5 + (3*x)/5;
N(2) = 2*x - 1;
N(3) = (1/2)*(3*(2*x - 1)^2 - 1);
%N(4) = (1/2)*()
N(1) = 1;
for i = 1:length(N)
    F(i) = double(int(f_x*N(i),0,1));
end
for i = 1:degree+1
    for j = 1:degree+1
        M(i,j) = double(int(N(i)*N(j),0,1))
    end
end



% % %% Solve with everything in terms of y and across -1 to 1
% % x = (1/2)*(y + 1);
% % f_x = x^3 - (8*x^2)/5 + (3*x)/5;
% % N0 = 1;
% % N1 = 2*x - 1;
% % N2 = (1/2)*(3*(2*x - 1)^2 - 1);
% % F0y = double(int(f_x*N0,-1,1)/2);
% % F1y = double(int(f_x*N1,-1,1)/2);
% % F2y = double(int(f_x*N2,-1,1)/2);