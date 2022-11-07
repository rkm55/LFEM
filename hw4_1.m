%% 16
clc; close all; clear;

a = [1; 4; 2];
b = [3; 6; 0];
i = [1; 0; 0];
j = [0; 1; 0];
k = [0; 0; 1];

C = crpr(a,k);

%% 20
clc; close all; clear;

E = 100000; %MPa
v = 0.25;
u = E/(2*(1+v));

sigma = [100; 0; 0; -50; 0; 0;];
S = [1/E -v/E -v/E 0 0 0;...
    -v/E 1/E -v/E 0 0 0;...
    -v/E -v/E 1/E 0 0 0;...
    0 0 0 1/u 0 0;...
    0 0 0 0 1/u 0;...
    0 0 0 0 0 1/u];

epsilon = S*sigma;

% 21
% S = compliance matrix, C = stiffness Matrix
C = inv(S)/1000;


%% Function(s)

function C = crpr(a,b)
    C = zeros(3,1);
    C(1) = a(2)*b(3) - a(3)*a(2);
    C(2) = a(3)*b(1) - a(1)*a(3);
    C(3) = a(1)*b(2) - a(2)*a(1);
end

