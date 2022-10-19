clc; close all; clear;
% Questions 65 and 66
% Splines 3.3 Hw

%% Question 65 
% C2 continuity means R matrices have 3 rows

RI1 = [0 0 0 1 -1 0 0 0 0 0 0 0;...
    0 0 -3 3 3 -3 0 0 0 0 0 0;...
    0 6 -12 6 -6 12 -6 0 0 0 0 0];

RI2 = [0 0 0 0 0 0 0 1 -1 0 0 0;...
    0 0 0 0 0 0 -3 3 3 -3 0 0;...
    0 0 0 0 0 6 -12 6 -6 12 -6 0];
RB = [RI1; RI2];
disp(' ')
disp('RB = ')
disp(' ')
disp(RB)

A = abs(RI1);
B = abs(RI2);
D = A.*B;   % 
disp(' D =   *non-zero entries in D mean the constraint systems interact')
disp(' ')
disp(D)


%% Question 66 

Rf = RB(:, 4:9);
if diff(size(Rf)) == 0
    disp(' ')
    disp('Rf = RB{3,4,5,6,7,8}     *if first col is idx 0')
    disp('     *Rf is a square matrix with rank 0 nullspace')
    disp(' ')
    disp(Rf)
else
    disp('Rf is not square')
end