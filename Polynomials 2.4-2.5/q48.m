% hw ME 507 polynomials 2.4 - 2.5
% problem 48
% linear moment fitting system to determine weights

% degree 2
I = [2; 0];
F = [1 1; -1/sqrt(3) 1/sqrt(3)];
w = F\I;
disp(w)

% degree 3
I = [2; 0; 0];
F = [1 1 1; -sqrt(3/5) 0 sqrt(3/5); 2/5 -1/2 2/5];
w = F\I;
disp(w)