
function f = ecdif(t,theta)
% resuelve sistemas de ecuaciones clasicas DM

L = length(theta); D = 1;
f(1) = 0.5*D*cos(theta(2));
f(2:L-1) = 0.5*D*(cos(theta(3:L))-cos(theta(1:L-2)));
f(L) = -0.5*D*cos(theta(L-1));
f = f';