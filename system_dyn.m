% ode solver for the system in ECE HW 1

function dM = system_dyn(~,M,theta)

dM = zeros(size(M));
dM(1) = 1;
dM(2) = theta(2)*M(3);
dM(3) = -1*theta(1)*M(1)*M(2);

end