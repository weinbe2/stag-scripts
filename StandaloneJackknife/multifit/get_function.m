function [yout] = get_function(xin, nt, type, coeff)
    
	% Get the global dircosh, modcosh functions.
	global dircosh;
	global modcosh;
	
    % coeff:
    % cosh: amp_1 mass_1
    % cosh: amp_2 mass_2
    % osc: amp_3 mass_3
    % osc: amp_4 mass_4

	switch type
         case 1 % cosh 1 oscil 0
            yout = +coeff(1)*dircosh(coeff(2),xin,nt);
         case 2 % cosh 2 oscil 0
            yout = +coeff(1)*dircosh(coeff(2),xin,nt)+coeff(3)*dircosh(coeff(4),xin,nt);
         case 3 % cosh 3 oscil 0
            yout = +coeff(1)*dircosh(coeff(2),xin,nt)+coeff(3)*dircosh(coeff(4),xin,nt)+coeff(5)*dircosh(coeff(6),xin,nt);
         case 4 % cosh 0 oscil 1
            yout = +coeff(7)*modcosh(coeff(8),xin,nt);
         case 5 % cosh 1 oscil 1
            yout = +coeff(1)*dircosh(coeff(2),xin,nt)+coeff(7)*modcosh(coeff(8),xin,nt);
         case 6 % cosh 2 oscil 1
            yout = +coeff(1)*dircosh(coeff(2),xin,nt)+coeff(3)*dircosh(coeff(4),xin,nt)+coeff(7)*modcosh(coeff(8),xin,nt);
         case 7 % cosh 3 oscil 1
            yout = +coeff(1)*dircosh(coeff(2),xin,nt)+coeff(3)*dircosh(coeff(4),xin,nt)+coeff(5)*dircosh(coeff(6),xin,nt)+coeff(7)*modcosh(coeff(8),xin,nt);
         case 8 % cosh 0 oscil 2
            yout = +coeff(7)*modcosh(coeff(8),xin,nt)+coeff(9)*modcosh(coeff(10),xin,nt);
         case 9 % cosh 1 oscil 2
            yout = +coeff(1)*dircosh(coeff(2),xin,nt)+coeff(7)*modcosh(coeff(8),xin,nt)+coeff(9)*modcosh(coeff(10),xin,nt);
         case 10 % cosh 2 oscil 2
            yout = +coeff(1)*dircosh(coeff(2),xin,nt)+coeff(3)*dircosh(coeff(4),xin,nt)+coeff(7)*modcosh(coeff(8),xin,nt)+coeff(9)*modcosh(coeff(10),xin,nt);
         case 11 % cosh 3 oscil 2
            yout = +coeff(1)*dircosh(coeff(2),xin,nt)+coeff(3)*dircosh(coeff(4),xin,nt)+coeff(5)*dircosh(coeff(6),xin,nt)+coeff(7)*modcosh(coeff(8),xin,nt)+coeff(9)*modcosh(coeff(10),xin,nt);
         case 12 % cosh 0 oscil 3
            yout = +coeff(7)*modcosh(coeff(8),xin,nt)+coeff(9)*modcosh(coeff(10),xin,nt)+coeff(11)*modcosh(coeff(12),xin,nt);
         case 13 % cosh 1 oscil 3
            yout = +coeff(1)*dircosh(coeff(2),xin,nt)+coeff(7)*modcosh(coeff(8),xin,nt)+coeff(9)*modcosh(coeff(10),xin,nt)+coeff(11)*modcosh(coeff(12),xin,nt);
         case 14 % cosh 2 oscil 3
            yout = +coeff(1)*dircosh(coeff(2),xin,nt)+coeff(3)*dircosh(coeff(4),xin,nt)+coeff(7)*modcosh(coeff(8),xin,nt)+coeff(9)*modcosh(coeff(10),xin,nt)+coeff(11)*modcosh(coeff(12),xin,nt);
         case 15 % cosh 3 oscil 3
            yout = +coeff(1)*dircosh(coeff(2),xin,nt)+coeff(3)*dircosh(coeff(4),xin,nt)+coeff(5)*dircosh(coeff(6),xin,nt)+coeff(7)*modcosh(coeff(8),xin,nt)+coeff(9)*modcosh(coeff(10),xin,nt)+coeff(11)*modcosh(coeff(12),xin,nt);

	end
	
end