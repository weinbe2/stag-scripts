function [yout] = get_function(xin, nt, type, coeff)
    
    % coeff:
    % cosh: amp_1 mass_1
    % cosh: amp_2 mass_2
    % osc: amp_3 mass_3
    % osc: amp_4 mass_4

	switch type
         case 1 % cosh 1 oscil 0
            yout = +coeff(1)*cosh(coeff(2)*(nt/2-xin));
         case 2 % cosh 2 oscil 0
            yout = +coeff(1)*cosh(coeff(2)*(nt/2-xin))+coeff(3)*cosh(coeff(4)*(nt/2-xin));
         case 3 % cosh 3 oscil 0
            yout = +coeff(1)*cosh(coeff(2)*(nt/2-xin))+coeff(3)*cosh(coeff(4)*(nt/2-xin))+coeff(5)*cosh(coeff(6)*(nt/2-xin));
         case 4 % cosh 0 oscil 1
            yout = +coeff(7)*(1-2*mod(xin,2)).*cosh(coeff(8)*(nt/2-xin));
         case 5 % cosh 1 oscil 1
            yout = +coeff(1)*cosh(coeff(2)*(nt/2-xin))+coeff(7)*(1-2*mod(xin,2)).*cosh(coeff(8)*(nt/2-xin));
         case 6 % cosh 2 oscil 1
            yout = +coeff(1)*cosh(coeff(2)*(nt/2-xin))+coeff(3)*cosh(coeff(4)*(nt/2-xin))+coeff(7)*(1-2*mod(xin,2)).*cosh(coeff(8)*(nt/2-xin));
         case 7 % cosh 3 oscil 1
            yout = +coeff(1)*cosh(coeff(2)*(nt/2-xin))+coeff(3)*cosh(coeff(4)*(nt/2-xin))+coeff(5)*cosh(coeff(6)*(nt/2-xin))+coeff(7)*(1-2*mod(xin,2)).*cosh(coeff(8)*(nt/2-xin));
         case 8 % cosh 0 oscil 2
            yout = +coeff(7)*(1-2*mod(xin,2)).*cosh(coeff(8)*(nt/2-xin))+coeff(9)*(1-2*mod(xin,2)).*cosh(coeff(10)*(nt/2-xin));
         case 9 % cosh 1 oscil 2
            yout = +coeff(1)*cosh(coeff(2)*(nt/2-xin))+coeff(7)*(1-2*mod(xin,2)).*cosh(coeff(8)*(nt/2-xin))+coeff(9)*(1-2*mod(xin,2)).*cosh(coeff(10)*(nt/2-xin));
         case 10 % cosh 2 oscil 2
            yout = +coeff(1)*cosh(coeff(2)*(nt/2-xin))+coeff(3)*cosh(coeff(4)*(nt/2-xin))+coeff(7)*(1-2*mod(xin,2)).*cosh(coeff(8)*(nt/2-xin))+coeff(9)*(1-2*mod(xin,2)).*cosh(coeff(10)*(nt/2-xin));
         case 11 % cosh 3 oscil 2
            yout = +coeff(1)*cosh(coeff(2)*(nt/2-xin))+coeff(3)*cosh(coeff(4)*(nt/2-xin))+coeff(5)*cosh(coeff(6)*(nt/2-xin))+coeff(7)*(1-2*mod(xin,2)).*cosh(coeff(8)*(nt/2-xin))+coeff(9)*(1-2*mod(xin,2)).*cosh(coeff(10)*(nt/2-xin));
         case 12 % cosh 0 oscil 3
            yout = +coeff(7)*(1-2*mod(xin,2)).*cosh(coeff(8)*(nt/2-xin))+coeff(9)*(1-2*mod(xin,2)).*cosh(coeff(10)*(nt/2-xin))+coeff(11)*(1-2*mod(xin,2)).*cosh(coeff(12)*(nt/2-xin));
         case 13 % cosh 1 oscil 3
            yout = +coeff(1)*cosh(coeff(2)*(nt/2-xin))+coeff(7)*(1-2*mod(xin,2)).*cosh(coeff(8)*(nt/2-xin))+coeff(9)*(1-2*mod(xin,2)).*cosh(coeff(10)*(nt/2-xin))+coeff(11)*(1-2*mod(xin,2)).*cosh(coeff(12)*(nt/2-xin));
         case 14 % cosh 2 oscil 3
            yout = +coeff(1)*cosh(coeff(2)*(nt/2-xin))+coeff(3)*cosh(coeff(4)*(nt/2-xin))+coeff(7)*(1-2*mod(xin,2)).*cosh(coeff(8)*(nt/2-xin))+coeff(9)*(1-2*mod(xin,2)).*cosh(coeff(10)*(nt/2-xin))+coeff(11)*(1-2*mod(xin,2)).*cosh(coeff(12)*(nt/2-xin));
         case 15 % cosh 3 oscil 3
            yout = +coeff(1)*cosh(coeff(2)*(nt/2-xin))+coeff(3)*cosh(coeff(4)*(nt/2-xin))+coeff(5)*cosh(coeff(6)*(nt/2-xin))+coeff(7)*(1-2*mod(xin,2)).*cosh(coeff(8)*(nt/2-xin))+coeff(9)*(1-2*mod(xin,2)).*cosh(coeff(10)*(nt/2-xin))+coeff(11)*(1-2*mod(xin,2)).*cosh(coeff(12)*(nt/2-xin));

	end
	
end