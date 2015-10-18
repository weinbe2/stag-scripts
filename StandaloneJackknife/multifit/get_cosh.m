% Depending on the fit settings, prepare a cosh and oscillating function
% that has the appropriate shifts or finite differences. 
% These get stored in the anonymous functions 'dircosh' and 'modcosh'
% for "direct cosh" and "modulated cosh," respectively.

global dircosh;
global modcosh;


if (is_baryon == 0 && is_sinh == 0)
	if (func_diff == 1) % Is it a finite diff?
		dircosh = @(m,t,nt)(cosh(m*(nt/2-t))-cosh(m*(nt/2-(t+1))));
		modcosh = @(m,t,nt)((1-2*mod(t,2)).*(cosh(m*(nt/2-t)))-(1-2*mod(t+1,2)).*(cosh(m*(nt/2-(t+1)))));
	elseif (func_zshift == 1) % Is it a zeroed cosh?
		dircosh = @(m,t,nt)(cosh(m*(nt/2-t))-1);
		modcosh = @(m,t,nt)((1-2*mod(t,2)).*(cosh(m*(nt/2-t))-1));
	else
		% Regular fit functions.
		dircosh = @(m,t,nt)(cosh(m*(nt/2-t)));
		modcosh = @(m,t,nt)((1-2*mod(t,2)).*cosh(m*(nt/2-t)));
	end
elseif (is_sinh == 1) % one operator is T negative
	if (func_diff == 1) % Is it a finite diff?
		dircosh = @(m,t,nt)(sinh(m*(nt/2-t))-sinh(m*(nt/2-(t+1))));
		modcosh = @(m,t,nt)((1-2*mod(t,2)).*(sinh(m*(nt/2-t)))-(1-2*mod(t+1,2)).*(sinh(m*(nt/2-(t+1)))));
	elseif (func_zshift == 1) % Is it a zeroed sinh? So a sinh...
		dircosh = @(m,t,nt)(sinh(m*(nt/2-t)));
		modcosh = @(m,t,nt)((1-2*mod(t,2)).*sinh(m*(nt/2-t)));
	else
		% Regular fit functions.
		dircosh = @(m,t,nt)(sinh(m*(nt/2-t)));
		modcosh = @(m,t,nt)((1-2*mod(t,2)).*sinh(m*(nt/2-t)));
	end
else % it's a baryon!
	% Baryon cosh function!
	bcosh = @(x,t)((exp(x)-(1-2*mod(t,2)).*exp(-x))/2.0);
	
	if (func_diff == 1) % Is it a finite diff?
		dircosh = @(m,t,nt)(bcosh(m*(nt/2-t),t)-bcosh(m*(nt/2-(t+1)),t));
		modcosh = @(m,t,nt)((1-2*mod(t,2)).*(bcosh(m*(nt/2-t),t))-(1-2*mod(t+1,2)).*(bcosh(m*(nt/2-(t+1)),t)));
	elseif (func_zshift == 1) % Is it a zeroed cosh?
		dircosh = @(m,t,nt)(bcosh(m*(nt/2-t),t)-bcosh(0,t));
		modcosh = @(m,t,nt)((1-2*mod(t,2)).*(bcosh(m*(nt/2-t),t)-bcosh(0,t)));
	else
		% Regular fit functions.
		dircosh = @(m,t,nt)(bcosh(m*(nt/2-t),t));
		modcosh = @(m,t,nt)((1-2*mod(t,2)).*bcosh(m*(nt/2-t),t));
	end
	
end