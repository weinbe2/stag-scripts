% Depending on the fit settings, prepare a cosh and oscillating function
% that has the appropriate shifts or finite differences. 
% These get stored in the anonymous functions 'dircosh' and 'modcosh'
% for "direct cosh" and "modulated cosh," respectively.

global dircosh;
global modcosh;

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