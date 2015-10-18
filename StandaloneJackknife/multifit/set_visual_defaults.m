% Set defaults for various settings.


	% Options: tmin, tmax
	% Fold? y/n
	% Positive Parity Project? y/n
	% Fit to only even time data? y/n
	% Full or diagonal correlator matrix? 

	fit_minimum = 10;
	fit_maximum = 24;
	fit_ppp = 0; % don't positive parity project
	fit_zero = 0; % don't normalize center to zero or "gap" it.
	fit_even = 0; % fit all
	fit_diag = 0; % fit full correlation matrix instead of diagonal.
	fit_diff = 0; % fit original data, not finite diffs. 
	fit_x_prec = 1e-20; % set the x precision of the fit.
	
	fit_cut = 0; % pick how many singular values to cut.
				 % Note that this also modifies the dof.
	
	% Modify fit function: 
	func_fold = 1; % fold 
	func_zshift = 0; % don't use zero-shifted cosh.
	func_diff = 0; % use original cosh.