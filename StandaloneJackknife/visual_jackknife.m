function visual_jackknife(ensemble, state, blockval)
	% Load some useful directories!
	if (~isdeployed)
		addpath('.\process', '-end');
		addpath('.\multifit', '-end');
	end	%addpath('/projectnb/qcd/Staggered8f/T2015-09-03CompileMatlab/process', '-end');
	%addpath('/projectnb/qcd/Staggered8f/T2015-09-03CompileMatlab/multifit', '-end');
	rel_path = '/projectnb/qcd/Staggered/';
    %rel_path = '..\2014-12-09VisualMultiExp\';

	if (~exist(strcat([rel_path, ensemble]),'dir'))
		disp('Ensemble does not exist.')
		exit;
	end
	
	[fl_l, fl_h, parse_Ns, parse_Nt, beta, m_l, m_h] = load_info_file(strcat([rel_path ensemble]));
	
	if (~exist(strcat([rel_path, ensemble, '/spectrum2/multifits/multifits.', state]),'file') || ~exist(strcat([rel_path, ensemble, '/spectrum2/corr/corr.', state]),'file'))
		disp('Multifit and corr does not exist.')
	end
	
	disp('Ready to load!')
	
	mass_l = m_l;
    blocksize = blockval;
    num_elim = 1;

	number_files = 1;
	volume = parse_Ns^3; %*48;
	fl_flavor = fl_l/4;
    
    % set what % of the data we look at.
    loc_begin = 0.0;
    loc_end = 1.0;
    
	% Import data.
    
    [scalar_sum, scalar_jack, scalar_cov_mat, ...
		scalar_err, num_blocks, scalar_jack_single] = ...
		get_correlator(strcat([rel_path, ensemble]), ...
		state, parse_Nt, parse_Ns, fl_flavor, blocksize, num_elim, ...
        loc_begin, loc_end);
	
	% Get rolling!
	xrange = (0:(parse_Nt-1))';
	
	rescale_y = scalar_sum;
	rescale_yerr = scalar_err;
	mass_answer = 0; % Effectively rescaling by a zero mass, which gives 1.
	vev_answer = 0; % subtract off this vev
	log_answer = 0; % plot on a log scale. 
	fit_form = 1; % One cosh. 
    


	%run set_visual_defaults;
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
	
	% See if it's a baryon...
	if (strcmp(state, 'nu') || strcmp(state, 'de'))
		is_baryon = 1;
	else
		is_baryon = 0;
	end
	
	% See if it's a negparity state.
	if (strcmp(state, 'plc'))
		is_sinh = 1;
	else
		is_sinh = 0;
	end
	
	% See if we want f_pi, or another f.
	if (strcmp(state, 'ps2') || strcmp(state, 'pll') || strcmp(state, 'plc') || strcmp(state, 'rcc'))
		is_fpi = 1;
	else
		is_fpi = 0;
	end
	
	% If it's 'dc_stoch' or 'sg_stoch', finite difference it!
	if (strcmp(state, 'dc_stoch') || strcmp(state, 'sg_stoch'))
		func_diff = 1;
		fit_diff = 1;
	end
	
	% Get correct cosh functions in place.
	%run get_cosh; 
	
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
	
	% Prepare the data.
	%run prepare_data;
	% Prepare data for a fit depending on flags.

% Throw the central value, jackknife blocks, and single
% elim jackknife blocks together.
rescale_corr = cat(2, scalar_sum, scalar_jack, scalar_jack_single);

% Recall the single elim jackknife blocks are used to define the
% var-covar matrix. 

% First---check if we fold!
if (func_fold == 1) % fold!
	% Use the fold_data function!
	rescale_corr = fold_data(rescale_corr, is_baryon, is_sinh);
end

% Next---check if we positive parity project it!
if (fit_ppp == 1) % positive parity project!
	% Take care of first data.
	tmp_first = 0.25*(rescale_corr(parse_Nt,:)+2*rescale_corr(1,:)+rescale_corr(2,:));
    tmp_meh = rescale_corr(1,:); % last term gets special treatment
    for tmp_i = 0:(parse_Nt-3)
        rescale_corr(tmp_i+1,:) = rescale_corr(tmp_i+1,:)+2*rescale_corr(tmp_i+2,:)+rescale_corr(tmp_i+3,:);
    end
    rescale_corr(parse_Nt-1,:) = rescale_corr(parse_Nt-1,:)+2*rescale_corr(parse_Nt,:)+tmp_meh;
    clear('tmp_meh');
    % Shift things forward.
    for tmp_i = (parse_Nt):-1:2
        rescale_corr(tmp_i,:) = 0.25*rescale_corr(tmp_i-1,:);
    end
    rescale_corr(1,:) = tmp_first; % we lose the very edge.
elseif (fit_ppp == -1) % negative parity project!
	% Take care of first data. Minus is to be consistent.
	tmp_first = -0.25*(-rescale_corr(parse_Nt,:)+2*rescale_corr(1,:)-rescale_corr(2,:));
    tmp_meh = rescale_corr(1,:); % last term gets special treatment
    for tmp_i = 0:(parse_Nt-3)
        rescale_corr(tmp_i+1,:) = (1-2*mod(tmp_i,2))*(-rescale_corr(tmp_i+1,:)+2*rescale_corr(tmp_i+2,:)-rescale_corr(tmp_i+3,:));
    end
    rescale_corr(parse_Nt-1,:) = rescale_corr(parse_Nt-1,:)-2*rescale_corr(parse_Nt,:)+tmp_meh;
    clear('tmp_meh');
    % Shift things forward.
    for tmp_i = (parse_Nt):-1:2
        rescale_corr(tmp_i,:) = 0.25*rescale_corr(tmp_i-1,:);
    end
    rescale_corr(1,:) = tmp_first; % restore end.
end

% Next---check if we zero shift it!
% Zero shifting just subtracts the value at the
% center of the correlator from the mean and
% from the jackknife blocks, so the value and error
% at the center is exactly zero.

% Gap: subtract the value at the center of the mean
% from the mean and all jackknife blocks, so there's
% still errors on the center (errors are, in fact,
% unchanged. If there's a constant in
% the fit, it'll exactly absorb that.
if (fit_zero == 1) % shift to zero
	shift_middle = rescale_corr((parse_Nt/2+1),:); 
    rep_shift_middle = repmat(shift_middle, [parse_Nt 1]);
    %shift_middle = rescale_corr((parse_Nt/2+1),1); 
    %rep_shift_middle = repmat(shift_middle, [parse_Nt size(rescale_corr, 2)]);
    rescale_corr = rescale_corr - rep_shift_middle;
    clear('shift_middle'); clear('rep_shift_middle');
elseif (fit_zero == -1) % 'gap' it.
	shift_middle = rescale_corr((parse_Nt/2+1),1); 
    rep_shift_middle = repmat(shift_middle, [parse_Nt size(rescale_corr, 2)]);
    rescale_corr = rescale_corr - rep_shift_middle;
    clear('shift_middle'); clear('rep_shift_middle');
end

% Finite difference the data!
if (fit_diff == 1) % finite differences!
	% First, compute finite difference of last data.
	tmp_last = rescale_corr(parse_Nt,:)-rescale_corr(1,:);
	for tmp_i=1:parse_Nt-1
		rescale_corr(tmp_i,:) = rescale_corr(tmp_i,:)-rescale_corr(tmp_i+1,:);
	end
	rescale_corr(parse_Nt,:) = tmp_last;
end

% Once all these things are done, get the error.

rescale_sum = rescale_corr(:,1);
rescale_jack = rescale_corr(:,2:(size(scalar_jack,2)+1));
rescale_jack_single = rescale_corr(:,(size(scalar_jack,2)+2):end);
[rescale_cov_mat, rescale_err] = ...
	errors_jackknife(rescale_sum, rescale_jack_single);

clear('rescale_jack_single');




	% Load multifit file.
	results_save = importdata(strcat([rel_path, ensemble, '/spectrum2/multifits/multifits.', state]));
	
	for s=2:size(results_save,1)
	
		% get the relevant parameters and initial guesses.
		fit_minimum = results_save(s, 1);
		fit_maximum = results_save(s, 2);
		fit_form = results_save(s, 3);
		coefficients = results_save(s, 4:2:26)';
		errors = results_save(s, 5:2:27);
		chisq_dof = results_save(s, 28);
		p_val = results_save(s, 29);
		constraints = results_save(s, 30:41);
		auxresult = results_save(s, 42:2:44);
		auxerrors = results_save(s, 43:2:45);
		
        disp(strcat(['Tmin ', num2str(fit_minimum)]))
        
		% If errors already exist, skip it!
		if (sum(abs(errors)) > 0)
			continue
		end

		% Whelp, here we go!
		the_fit_output = get_all_nlfit_multi(rescale_sum, rescale_cov_mat, fit_minimum, fit_maximum, parse_Nt, fit_even, fit_cut, func_diff, mod(fit_form,4), (fit_form-mod(fit_form,4))/4, fit_diag, fit_x_prec, coefficients, constraints);
		
		
		if (~(size(the_fit_output, 1) == 0))
		
			%coefficients(:) = the_fit_output(1, 3:14);
			
			chisq_dof = the_fit_output(1, 15);
			p_val = the_fit_output(1, 16);
			cond_num = the_fit_output(1, 17);
			
			% Into the rabbit hole...
			coefficients_blocks = zeros(num_blocks, 12);
			auxresult_blocks = zeros(num_blocks, 2);
			
			jack_flag = 1; is_blocking = 1;
			for b=1:num_blocks
				if (jack_flag == 1)
					block_fit_output = get_all_nlfit_multi(rescale_jack(:,b), rescale_cov_mat, fit_minimum, fit_maximum, parse_Nt, fit_even, fit_cut, func_diff, mod(fit_form,4), (fit_form-mod(fit_form,4))/4, fit_diag, fit_x_prec, coefficients, constraints);
				
					if (numel(block_fit_output) == 0)
						jack_flag = 0;
						continue;
					end
				
					coefficients_blocks(b, :) = block_fit_output(3:14);
					
					% Check auxillary results, such as F_pi.
					if (is_fpi == 1)
						if (strcmp(state, 'ps2'))
							auxresult_blocks(b,1) = 2*m_l*sqrt(coefficients_blocks(b,1)*cosh(parse_Nt*coefficients_blocks(b,2)/2)*((parse_Ns)^3)*3/(4*(coefficients_blocks(b,2)^3))); % f_pi
						elseif (strcmp(state, 'pll'))
							auxresult_blocks(b,1) = 2*m_l*sqrt(coefficients_blocks(b,1)*cosh(parse_Nt*coefficients_blocks(b,2)/2)*((parse_Ns)^3)/(4*(coefficients_blocks(b,2)^3))); % f_pi
						elseif (strcmp(state, 'plc'))
							auxresult_blocks(b,1) = sqrt(coefficients_blocks(b,1)*((parse_Ns)^3)*m_l*sinh(coefficients_blocks(b,2)*parse_Nt/2)/(4*(coefficients_blocks(b,2)^2)*cosh(coefficients_blocks(b,2)/2))); % f_pi
						elseif (strcmp(state, 'rcc')) % F_v
							auxresult_blocks(b,1) = sqrt(coefficients_blocks(b,1)*(parse_Ns^3)*cosh(coefficients_blocks(b,2)*parse_Nt/2)/(16*coefficients_blocks(b,2)));
						end
					end
				end
			end
			if jack_flag == 1 % it worked!
				coefficients_rep = repmat(coefficients', [num_blocks 1]);
				errors = sqrt(sum((coefficients_blocks - coefficients_rep).^2,1).*(num_blocks-1)./num_blocks);
				
				auxresult_rep = repmat(auxresult, [num_blocks 1]);
				auxerrors = sqrt(sum((auxresult_blocks - auxresult_rep).^2,1).*(num_blocks-1)./num_blocks);
				
			else % it failed
				errors = zeros(1,12);
				auxerrors = zeros(1,2);
			end
			
			is_blocking  = 0;
			
			% Push errors to the save stack.								
			results_save(s, 4:2:26) = coefficients(:);
			results_save(s, 5:2:27) = errors(:);
			results_save(s, 28) = chisq_dof;
			results_save(s, 29) = p_val;
			results_save(s, 42:2:44) = auxresult(:);
			results_save(s, 43:2:45) = auxerrors(:);
		
		else
			errors = zeros(12, 1);
		end
	
	
	end
	
    save(strcat([rel_path, ensemble, '/spectrum2/multifits/multifits.', state]), 'results_save', '-ascii', '-double');
end