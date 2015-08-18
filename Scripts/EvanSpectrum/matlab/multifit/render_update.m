
	% Get the global dircosh, modcosh functions.
	global dircosh;
	global modcosh;

	% First just load up prepare_data. This does all the projecting, etc.
	% This guarantees rescale_* are all set.
	run prepare_data; 

	% Do purely visual rescaling. 
	rescale_vals = cosh(mass_answer*(parse_Nt/2 - xrange));
	rescale_y = (rescale_sum-vev_answer) ./ rescale_vals; % sub vev here
	rescale_yerr = rescale_err ./ rescale_vals;

	if (log_answer == 1)
		rescale_yt = log(rescale_y);
		rescale_yerrt = rescale_yerr./rescale_y;
		
		imag_flag = (abs(imag(rescale_yt)) > 1e-7);
		rescale_yt(imag_flag) = -1e100; % don't try to plot complex things.
		
		rescale_y = rescale_yt;
		rescale_yerr = rescale_yerrt;
		
		clear('rescale_yt'); clear('rescale_yerrt');
    end
	
	%clear('rescale_sum'); clear('rescale_err');
	
	% Plot!
	figure(h);
	errorbar(xrange, rescale_y, rescale_yerr, '.k'); hold on;
	[funcdata] = get_function(xrange, parse_Nt, fit_form, coefficients);
	
    if (log_answer == 1)
        
        logfuncdata = log((funcdata-vev_answer)./rescale_vals);
        imag_flag = (abs(imag(logfuncdata)) > 1e-7);
		logfuncdata(imag_flag) = -1e100; % don't try to plot complex things.
        
        plot(xrange, logfuncdata); hold off;
        
        
        tmp_flag = (rescale_y < -1e90);
        tmp_dat = rescale_y;
        tmp_dat(tmp_flag) = [];
		axis([-1 parse_Nt (min(tmp_dat)-log(1.1)) (max(tmp_dat)+log(1.1))]);
        clear('tmp_flag'); clear('tmp_dat');
    else
        
        plot(xrange, (funcdata-vev_answer)./rescale_vals); hold off;
        
		axis([-1 parse_Nt -max(abs(rescale_y))*1.1 max(abs(rescale_y))*1.1]);
	end

	
	% And clean up the info window.
	fold_string = '';
	if (func_fold == 0)
		fold_string = 'No';
	else
		fold_string = 'Yes';
	end
	
	ppp_string = '';
	if (fit_ppp == 0)
		ppp_string = 'No';
	elseif (fit_ppp == 1)
		ppp_string = 'Positive';
	else
		ppp_string = 'Negative';
	end
	
	zero_string = '';
	if (fit_zero == 0)
		zero_string = 'No';
	elseif (fit_zero == 1)
		zero_string = 'Yes';
	elseif (fit_zero == -1)
		zero_string = 'Gap';
	end
	
	
	
	even_string = '';
	if (fit_even == 0)
		even_string = 'All t';
	else
		even_string = 'Even t';
	end
	
	diag_string = '';
	if (fit_diag == 0)
		diag_string = 'Full Var-Covar';
	else
		diag_string = 'Diagonal Var-Covar';
    end
	
	cosh_string = '';
	if (is_baryon == 0 && is_sinh == 0)
		if (func_zshift == 0 && func_diff == 0)
			cosh_string = 'Normal';
		elseif (func_diff == 1) % zero shift doesn't matter.
			cosh_string = 'Finite Difference';
		else % must be zero shifted.
			cosh_string = 'Zero Shifted';
		end
	elseif (is_sinh == 1)
		if (func_zshift == 0 && func_diff == 0)
			cosh_string = 'Sinh';
		elseif (func_diff == 1) % zero shift doesn't matter.
			cosh_string = 'Sinh Finite Difference';
		else % must be zero shifted.
			cosh_string = 'Sinh Zero Shifted';
		end
	else % is baryon
		if (func_zshift == 0 && func_diff == 0)
			cosh_string = 'Baryon';
		elseif (func_diff == 1) % zero shift doesn't matter.
			cosh_string = 'Baryon Finite Difference';
		else % must be zero shifted.
			cosh_string = 'Baryon Zero Shifted';
		end
	end
	
	diff_string = '';
	if (fit_diff == 0)
		diff_string = 'No';
	else
		diff_string = 'By One';
	end
    
    saving_string = '';
    if (saving == 0)
        saving_string = 'Off';
    else
        saving_string = 'On';
    end
    
    blocking_string = ''; 
    if (is_blocking == 0)
        blocking_string = strcat(['--/', num2str(num_blocks)]);
    else
        blocking_string = strcat([num2str(b), '/', num2str(num_blocks)]);
    end
	
	subset_string = '';
	if (loc_begin == 0 && loc_end == 1)
		subset_string = 'All';
	else
		subset_string = strcat(['Fraction ', num2str(loc_begin), ' to ', num2str(loc_end)]);
	end
	
	% Build up string for the masses and whatnot in advance.
	
	vstr{1} = strcat(['Cosh: Ampl 1: ' num2str(coefficients(1))]);
	vstr{2} = strcat(['Cosh: Mass 1: ' num2str(coefficients(2))]);
	vstr{3} = strcat(['Cosh: Ampl 2: ' num2str(coefficients(3))]);
	vstr{4} = strcat(['Cosh: Mass 2: ' num2str(coefficients(4))]);
    vstr{5} = strcat(['Cosh: Ampl 3: ' num2str(coefficients(5))]);
	vstr{6} = strcat(['Cosh: Mass 3: ' num2str(coefficients(6))]);
	vstr{7} = strcat(['Oscl: Ampl 4: ' num2str(coefficients(7))]);
	vstr{8} = strcat(['Oscl: Mass 4: ' num2str(coefficients(8))]);
	vstr{9} = strcat(['Oscl: Ampl 5: ' num2str(coefficients(9))]);
	vstr{10} = strcat(['Oscl: Mass 5: ' num2str(coefficients(10))]);
	vstr{11} = strcat(['Oscl: Ampl 6: ' num2str(coefficients(11))]);
	vstr{12} = strcat(['Oscl: Mass 6: ' num2str(coefficients(12))]);
	
	
	% Add any constraints!
	for param_num = 1:12
		if (errors(param_num) ~= 0)
			vstr{param_num} = strcat([vstr{param_num} ' (' num2str(errors(param_num)) ')']);
		end
		if (constraints(param_num) ~= 0)
			vstr{param_num} = strcat([vstr{param_num} ' [Cnstr: ' num2str(constraints(param_num)) ']']);
		end
	end
	
	% Build up aux strings.
	if (is_fpi == 1)
		astr{1} = strcat(['F_pi: ' num2str(auxresult(1))]);
		astr{2} = strcat(['No auxiliary param 2.']);
		for param_num = 1:2
			if (auxerrors(param_num) ~= 0)
				astr{param_num} = strcat([astr{param_num} ' (' num2str(auxerrors(param_num)) ')']);
			end
		end
	else % no aux string.
		astr{1} = strcat(['No auxiliary param 1.']);
		astr{2} = strcat(['No auxiliary param 2.']);
	end
		
	set(th, 'string', {strcat(['Directory: ' directory_answer{1}]), ...
						strcat(['State: ' spectrum_list{spectrum_answer}]), ...
						'', ...
						strcat(['Rescale mass: ' num2str(mass_answer)]), ...
						strcat(['Rescale vev: ' num2str(vev_answer)]), ...
						'', ...
						strcat(['Fit Form: ' form_list{fit_form}]), ...
						strcat(['Cosh Type: ' cosh_string]), ...
						strcat(['Fit Range: [' num2str(fit_minimum) ':' num2str(fit_maximum) ']']), ...
						strcat(['Fit X Precision: ' num2str(fit_x_prec)]), ...
						strcat(['Fit Points: ' even_string]), ...
						strcat(['Cut Singular Values: ' num2str(fit_cut)]), ...
						strcat(['Variance-Covariance: ' diag_string]), ...
						'', ...
						strcat(['Fold: ' fold_string]), ...
						strcat(['Parity Project: ' ppp_string]), ...
						strcat(['Zero Center: ' zero_string]), ...
						strcat(['Finite Differences: ' diff_string]), ...
						'', ...
						strcat(['Subset: ' subset_string]), ...
						strcat(['Binsize: ' num2str(blocksize)]), ...
						strcat(['Eliminate on Jackknife: ' num2str(num_elim)]), ...
                        strcat(['Saving: ' saving_string]), ...
                        strcat(['Blocking status: ' blocking_string]), ...
                        '', ...
						vstr{1}, ...
						vstr{2}, ...
						vstr{3}, ...
						vstr{4}, ...
						vstr{5}, ...
						vstr{6}, ...
						vstr{7}, ...
						vstr{8}, ...
						vstr{9}, ...
						vstr{10}, ...
						vstr{11}, ...
						vstr{12}, ...
						astr{1}, ...
						astr{2}, ...
                        '', ...
						strcat(['Chisq DoF: ' num2str(chisq_dof)]), ...
						strcat(['P-value: ' num2str(p_val)]), ...
						strcat(['Condition Num: ' num2str(cond_num)])});

	clear('fold_string'); clear('ppp_string'); clear('zero_string'); clear('even_string'); clear('diag_string'); clear('vstr'); clear('saving_string'); clear('blocking_string'); clear('cosh_string'); clear('subset_string');
