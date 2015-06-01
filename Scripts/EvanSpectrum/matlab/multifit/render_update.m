
	% First just load up prepare_data. This does all the projecting, etc.
	% This guarantees rescale_* are all set.
	run prepare_data; 

	% Do purely visual rescaling. 
	rescale_vals = cosh(mass_answer*(parse_Nt/2 - xrange));
	rescale_y = (rescale_sum-vev_answer) ./ rescale_vals; % sub vev here
	rescale_yerr = rescale_err ./ rescale_vals;

	if (log_answer == 1)
		rescale_yt = log(rescale_y);
		rescale_yerrt = rescale_err/rescale_y;
		
		imag_flag = (abs(imag(rescale_yt)) > 1e-7);
		rescale_yt(flag) = -1; % don't try to plot complex things.
		
		rescale_y = rescale_yt;
		rescale_yerr = rescale_yerrt;
		
		clear('rescale_yt'); clear('rescale_yerr');
	end
	
	clear('rescale_sum'); clear('rescale_err');
	
	% Plot!
	figure(h);
	errorbar(xrange, rescale_y, rescale_yerr, '.k'); hold on;
	[funcdata] = get_function(xrange, parse_Nt, fit_form, coefficients);
	plot(xrange, (funcdata-vev_answer)./rescale_vals); hold off;
	axis([-1 parse_Nt -max(abs(rescale_y))*1.1 max(abs(rescale_y))*1.1]);

	
	% And clean up the info window.
	fold_string = '';
	if (fit_fold == 0)
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
	else
		zero_string = 'Yes';
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

	set(th, 'string', {strcat(['Directory: ' directory_answer{1}]), ...
						strcat(['State: ' spectrum_list{spectrum_answer}]), ...
						'', ...
						strcat(['Rescale mass: ' num2str(mass_answer)]), ...
						strcat(['Rescale vev: ' num2str(vev_answer)]), ...
						'', ...
						strcat(['Fit Form: ' form_list{fit_form}]), ...
						strcat(['Fit Range: [' num2str(fit_minimum) ':' num2str(fit_maximum) ']']), ...
						strcat(['Fit X Precision: ' num2str(fit_x_prec)]), ...
						strcat(['Fold: ' fold_string]), ...
						strcat(['Parity Project: ' ppp_string]), ...
						strcat(['Zero Center: ' zero_string]), ...
						strcat(['Fit Points: ' even_string]), ...
						strcat(['Variance-Covariance: ' diag_string]), ...
						'', ...
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
                        '', ...
						strcat(['Chisq DoF: ' num2str(chisq_dof)]), ...
						strcat(['P-value: ' num2str(p_val)]), ...
						strcat(['Condition Num: ' num2str(cond_num)])});

	clear('fold_string'); clear('ppp_string'); clear('zero_string'); clear('even_string'); clear('diag_string'); clear('vstr'); clear('saving_string'); clear('blocking_string');
