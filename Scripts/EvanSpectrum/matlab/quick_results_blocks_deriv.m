% Get results from jackknife blocks being passed in, instead
% of fnames and states.
% Takes optional arguments of blocked values.
function [wall_results, varargout] = quick_results_blocks_deriv(wall_sum, wall_jack, wall_cov_mat, fit_type, parse_Nt, tmin, tmax, diag, fold, noerrors)

    % Find out if we're returning additional info.
	nout = max(nargout,1) - 1; % Number of additional returns.
    
	% Check if we're ignoring error analysis.
	check_errs = 0;
	
	if (exist('noerrors', 'var'))
		check_errs = noerrors;
	end
	
    % Prepare basic information.
    num_blocks = size(wall_jack, 2);
    
    if (fit_type == 3 || fit_type == 13)
        oscil = 1;
    else
        oscil = 0;
    end

	% Now let's get some initial guesses. This is all heuristic, but works well enough!
	[masses, roots, amps] = effective_mass_utility(wall_sum, parse_Nt, 1+oscil, 3+2*oscil, 0);
    
    tguess = max([tmin-2, 1]);
    
    cosh_mass_guess = abs(real(masses(tguess, 1+oscil)));
	cosh_amp_guess = real(amps(tguess, 1+oscil));
    oscil_mass_guess = 0;
    oscil_amp_guess = 0;
    if (oscil == 1)
        oscil_mass_guess = abs(real(masses(tguess, 1)));
        oscil_amp_guess = -sign(cosh_amp_guess)*abs(real(amps(tguess, 1)));
    end
    
	% Let's try getting a central value for the wall.
	coefficients = zeros(1, 8);
	coefficients(1) = -cosh_amp_guess;
	coefficients(2) = cosh_mass_guess;
	coefficients(5) = oscil_amp_guess;
	coefficients(6) = oscil_mass_guess;

	% 0 means fit all, 1 means cosh or 3 means cosh+oscil.
    if (fold == 1)
    	wall_output = get_all_nlfit(wall_sum, wall_cov_mat, tmin, tmax, parse_Nt, 0, fit_type, diag, coefficients);
    else
        wall_output = get_all_nlfit(wall_sum, wall_cov_mat, tmin, parse_Nt-tmin, parse_Nt, 0, fit_type, diag, coefficients);
    end
	
	wall_results = zeros(1, 25); % We have four extra values because of the original coefficients.
	wall_coefficients = zeros(1, 8);

	if (numel(wall_output) == 0) % central fit failed.
		disp(strcat('Failed on ', num2str(tmin)));
		return;
	end
	
	wall_results(1, 1) = tmin; %tmin
    if fold == 1
        wall_results(1, 2) = tmax;
    else
        wall_results(1,2) = parse_Nt-tmin;
    end
	wall_results(1, 3) = wall_output(3);
	wall_results(1, 5) = wall_output(4); 
	wall_results(1, 7) = wall_output(5);
	wall_results(1, 9) = wall_output(6);
	wall_results(1, 11) = wall_output(7);
	wall_results(1, 13) = wall_output(8); 
	wall_results(1, 15) = wall_output(9);
	wall_results(1, 17) = wall_output(10);
	wall_coefficients(:) = wall_output(3:10);
	wall_results(1, 19) = wall_output(11); %chisqdof
	wall_results(1, 20) = wall_output(12); %p-value
	wall_results(1, 21) = wall_output(13); %condition #
	
	wall_results(1, 22) = wall_output(3)/(2*sinh(wall_output(4)/2)); % cosh.
	wall_results(1, 24) = wall_output(7)/(2*cosh(wall_output(8)/2)); % oscil.
	

	% Jackknife, if needed.
	if check_errs == 0
		wall_coefficients_blocks = zeros(num_blocks, 10); % because 8 coefficients, in general, +2 for new coeffs.
		flag = 1;
		tic
		for b=1:num_blocks
			if (flag == 0)
				continue;
			end
			
			if (fold == 1)
				fit_output = get_all_nlfit(wall_jack(:,b), wall_cov_mat, tmin, tmax, parse_Nt, 0, fit_type, diag, coefficients);
			else
				fit_output = get_all_nlfit(wall_jack(:,b), wall_cov_mat, tmin, parse_Nt-tmin+1, parse_Nt, 0, fit_type, diag, coefficients);
			end
		
		

			if numel(fit_output) == 0
				flag = 0;
				continue;
			end

			wall_coefficients_blocks(b, 1:8) = fit_output(3:10);
			wall_coefficients_blocks(b, 9) = fit_output(3)/(2*sinh(fit_output(4)/2)); % cosh.
			wall_coefficients_blocks(b, 10) = fit_output(7)/(2*cosh(fit_output(8)/2)); % oscil.
		end
		toc
		if flag == 0
			disp(strcat('Failed on ', num2str(tmin)));
			return;
		end
		
		wall_coefficients_extend = zeros(1,10);
		wall_coefficients_extend(1:8) = wall_coefficients;
		wall_coefficients_extend(9) = wall_results(1, 22);
		wall_coefficients_extend(10) = wall_results(1, 24);
		
		wall_coefficients_rep = repmat(wall_coefficients_extend, [num_blocks, 1]);
		wall_coefficients_err = sqrt(sum((wall_coefficients_blocks-wall_coefficients_rep).^2, 1).*(num_blocks-1)./num_blocks);
		wall_results(1,4) = wall_coefficients_err(1, 1);
		wall_results(1,6) = wall_coefficients_err(1, 2);
		wall_results(1,8) = wall_coefficients_err(1, 3);
		wall_results(1,10) = wall_coefficients_err(1, 4);
		wall_results(1,12) = wall_coefficients_err(1, 5);
		wall_results(1,14) = wall_coefficients_err(1, 6);
		wall_results(1,16) = wall_coefficients_err(1, 7);
		wall_results(1,18) = wall_coefficients_err(1, 8);
		wall_results(1,23) = wall_coefficients_err(1, 9);
		wall_results(1,25) = wall_coefficients_err(1, 10);

		% Check if we also send back more info.
		if (nout == 1)
			% Just send back the blocked values.
			varargout{1} = wall_coefficients_blocks;
		elseif (nout == 2)
			% Send back center and blocks.
			varargout{1} = wall_coefficients_extend;
			varargout{2} = wall_coefficients_blocks;
		end
		
	end

end