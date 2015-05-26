% Get results from jackknife blocks being passed in, instead
% of fnames and states.
function wall_results = quick_results_blocks(wall_sum, wall_jack, wall_cov_mat, oscil, parse_Nt, tmin, tmax, diag, fold, noerrors)

	% Check if we're ignoring error analysis.
	check_errs = 0;
	
	if (exist('noerrors', 'var'))
		check_errs = noerrors;
	end

    % Prepare basic information.
    num_blocks = size(wall_jack, 2);

	% Now let's get some initial guesses. This is all heuristic, but works well enough!
	if (oscil == 1 || oscil == 2)
        state_count = 2;
        guess_count = 5;
    else
        state_count = 1;
        guess_count = 3;
    end
	[masses, roots, amps] = effective_mass_utility(wall_sum, parse_Nt, state_count, guess_count, 0);
    
    tguess = max([tmin-2, 1]);
    
    cosh_mass_guess = abs(real(masses(tguess, state_count)));
	cosh_amp_guess = real(amps(tguess, state_count));
    oscil_mass_guess = 0;
    oscil_amp_guess = 0;
    if (oscil == 1 || oscil == 2)
        oscil_mass_guess = abs(real(masses(tguess, 1)));
        oscil_amp_guess = -sign(cosh_amp_guess)*abs(real(amps(tguess, 1)));
    end
    	
	% Fold it, if we want.
    if fold == 1
		% We fold it.
		wall_sum = fold_data(wall_sum, 0); % 0 b/c not a baryon.
		wall_jack = fold_data(wall_jack, 0); 
	
    end
    
	% Let's try getting a central value for the wall.
	coefficients = zeros(1, 8);
	coefficients(1) = cosh_amp_guess;
	coefficients(2) = cosh_mass_guess;
	coefficients(5) = oscil_amp_guess;
	coefficients(6) = oscil_mass_guess;

	% 0 means fit all, 1 means cosh or 3 means cosh+oscil or 7 means cosh+oscil+const
	if fold == 0
		wall_output = get_all_nlfit(wall_sum, wall_cov_mat, tmin, parse_Nt-tmin, parse_Nt, 0, 2^(oscil+1)-1, diag, coefficients);
	else
		wall_output = get_all_nlfit(wall_sum, wall_cov_mat, tmin, tmax, parse_Nt, 0, 2^(oscil+1)-1, diag, coefficients);
	end

	wall_results = zeros(1, 21);
	wall_coefficients = zeros(1, 8);

	if (numel(wall_output) == 0) % central fit failed.
		disp(strcat('Failed on ', num2str(tmin)));
		return;
	end
	
	wall_results(1, 1) = tmin; %tmin
	if fold == 0
		wall_results(1, 2) = parse_Nt-tmin;
	else
		wall_results(1, 2) = tmax;
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

	% Jackknife, if needed.
	if check_errs == 0
			
		wall_coefficients_blocks = zeros(num_blocks, 8); % because 8 coefficients, in general.
		flag = 1;
		tic
		for b=1:num_blocks
			if (flag == 0)
				continue;
			end
			if fold == 0
				fit_output = get_all_nlfit(wall_jack(:,b), wall_cov_mat, tmin, parse_Nt-tmin, parse_Nt, 0, 1+2*oscil, diag, wall_coefficients);
			else
				fit_output = get_all_nlfit(wall_jack(:,b), wall_cov_mat, tmin, tmax, parse_Nt, 0, 1+2*oscil, diag, wall_coefficients);
			end

			if numel(fit_output) == 0
				flag = 0;
				continue;
			end

			wall_coefficients_blocks(b, :) = fit_output(3:10);

		end
		toc
		if flag == 0
			disp(strcat('Failed on ', num2str(tmin)));
			return;
		end
		wall_coefficients_rep = repmat(wall_coefficients, [num_blocks, 1]);
		wall_coefficients_err = sqrt(sum((wall_coefficients_blocks-wall_coefficients_rep).^2, 1).*(num_blocks-1)./num_blocks);
		wall_results(1,4) = wall_coefficients_err(1, 1);
		wall_results(1,6) = wall_coefficients_err(1, 2);
		wall_results(1,8) = wall_coefficients_err(1, 3);
		wall_results(1,10) = wall_coefficients_err(1, 4);
		wall_results(1,12) = wall_coefficients_err(1, 5);
		wall_results(1,14) = wall_coefficients_err(1, 6);
		wall_results(1,16) = wall_coefficients_err(1, 7);
		wall_results(1,18) = wall_coefficients_err(1, 8);

	end
	

end