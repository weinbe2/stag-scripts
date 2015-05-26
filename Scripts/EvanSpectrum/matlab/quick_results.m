% Get some scalar stuff!
function wall_results = quick_results(fname, state, parse_Ns, parse_Nt, blocksize, tmin, tmax, diag, fold, noerrors)

	% Check if we're ignoring error analysis.
	check_errs = 0;
	
	if (exist('noerrors', 'var'))
		check_errs = noerrors;
	end
	

	% We need to get some guesses. This depends on if we have an oscillating mode or not.
	oscil = 0;
	if (strcmp(state, 'ps') || strcmp(state, 'ps2') || strcmp(state, 'i5') || strcmp(state, 'ps_wt0_10') || strcmp(state, 'ps2_wt0_10') || strcmp(state, 'i5_wt0_10') || strcmp(state, 'ps_wt0_30') || strcmp(state, 'ps2_wt0_30') || strcmp(state, 'i5_wt0_30'))
		oscil = 1;
	end
	
	% Now let's get some initial guesses. 
	[roots, masses, amps] = effective_mass(fname, state, parse_Nt, 2-oscil, 5-2*oscil, 0);
    
    tguess = max([tmin-2, 1]);
    
    cosh_mass_guess = abs(real(masses(tguess, 2-oscil)));
	cosh_amp_guess = real(amps(tguess, 2-oscil));
    oscil_mass_guess = 0;
    oscil_amp_guess = 0;
    if (oscil == 0)
        oscil_mass_guess = abs(real(masses(tguess, 1)));
        oscil_amp_guess = -sign(cosh_amp_guess)*abs(real(amps(tguess, 1)));
    end
    
    % Load things up!
	
	% Load the correlator.
	wall_corr = load_correlator(fname, state, parse_Nt);
	
	% Fold it, if we want.
    if fold == 1
		% Check if it's a Baryon. This only effects folding.
		is_baryon = 0;
		if (strcmp(state,'nu') || strcmp(state,'de'))
			is_baryon = 1;
		end

		% We fold it.
		wall_corr = fold_data(wall_corr, is_baryon);
	
    end

	% Block data!
	[wall_blocks num_blocks] = block_data(wall_corr, 2, blocksize);
    
	
	
    % Now that we have all of this in line, form our jackknife blocks and cov matrix.
    wall_sum = mean(wall_blocks, 2); % Use blocks!
	[wall_jack wall_cov_mat wall_err] = jackknife_from_blocks(wall_blocks);
    
    %disp([num2str(wall_sum(tmin+1)) '+/-' num2str(wall_err(tmin+1))])

        % Let's try getting a central value for the wall.
        coefficients = zeros(1, 8);
        coefficients(1) = cosh_amp_guess;
        coefficients(2) = cosh_mass_guess;
        coefficients(5) = oscil_amp_guess;
        coefficients(6) = oscil_mass_guess;

        % 0 means fit all, 1 means cosh or 3 means cosh+oscil.
        if fold == 0
            wall_output = get_all_nlfit(wall_sum, wall_cov_mat, tmin, parse_Nt-tmin, parse_Nt, 0, 3-2*oscil, diag, coefficients);
        else
            wall_output = get_all_nlfit(wall_sum, wall_cov_mat, tmin, tmax, parse_Nt, 0, 3-2*oscil, diag, coefficients);
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
					fit_output = get_all_nlfit(wall_jack(:,b), wall_cov_mat, tmin, parse_Nt-tmin, parse_Nt, 0, 3-2*oscil, diag, wall_coefficients);
				else
					fit_output = get_all_nlfit(wall_jack(:,b), wall_cov_mat, tmin, tmax, parse_Nt, 0, 3-2*oscil, diag, wall_coefficients);
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
