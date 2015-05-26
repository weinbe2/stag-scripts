% Get some scalar stuff!
function wall_results = quick_results_baryon(fname, state, parse_Ns, parse_Nt, blocksize, tmin, tmax, diag, fold, noerrors)

	% Check if we're ignoring error analysis.
	check_errs = 0;
	
	if (exist('noerrors', 'var'))
		check_errs = noerrors;
	end

	% Check if we're actually looking at a baryon.
	if (~(strcmp(state, 'nu') || strcmp(state, 'de')))
		disp('Not actually a baryon!')
		return;
    end
	
	% Now let's get some initial guesses. 
	[roots, roots_err, masses, masses_err, amps, amps_err] = effective_mass_err(fname, state, parse_Nt, 2, 5, 0, 0, 0, 2, 0);
    
	% Alright... this is painful.
	% Baryons are annoying. The cosh correction is because the fit assumes
	% an exponential, but the effective mass code assumes a cosh. The lack
	% of a factor of two is because at the center of the lattice, the
	% decaying and oscillating mode of a state contribute equally.
	
    tguess = max([tmin-2, 1]);
    
	pos_mass_guess = abs(real(masses(tguess, 2)));
	pos_amp_guess = real(amps(tguess, 2))*cosh(pos_mass_guess*parse_Nt/2);
	neg_mass_guess = abs(real(masses(tguess, 1)));
	neg_amp_guess = abs(real(amps(tguess, 1)))*cosh(neg_mass_guess*parse_Nt/2);
    
    % Load things up!
	
	% Load the correlator.
	wall_corr = load_correlator(fname, state, parse_Nt);
	
	% Fold it, if we want.
    if fold == 1
		% Default: We fold it.
		wall_corr = fold_data(wall_corr, 1); % It's a baryon!
	
    end
    
	% Block data!
	[wall_blocks num_blocks] = block_data(wall_corr, 2, blocksize);

	
    % Now that we have all of this in line, form our jackknife blocks and cov matrix.
    wall_sum = mean(wall_blocks, 2); % Use blocks!
	[wall_jack wall_cov_mat wall_err] = jackknife_from_blocks(wall_blocks);
    
    %disp([num2str(wall_sum(tmin+1)) '+/-' num2str(wall_err(tmin+1))])

        % Let's try getting a central value for the wall.
        coefficients = zeros(1, 8);
        coefficients(1) = pos_amp_guess;
        coefficients(2) = pos_mass_guess;
        coefficients(5) = neg_amp_guess;
        coefficients(6) = neg_mass_guess;

        % 0 means fit all, 12 means baryon!
        if fold == 1
            wall_output = get_all_nlfit(wall_sum, wall_cov_mat, tmin, tmax, parse_Nt, 0, 12, diag, coefficients);
        else
            wall_output = get_all_nlfit(wall_sum, wall_cov_mat, tmin, parse_Nt-tmin, parse_Nt, 0, 12, diag, coefficients);
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
				   
				% 0 means fit all, 12 means baryon!
				if fold == 1
					fit_output = get_all_nlfit(wall_jack(:,b), wall_cov_mat, tmin, tmax, parse_Nt, 0, 12, diag, coefficients);
				else
					fit_output = get_all_nlfit(wall_jack(:,b), wall_cov_mat, tmin, parse_Nt-tmin, parse_Nt, 0, 12, diag, coefficients);
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