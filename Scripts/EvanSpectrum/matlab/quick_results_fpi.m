% Get some scalar stuff!
function wall_results = quick_results_fpi(fname, state, parse_Ns, parse_Nt, blocksize, tmin, tmax, diag, fold, m_l, noerrors)

	% Check if we're ignoring error analysis.
	check_errs = 0;
	
	if (exist('noerrors', 'var'))
		check_errs = noerrors;
	end

	if (~strcmp(state, 'ps2'))
		error('The fpi analysis function should be looking at the ps2 state!');
	end
	
	% Now let's get some initial guesses. 
	[roots, masses, amps] = effective_mass(fname, state, parse_Nt, 1, 3, 0);
    
    tguess = max([tmin-2, 1]);
    
    cosh_mass_guess = abs(real(masses(tguess, 1)));
	cosh_amp_guess = real(amps(tguess, 1));
    
    % Load things up!
	
	% Load the correlator.
	wall_corr = load_correlator(fname, state, parse_Nt);
	
	% Fold it, if we want.
    if fold == 1
		% We fold it.
		wall_corr = fold_data(wall_corr, 0); % not a baryon!
	
    end
    
    %[mean_v(1) mean_v(2) mean_v(3) mean_v(4) mean_v(5)] = UWerr(transpose(wall_corr), 1.5, size(wall_corr,2), 'Name', tmin+1);
    
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
        
        % 0 means fit all, 1 means cosh or 3 means cosh+oscil.
        if fold == 0
            wall_output = get_all_nlfit(wall_sum, wall_cov_mat, tmin, parse_Nt-tmin, parse_Nt, 0, 1, diag, coefficients);
        else
            wall_output = get_all_nlfit(wall_sum, wall_cov_mat, tmin, tmax, parse_Nt, 0, 1, diag, coefficients);
        end

        wall_results = zeros(1, 23); % We have two extra values because of fpi.
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
		wall_results(1, 22) = 2*m_l*sqrt(wall_output(3)*cosh(parse_Nt*wall_output(4)/2)*((parse_Nt/2)^3)*3/(4*(wall_output(4)^3))); % f_pi

        % Jackknife, if needed.
		if check_errs == 0
		
			wall_coefficients_blocks = zeros(num_blocks, 9); % because 8 coefficients, in general, +1 for fpi.
			flag = 1;
			tic
			for b=1:num_blocks
				if (flag == 0)
					continue;
				end
				if fold == 0
					fit_output = get_all_nlfit(wall_jack(:,b), wall_cov_mat, tmin, parse_Nt-tmin, parse_Nt, 0, 1, diag, wall_coefficients);
				else
					fit_output = get_all_nlfit(wall_jack(:,b), wall_cov_mat, tmin, tmax, parse_Nt, 0, 1, diag, wall_coefficients);
				end

				if numel(fit_output) == 0
					flag = 0;
					continue;
				end

				wall_coefficients_blocks(b, 1:8) = fit_output(3:10);
				wall_coefficients_blocks(b, 9) = 2*m_l*sqrt(fit_output(3)*cosh(parse_Nt*fit_output(4)/2)*((parse_Nt/2)^3)*3/(4*(fit_output(4)^3))); % f_pi

			end
			toc
			if flag == 0
				disp(strcat('Failed on ', num2str(tmin)));
				return;
			end
			
			wall_coefficients_extend = zeros(1, 9);
			wall_coefficients_extend(1:8) = wall_coefficients;
			wall_coefficients_extend(9) = wall_results(1, 22);
			
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
			
		end
	

end