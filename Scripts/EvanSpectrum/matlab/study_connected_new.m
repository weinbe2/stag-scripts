% Run analysis on a state. blocksize, tmin_min and tmin_max are optional arguments!
% If the optional arguments aren't supplied, it looks for a fitparams file.
% 08-20-2014: Load information from info.dat.
% 01-12-2015: Added 'noerrors' flag, which much be set to 0 or 1 if one wants to specify maxes and such.
%             Default behavior of noerrors is 0.
% 03-18-2015: New version. Runs faster, less io, everyone's happy.
% 04-17-2015: noerrors now is a bitmap. (diag only)(no errors) as above.
function outform = study_connected_new(fname, state, noerrors, blockval, tmin_min, tmin_max, tmax)
	% Load ensemble info.
	[fl_l, fl_h, parse_Ns, parse_Nt, beta, m_l, m_h] = load_info_file(fname);
    %[~, ~, parse_Ns, parse_Nt, ~, m_l, ~] = load_info_file(fname);
	
	min_search = 0; max_search = 0; tmax_search = 0; blocksize = 0;
	check_errs = 0; diag_only = 0;
	
	if (exist('noerrors', 'var'))
		if (mod(noerrors, 2) == 1)
			check_errs = 1;
			noerrors = noerrors - 1;
		end
		noerrors = noerrors / 2;
		if (noerrors == 1)
			diag_only = 1;
		end
	end
	
	% If we didn't get range arguments...
	if (~exist('blockval', 'var') && ~exist('tmin_min', 'var') && ~exist('tmin_max', 'var') && ~exist('tmax', 'var'))
		% See if a fitparams file exists.
		if (exist(fullfile(cd, strcat([fname '/spectrum2/fitparams/fitparam.' state])), 'file'))
			fitparam = importdata(strcat([fname '/spectrum2/fitparams/fitparam.' state]));
			min_search = fitparam(1);
			max_search = fitparam(3);
			tmax_search = fitparam(4);
			blocksize = fitparam(5);
			
		else
			error('A fitparams file is expected but does not exist!');
		end
	else
		min_search = tmin_min;
		max_search = tmin_max;
		tmax_search = tmax;
		blocksize = blockval;
	end
		
	mass_l = m_l;
  
	% Figure out what fit type we need.
	fit_type = 3; % default cosh+oscil
  
	% Learn about our measurement from the state name.
	is_baryon = 0;
	if (strcmp(state, 'nu') || strcmp(state, 'de'))
		is_baryon = 1;
		fit_type = 12; % baryon
	end
	
	is_fpi = 0;
	if (strcmp(state, 'ps2'))
		is_fpi = 1;
    end
	
	is_oscil = 1;
	if (strcmp(state, 'ps') || strcmp(state, 'ps2') || strcmp(state, 'i5'))
		is_oscil = 0;
		fit_type = 1; % single cosh
	end
	
	% Load the correlator, run autocorrelation, block it.
	connected = load_correlator(fname, state, parse_Nt);
	
	% Run autocorrelation and save it.
    
    % Standard things we output.
    wall_corr = connected;
    [the_mean the_error the_error_error the_tau the_tau_err] = UWerr(transpose(wall_corr), 1.5, size(wall_corr,2), 'Name', 1);
    disp(strcat(['1 ' num2str(the_mean) ' ' num2str(the_error) ' ' num2str(the_error_error) ' ' num2str(the_tau) ' ' num2str(the_tau_err)]))
    [the_mean the_error the_error_error the_tau the_tau_err] = UWerr(transpose(wall_corr), 1.5, size(wall_corr,2), 'Name', floor(parse_Nt/4)+1);
    disp(strcat([num2str(floor(parse_Nt/4)+1) ' ' num2str(the_mean) ' ' num2str(the_error) ' ' num2str(the_error_error) ' ' num2str(the_tau) ' ' num2str(the_tau_err)]))
    [the_mean the_error the_error_error the_tau the_tau_err] = UWerr(transpose(wall_corr), 1.5, size(wall_corr,2), 'Name', floor(parse_Nt/2)+1);
    disp(strcat([num2str(floor(parse_Nt/2)+1) ' ' num2str(the_mean) ' ' num2str(the_error) ' ' num2str(the_error_error) ' ' num2str(the_tau) ' ' num2str(the_tau_err)]))

    output_autocorr_info = zeros(parse_Nt, 6);
    for i=1:parse_Nt
        [the_mean the_error the_error_error the_tau the_tau_err] = UWerr(transpose(wall_corr), 1.5, size(wall_corr,2), 'Name', i);
        output_autocorr_info(i,1) = i-1; 
        output_autocorr_info(i,2) = the_mean;
        output_autocorr_info(i,3) = the_error;
        output_autocorr_info(i,4) = the_error_error;
        output_autocorr_info(i,5) = the_tau;
        output_autocorr_info(i,6) = the_tau_err;
    end
    full_fname = strcat(fname, '/spectrum2/uwerr/uwerr.', state);
    save(full_fname, 'output_autocorr_info', '-ascii');
		
	% End outputting autocorrelation info.
		
	% Fold the correlator, appropriately handling if it's a baryon.
    connected = fold_data(connected, is_baryon);

	% Block data!
	[connected_blocks, num_blocks] = block_data(connected, 2, blocksize);
		
	% Get jackknife blocks and covariance matrix.
	connected_sum = mean(connected_blocks, 2);
	[connected_jack, connected_cov_mat, connected_err] = jackknife_from_blocks(connected_blocks);
		
	% Save it!
    data = zeros((parse_Nt/2)+1,3);
	for i=1:((parse_Nt/2)+1)
		data(i,1) = i-1;
		data(i,2) = connected_sum(i);
		data(i,3) = connected_err(i);
	end
	full_fname = strcat(fname, '/spectrum2/sum/sum.', state);	
	save(full_fname, 'data', '-ascii');
	
	% Spit out some effective mass data.
	% This function computes and saves it.
	effective_mass_err(fname, state, parse_Nt, 1+is_oscil, 2+2*is_oscil, 0, 1, 1, blocksize, 1);
	% Let's prepare ourselves to get good initial guesses to fits, too.
	% The difference of 3 instead of 2 makes it an overcontrained effective mass, which tends to work better for guesses.
	[ roots, roots_err, masses, masses_err, amps, amps_err ] = effective_mass_err(fname, state, parse_Nt, 1+is_oscil, 3+2*is_oscil, 0, 1, 1, blocksize, 0);
	%[ roots, ~, masses, masses_err, amps, amps_err ] = effective_mass_err(fname, state, parse_Nt, 1+is_oscil, 3+2*is_oscil, 0, 1, 1, blocksize, 0);

    
    
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Get smarter initial guesses by looking between 5 and 10. %
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	cosh_mass_guess = 0;
	cosh_amp_guess = 0;
	cosh_amp_sign = 0; % save the sign for now.
	oscil_mass_guess = 0;
	oscil_amp_guess = 0;
	
	tguess_min = 5; tguess_max = 10;
	cm_mean = 0; cm_weight = 0;
	ca_mean = 0; ca_weight = 0; % get log amplitudes.
	om_mean = 0; om_weight = 0;
	oa_mean = 0; oa_weight = 0;
	num_mean = 0; num_weight = 0; % figure out where to start.
	for i=tguess_min:tguess_max
		if (is_oscil == 0) % Is it just a cosh?
			
			% Trust the earliest sign.
			if (cosh_amp_sign == 0)
				cosh_amp_sign = sign(amps(i, 1));
			end
			
			% We trust everything is sane.
			cm_mean = cm_mean + abs(real(masses(i,1)))/(masses_err(i,1).^2);
			cm_weight = cm_weight + 1/(masses_err(i,1).^2);
			ca_mean = ca_mean + log(abs(real(amps(i,1))))*abs(real(amps(i,1))).^2/(abs(real(amps_err(i,1))).^2);
			ca_weight = ca_weight + abs(real(amps(i,1))).^2/(abs(real(amps_err(i,1))).^2);
			
		else % The rest works the same, for now, whether
			 % or not it's a cosh+oscil or a baryon.
			
			% Trust the earliest sign.
			if (cosh_amp_sign == 0)
				cosh_amp_sign = sign(amps(i, 2));
			end
			
			% Make sure everything is sane, otherwise skip!
			if (real(roots(i,1)) >= -1 || abs(real(masses(i,1))) <= 1e-16 || abs(real(masses(i,2))) <= 1e-16)
				continue;
			end
			
			cm_mean = cm_mean + abs(real(masses(i,2)))/(masses_err(i,2).^2);
			cm_weight = cm_weight + 1/(masses_err(i,2).^2);
			ca_mean = ca_mean + log(abs(real(amps(i,2))))*abs(real(amps(i,2))).^2/(abs(real(amps_err(i,2))).^2);
			ca_weight = ca_weight + abs(real(amps(i,2))).^2/(abs(real(amps_err(i,2))).^2);
			
			om_mean = om_mean + abs(real(masses(i,1)))/(masses_err(i,1).^2);
			om_weight = om_weight + 1/(masses_err(i,1).^2);
			oa_mean = oa_mean + log(abs(real(amps(i,1))))*abs(real(amps(i,1))).^2/(abs(real(amps_err(i,1))).^2);
			oa_weight = oa_weight + abs(real(amps(i,1))).^2/(abs(real(amps_err(i,1))).^2);
			
		end
		
		num_mean = num_mean + i;
		num_weight = num_weight + 1;
	end
	
	% Make sure we got something!
	if (abs(cm_mean) < 1e-16) % Uh oh!
		disp('No valid initial guesses existed!');
		return;
	end
	
	% Get our initial guesses!
	cosh_mass_guess = cm_mean/cm_weight;
	cosh_amp_guess = cosh_amp_sign*exp(ca_mean/ca_weight);
	if (is_oscil == 1)
		oscil_mass_guess = om_mean/om_weight;
		oscil_amp_guess = -cosh_amp_sign*exp(oa_mean/oa_weight);
	end
	if (is_baryon == 1) % There's a normalization fix.
	% The cosh correction is because the fit assumes
	% an exponential, but the effective mass code assumes a cosh. The lack
	% of a factor of two is because at the center of the lattice, the
	% decaying and oscillating mode of a state contribute equally.
		cosh_amp_guess = cosh_amp_guess*cosh(cosh_mass_guess*parse_Nt/2);
		oscil_amp_guess = abs(oscil_amp_guess)*cosh(oscil_mass_guess*parse_Nt/2);
	end
	
	tguess = floor(num_mean/num_weight); % Start with tmin here!
	% There's never a reason this shouldn't be a real, but sometimes weird cases happen.
	guess_base = real([cosh_amp_guess, cosh_mass_guess, 0, 0, oscil_amp_guess, oscil_mass_guess, 0, 0]);
	guess = guess_base; 
	
	% debug
	guess_base 
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Start looking at tguess with the initial guesses, then %
	% work backwards first, then work forwards.              %
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	% Count gets computed a little bit differently now.
    % Start at tmin = 1, go until tmax - dof-1.
    %       tmax          [dof]              [tmin + 1]
    count = tmax_search - (2 + 2*is_oscil);
    
    output_ps = [];
	if (is_fpi == 1)
		output_ps = zeros(count, 23);
	else
		output_ps = zeros(count, 21);
	end
	
	for i=tguess:(-1):1 % Go all the way back to zero.
		
		% get_all_nlfit(corr_fcn, corr_mat, tmin, tmax, nt, ...
		% fiteven, fitfunc, fitdiag, coeff)
		fit_output = get_all_nlfit(connected_sum, connected_cov_mat, ...
						i, tmax_search, parse_Nt, 0, fit_type, diag_only, guess);
		
		if (numel(fit_output) == 0) % central fit failed.
			disp(strcat('Failed on ', num2str(i)));
			continue
		end
		
		% Save all the central value information.
		output_ps(i,1) = i;
		output_ps(i,2) = tmax_search;
		output_ps(i, 3) = fit_output(3);
		output_ps(i, 5) = fit_output(4); 
		output_ps(i, 7) = fit_output(5);
		output_ps(i, 9) = fit_output(6);
		output_ps(i, 11) = fit_output(7);
		output_ps(i, 13) = fit_output(8); 
		output_ps(i, 15) = fit_output(9);
		output_ps(i, 17) = fit_output(10);
		output_ps(i, 19) = fit_output(11); %chisqdof
		output_ps(i, 20) = fit_output(12); %p-value
		output_ps(i, 21) = fit_output(13); %condition #
		if (is_fpi == 1)
			output_ps(i, 22) = 2*mass_l*sqrt(fit_output(3)*cosh(parse_Nt*fit_output(4)/2)*((parse_Ns)^3)*3/(4*(fit_output(4)^3))); % f_pi
		end
		
		fit_coefficients = zeros(1, 9); % 9 in case of fpi.
		fit_coefficients(1:8) = fit_output(3:10);
		if (is_fpi == 1)
			fit_coefficients(9) = output_ps(i,22);
		end
		
		temp_guess = fit_coefficients(1:8);
		guess = temp_guess; 
		
		% Check and see if we're taking errors or not.
		if (check_errs == 0) % check them! 
			fit_coefficients_blocks = zeros(num_blocks, 9);
			flag = 1; % to check if a fit failed along the way.
			tic
			for b=1:num_blocks
				if (flag == 0)
					continue;
				end
				fit_output = get_all_nlfit(connected_jack(:,b), ...
				connected_cov_mat, i, tmax_search, parse_Nt, ...
				0, fit_type, diag_only, guess);
				
				if (numel(fit_output) == 0)
					flag = 0;
					continue;
				end
				
				fit_coefficients_blocks(b,1:8) = fit_output(3:10);
				if (is_fpi == 1)
					fit_coefficients_blocks(b,9) = 2*mass_l*sqrt(fit_output(3)*cosh(parse_Nt*fit_output(4)/2)*((parse_Ns)^3)*3/(4*(fit_output(4)^3))); % f_pi
				end
			end
			toc
			
			if (flag == 0)
				disp(strcat('Failed on ', num2str(i)));
				continue;
			end
			
			% Perform the jackknife.
			fit_coefficients_err = sqrt(sum((fit_coefficients_blocks-repmat(fit_coefficients, [num_blocks, 1])).^2, 1).*(num_blocks-1)./num_blocks);
			
			% And store the results in place.
			output_ps(i,4) = fit_coefficients_err(1, 1);
			output_ps(i,6) = fit_coefficients_err(1, 2);
			output_ps(i,8) = fit_coefficients_err(1, 3);
			output_ps(i,10) = fit_coefficients_err(1, 4);
			output_ps(i,12) = fit_coefficients_err(1, 5);
			output_ps(i,14) = fit_coefficients_err(1, 6);
			output_ps(i,16) = fit_coefficients_err(1, 7);
			output_ps(i,18) = fit_coefficients_err(1, 8);
			if (is_fpi == 1)
				output_ps(i,23) = fit_coefficients_err(1, 9);
			end
				
		end
	end
	
	% Okay, we got to the start, now let's run off to the end!
	% Get the right guess in place again.
	guess = [output_ps(tguess, 3), output_ps(tguess, 5), ...
		output_ps(tguess, 7), output_ps(tguess, 9), ...
		output_ps(tguess, 11), output_ps(tguess, 13), ...
		output_ps(tguess, 15), output_ps(tguess, 17)];
	
	for i=(tguess+1):count 
		
		% get_all_nlfit(corr_fcn, corr_mat, tmin, tmax, nt, ...
		% fiteven, fitfunc, fitdiag, coeff)
		fit_output = get_all_nlfit(connected_sum, connected_cov_mat, ...
						i, tmax_search, parse_Nt, 0, fit_type, diag_only, guess);
		
		if (numel(fit_output) == 0) % central fit failed.
			disp(strcat('Failed on ', num2str(i)));
			continue;
		end
		
		% Save all the central value information.
		output_ps(i,1) = i;
		output_ps(i,2) = tmax_search;
		output_ps(i, 3) = fit_output(3);
		output_ps(i, 5) = fit_output(4); 
		output_ps(i, 7) = fit_output(5);
		output_ps(i, 9) = fit_output(6);
		output_ps(i, 11) = fit_output(7);
		output_ps(i, 13) = fit_output(8); 
		output_ps(i, 15) = fit_output(9);
		output_ps(i, 17) = fit_output(10);
		output_ps(i, 19) = fit_output(11); %chisqdof
		output_ps(i, 20) = fit_output(12); %p-value
		output_ps(i, 21) = fit_output(13); %condition #
		if (is_fpi == 1)
			output_ps(i, 22) = 2*mass_l*sqrt(fit_output(3)*cosh(parse_Nt*fit_output(4)/2)*((parse_Ns)^3)*3/(4*(fit_output(4)^3))); % f_pi
		end
		
		fit_coefficients = zeros(1, 9); % 9 in case of fpi.
		fit_coefficients(1:8) = fit_output(3:10);
		if (is_fpi == 1)
			fit_coefficients(9) = output_ps(i,22);
		end
		
		temp_guess = fit_coefficients(1:8);
		guess = temp_guess; 
		
		% Check and see if we're taking errors or not.
		if (check_errs == 0) % check them! 
			fit_coefficients_blocks = zeros(num_blocks, 9);
			flag = 1; % to check if a fit failed along the way.
			tic
			for b=1:num_blocks
				if (flag == 0)
					continue;
				end
				fit_output = get_all_nlfit(connected_jack(:,b), ...
				connected_cov_mat, i, tmax_search, parse_Nt, ...
				0, fit_type, diag_only, guess);
				
				if (numel(fit_output) == 0)
					flag = 0;
					continue;
				end
				
				fit_coefficients_blocks(b,1:8) = fit_output(3:10);
				if (is_fpi == 1)
					fit_coefficients_blocks(b,9) = 2*mass_l*sqrt(fit_output(3)*cosh(parse_Nt*fit_output(4)/2)*((parse_Ns)^3)*3/(4*(fit_output(4)^3))); % f_pi
				end
			end
			toc
			
			if (flag == 0)
				disp(strcat('Failed on ', num2str(i)));
				continue;
			end
			
			% Perform the jackknife.
			fit_coefficients_err = sqrt(sum((fit_coefficients_blocks-repmat(fit_coefficients, [num_blocks, 1])).^2, 1).*(num_blocks-1)./num_blocks);
			
			% And store the results in place.
			output_ps(i,4) = fit_coefficients_err(1, 1);
			output_ps(i,6) = fit_coefficients_err(1, 2);
			output_ps(i,8) = fit_coefficients_err(1, 3);
			output_ps(i,10) = fit_coefficients_err(1, 4);
			output_ps(i,12) = fit_coefficients_err(1, 5);
			output_ps(i,14) = fit_coefficients_err(1, 6);
			output_ps(i,16) = fit_coefficients_err(1, 7);
			output_ps(i,18) = fit_coefficients_err(1, 8);
			if (is_fpi == 1)
				output_ps(i,23) = fit_coefficients_err(1, 9);
			end
				
		end
	end
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Save everything we got out to file! %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    full_fname = strcat(fname, '/spectrum2/fits/fit_new.', state);
    save(full_fname, 'output_ps', '-ascii');
    
    % Also just save fit masses (and fpi, in the case of fpi. Man fpi is
    % annoying.
    
    % Cosh state.
    raw_output = zeros(size(output_ps,1), 3);
    for j=1:(size(output_ps,1))
        raw_output(j,1) = real(output_ps(j,1));
        raw_output(j,2) = real(output_ps(j, 5));
        raw_output(j,3) = real(output_ps(j, 6));
    end
    full_fname = strcat(fname, '/spectrum2/fits/fit.', state, num2str(1));
    save(full_fname, 'raw_output', '-ascii');

    % Oscil state, if it exists.
    if (is_oscil == 1)
        raw_output = zeros(size(output_ps,1), 3);
        for j=1:(size(output_ps,1))
            raw_output(j,1) = real(output_ps(j,1));
            raw_output(j,2) = real(output_ps(j, 13));
            raw_output(j,3) = real(output_ps(j, 14));
        end
        full_fname = strcat(fname, '/spectrum2/fits/fit.', state, num2str(2));
        save(full_fname, 'raw_output', '-ascii');
    end
    
    % And don't forget fpi.
    if (is_fpi == 1)
        raw_output = zeros(size(output_ps,1), 3);
        for j=1:(size(output_ps,1))
            raw_output(j,1) = real(output_ps(j,1));
            raw_output(j,2) = real(output_ps(j, 22));
            raw_output(j,3) = real(output_ps(j, 23));
        end
        full_fname = strcat(fname, '/spectrum2/fits/fit.', state, num2str(2));
        save(full_fname, 'raw_output', '-ascii');
    end
    
	outform = output_ps;

end
