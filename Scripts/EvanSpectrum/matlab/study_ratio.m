% Study the ratio of two mass values.
% Created 10/28/2014
function outform = study_ratio(fname, state1, statenum1, state2, statenum2, blockval1, tmin1, tmax1, blockval2, tmin2, tmax2)
	% Load ensemble info.
	%[fl_l, fl_h, parse_Ns, parse_Nt, beta, m_l, m_h] = load_info_file(fname);
    [~, ~, parse_Ns, parse_Nt, ~, m_l, ~] = load_info_file(fname);
	
	blocksize1 = 0; tmin_val1 = 0; tmax_val1 = 0;
	blocksize2 = 0; tmin_val2 = 0; tmax_val2 = 0;
	
	% If we didn't get range arguments...
	if (~exist('blockval1', 'var') && ~exist('tmin1', 'var') && ~exist('tmax1', 'var') && ~exist('blockval2', 'var') && ~exist('tmin2', 'var') && ~exist('tmax2', 'var'))
		% See if a fitparams file exists.
		if (exist(fullfile(cd, strcat([fname '/spectrum2/fitparams/fitparam.' state1])), 'file'))
			fitparam = importdata(strcat([fname '/spectrum2/fitparams/fitparam.' state1]));
			tmin_val1 = fitparam(2);
			tmax_val1 = fitparam(4);
			blocksize1 = fitparam(5);
			
		else
			error('A fitparams file is expected but does not exist!');
		end
		
		if (exist(fullfile(cd, strcat([fname '/spectrum2/fitparams/fitparam.' state2])), 'file'))
			fitparam = importdata(strcat([fname '/spectrum2/fitparams/fitparam.' state2]));
			tmin_val2 = fitparam(2);
			tmax_val2 = fitparam(4);
			blocksize2 = fitparam(5);
			
		else
			error('A fitparams file is expected but does not exist!');
		end
	else
		blocksize1 = blockval1;
		tmin_val1 = tmin1;
		tmax_val1 = tmax1;
		
		blocksize2 = blockval2;
		tmin_val2 = tmin2;
		tmax_val2 = tmax2;
	end
		
	tmin1 = tmin_val1;
	tmin2 = tmin_val2;
	tmax1 = tmax_val1;
	tmax2 = tmax_val2;
	
	mass_l = m_l;
  
	% Learn about our measurement from the state name.
	% We're not dealing with baryons here.
	if (strcmp(state1, 'nu') || strcmp(state1, 'de') || strcmp(state2, 'nu') || strcmp(state2, 'de'))
		error('This function does not support baryons.');
	end
	
	% Next, for the tmin and tmax, we get fit values.
	% This part is borrowed from quick_results.
	
	% We need to get some guesses. This depends on if we have an oscillating mode or not.
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% State 1
	oscil = 0;
	if (strcmp(state1, 'ps') || strcmp(state1, 'ps2') || strcmp(state1, 'i5'))
		oscil = 1;
	end
	
	% Now let's get some initial guesses. 
	[roots, masses, amps] = effective_mass(fname, state1, parse_Nt, 2-oscil, 5-2*oscil, 0);
    
    tguess = max([tmin1-2, 1]);
    
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
	wall_corr = load_correlator(fname, state1, parse_Nt);
	
	% Fold it
	wall_corr = fold_data(wall_corr, 0); % Not a baryon
    
	% Block data!
	[wall_blocks num_blocks] = block_data(wall_corr, 2, blocksize1);
	
    % Now that we have all of this in line, form our jackknife blocks and cov matrix.
    wall_sum = mean(wall_blocks, 2); % Use blocks!
	[wall_jack wall_cov_mat wall_err] = jackknife_from_blocks(wall_blocks);

	% Let's try getting a central value for the wall.
	coefficients = zeros(1, 8);
	coefficients(1) = cosh_amp_guess;
	coefficients(2) = cosh_mass_guess;
	coefficients(5) = oscil_amp_guess;
	coefficients(6) = oscil_mass_guess;

	% 0 means fit all, 1 means cosh or 3 means cosh+oscil.
	wall_output = get_all_nlfit(wall_sum, wall_cov_mat, tmin1, tmax1, parse_Nt, 0, 3-2*oscil, 0, coefficients);

	wall_results = zeros(1, 21);
	wall_coefficients = zeros(1, 8);

	if (numel(wall_output) == 0) % central fit failed.
		disp(strcat('Failed on ', num2str(tmin1)));
		return;
	end
	
	wall_results(1, 1) = tmin1; %tmin
	wall_results(1, 2) = tmax1;
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

	% Jackknife. 
	wall_coefficients_blocks = zeros(num_blocks, 8); % because 8 coefficients, in general.
	flag = 1;
	tic
	for b=1:num_blocks
		if (flag == 0)
			continue;
		end
		fit_output = get_all_nlfit(wall_jack(:,b), wall_cov_mat, tmin1, tmax1, parse_Nt, 0, 3-2*oscil, 0, wall_coefficients);

		if numel(fit_output) == 0
			flag = 0;
			continue;
		end

		wall_coefficients_blocks(b, :) = fit_output(3:10);

	end
	toc
	if flag == 0
		disp(strcat('Failed on ', num2str(tmin1)));
		return;
	end
	
	state1_info = wall_coefficients;
	state1_blocks = wall_coefficients_blocks;
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% State 2
	
		% We need to get some guesses. This depends on if we have an oscillating mode or not.
	oscil = 0;
	if (strcmp(state2, 'ps') || strcmp(state2, 'ps2') || strcmp(state2, 'i5'))
		oscil = 1;
	end
	
	% Now let's get some initial guesses. 
	[roots, masses, amps] = effective_mass(fname, state2, parse_Nt, 2-oscil, 5-2*oscil, 0);
    
    tguess = max([tmin2-2, 1]);
    
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
	wall_corr = load_correlator(fname, state2, parse_Nt);
	
	% Fold it
	wall_corr = fold_data(wall_corr, 0); % Not a baryon
    
	% Block data!
	[wall_blocks num_blocks] = block_data(wall_corr, 2, blocksize2);
	
    % Now that we have all of this in line, form our jackknife blocks and cov matrix.
    wall_sum = mean(wall_blocks, 2); % Use blocks!
	[wall_jack wall_cov_mat wall_err] = jackknife_from_blocks(wall_blocks);

	% Let's try getting a central value for the wall.
	coefficients = zeros(1, 8);
	coefficients(1) = cosh_amp_guess;
	coefficients(2) = cosh_mass_guess;
	coefficients(5) = oscil_amp_guess;
	coefficients(6) = oscil_mass_guess;

	% 0 means fit all, 1 means cosh or 3 means cosh+oscil.
	wall_output = get_all_nlfit(wall_sum, wall_cov_mat, tmin2, tmax2, parse_Nt, 0, 3-2*oscil, 0, coefficients);

	wall_results = zeros(1, 21);
	wall_coefficients = zeros(1, 8);

	if (numel(wall_output) == 0) % central fit failed.
		disp(strcat('Failed on ', num2str(tmin2)));
		return;
	end
	
	wall_results(1, 1) = tmin2; %tmin
	wall_results(1, 2) = tmax2;
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

	% Jackknife. 
	wall_coefficients_blocks = zeros(num_blocks, 8); % because 8 coefficients, in general.
	flag = 1;
	tic
	for b=1:num_blocks
		if (flag == 0)
			continue;
		end
		fit_output = get_all_nlfit(wall_jack(:,b), wall_cov_mat, tmin2, tmax2, parse_Nt, 0, 3-2*oscil, 0, wall_coefficients);

		if numel(fit_output) == 0
			flag = 0;
			continue;
		end

		wall_coefficients_blocks(b, :) = fit_output(3:10);

	end
	toc
	if flag == 0
		disp(strcat('Failed on ', num2str(tmin2)));
		return;
	end
	
	state2_info = wall_coefficients;
	state2_blocks = wall_coefficients_blocks;
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Now, build up the ratio.
	
	mass1 = state1_info(1,-2+statenum1*4);
	mass2 = state2_info(1,-2+statenum2*4);
	
	mass1_blocks = state1_blocks(:,-2+statenum1*4);
	mass2_blocks = state2_blocks(:,-2+statenum2*4);
	
	ratio_central = mass1./mass2;
	ratio_blocks = mass1_blocks./mass2_blocks;
	
	ratio_rep = repmat(ratio_central, [num_blocks 1]);
	ratio_err = sqrt(sum((ratio_rep-ratio_blocks).^2,1).*(num_blocks-1)./num_blocks);
    
    full_fname = strcat(fname, '/spectrum2/central/central_ratio.', state1, num2str(statenum1), '_', state2, num2str(statenum2));
	
	ratio_vals = zeros(1, 2);
	ratio_vals(1) = ratio_central;
	ratio_vals(2) = ratio_err;
    save(full_fname, 'ratio_vals', '-ascii');
    
	outform = ratio_vals;
end
