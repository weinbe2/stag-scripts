% 2015-01-05 Test cosh+oscil, cosh+oscil+const, cosh+oscil diff fits. 
function test_zero_fit(fname, blocksize, diag, fold, just_d, tmax_in)
	% Load ensemble info.
	%[fl_l, fl_h, parse_Ns, parse_Nt, beta, m_l, m_h] = load_info_file(fname);
    [fl_l, fl_h, parse_Ns, parse_Nt, beta, m_l, m_h] = load_info_file(fname);
	fl_flavor = fl_l/4; % Staggered counting.
    
	if (just_d == 0)
		state = 'sg_stoch';
	else
		state = 'dc_stoch';
	end
	
	mass_l = m_l;
	number_bl = 6;
	
	% Well, look, this covariance matrix is wrong.
	%{
	% Load
	correlator = load_correlator(fname, state, parse_Nt);
	
	% Block, etc.
	[correlator_blocks, num_blocks] = block_data(correlator, 2, blocksize);
	
	% Get jackknife blocks and covariance matrix.
	correlator_sum = mean(correlator_blocks, 2);
	[correlator_jack, correlator_cov_mat, correlator_err] = jackknife_from_blocks(correlator_blocks);
	
	correlator_cov_mat
	%}
	
	% Steal the cov matrix from 'study_disconnected'
	[disc_sum, disc_jack, disc_cov_mat, disc_err] = load_correlator_scalar(fname, fl_flavor, parse_Nt, parse_Ns, number_bl, blocksize);
	num_blocks = size(disc_jack, 2);
	
	% Pull out sigma if we need to.
	if (strcmp(state, 'sg_stoch'))
		connected = load_correlator(fname, 'sc_stoch', parse_Nt);

		connected = fold_data(connected, 0); % 0 = not a baryon.

		% Block data!
		[connected_blocks, num_blocks] = block_data(connected, 2, blocksize);

		% Get jackknife blocks and covariance matrix.
		connected_sum = mean(connected_blocks, 2);
		[connected_jack, connected_cov_mat, connected_err] = jackknife_from_blocks(connected_blocks);
	
		disc_sum = disc_sum - connected_sum; % we put in the factor of N_f/4 earlier.
		disc_jack = disc_jack - connected_jack;
		disc_rep = repmat(disc_sum, [1 num_blocks]);
		
		disc_cov_mat = zeros(parse_Nt);
		disc_err = zeros(parse_Nt, 1);
		for t1 = 1:parse_Nt
			for t2 = 1:parse_Nt
				disc_cov_mat(t1,t2) = sum((disc_rep(t1,:)-disc_jack(t1,:)).*(disc_rep(t2,:)-disc_jack(t2,:)),2).*(num_blocks-1)./num_blocks;
				if t1 == t2
					disc_err(t1) = sqrt(disc_cov_mat(t1,t1));
				end
			end
		end
	
	end
	
	% And rename because I'm lazy.
	correlator_sum = disc_sum;
	correlator_jack = disc_jack;
	correlator_cov_mat = disc_cov_mat;
	correlator_err = disc_err;
	
	% Diagonalize the cov matrix if we don't care.
	if (diag == 1)
		correlator_cov_mat = diag(diag(correlator_cov_mat));
	end
	
    % Let's get some good guesses to work with. 
	[masses, root, amps] = effective_mass_utility(correlator_sum, parse_Nt, 2, 5, 0);
	
	% We assume the guess at 10 is good, and if it's not
	% we work back from there.
	tguess = 10;
	flag = 0;
	a_shift = 0;
	
	while (flag == 0 && tguess > 5) % don't go too far back.
		% Make sure '1' corresponds to an oscillating, and that
		% both masses are positive.
		if (real(root(tguess-a_shift, 1)) < -1 && abs(real(masses(tguess-a_shift, 1))) > 1e-16 && abs(real(masses(tguess-a_shift, 2))) > 1e-16)
			flag = 1;
		else
			tguess = tguess - 1;
		end
	end
	
	if (flag == 0) % if we never found a mass, well, eh.
		[masses, root, amps] = effective_mass_utility(correlator_sum, parse_Nt, 1, 3, 0); % Assume a cosh?
		tguess = 10;
		cosh_mass_guess = abs(real(masses(tguess-a_shift, 1)));
		cosh_amp_guess = real(amps(tguess-a_shift, 1));
		% aaaand make something up.
		oscil_mass_guess = cosh_mass_guess;
		oscil_amp_guess = -0.05*cosh_amp_guess;
	else % we found something!
		cosh_mass_guess = abs(real(masses(tguess-a_shift, 2)));
		cosh_amp_guess = real(amps(tguess-a_shift, 2));
		oscil_mass_guess = abs(real(masses(tguess-a_shift, 1)));
        oscil_amp_guess = -sign(cosh_amp_guess)*abs(real(amps(tguess-a_shift, 1)));
	end
	
	% Now that we have the guess, we can subtract the center of the mean correlator. 
	center_to_sub = correlator_sum(parse_Nt/2+1);
	correlator_sum = correlator_sum - center_to_sub;
	correlator_jack = correlator_jack - center_to_sub;
		
    % Alright, pick some min t and max t to look at.
	tmin = 2;
        if (exist('tmax_in', 'var'))
		tmax = tmax_in;
	else
		tmax = 16;
	end

	% Set up what we need to do the fits.
	% This is the anonymous function for a cosh+const+oscil fit.
	guess_base = [cosh_amp_guess, cosh_mass_guess, oscil_amp_guess, oscil_mass_guess, -cosh_amp_guess];
	guess = guess_base; 
	
	% debug
	guess_base
	
	% Prepare the output.
	output_wconst = zeros(tmax-tmin+1, 23);
	nt = parse_Nt; 
	
	% Alrighty! First start at tguess and work back to tmin.
	for i=tguess:(-1):tmin
		
		disp(num2str(i))
		
		t1 = i;
		if (fold == 1)
			t2 = parse_Nt/2;
		else
			t2 = parse_Nt-i;
		end
		%t2 = 12;
		
		xval = t1:t2;
		yval = correlator_sum((t1+1):(t2+1))';
		ycorr = correlator_cov_mat((t1+1):(t2+1),(t1+1):(t2+1));
		output_wconst(i-tmin+1, 21) = cond(ycorr);
		%log(eig(ycorr))/log(10)
		
		chisqfunc = @(x)(1.0/(size(yval,2)-5)*(yval(:)-x(1)*cosh(x(2)*(nt/2-xval(:)))-x(3)*(1-2*mod(xval(:),2)).*cosh(x(4)*(nt/2-xval(:)))-x(5))'*(ycorr\(yval(:)-x(1)*cosh(x(2)*(nt/2-xval(:)))-x(3)*(1-2*mod(xval(:),2)).*cosh(x(4)*(nt/2-xval(:)))-x(5))));
		
		% Fit the center.
		% #19 is the chisqperdof.
		[temp_guess, output_wconst(i-tmin+1, 19), succ_code] = fminsearch(chisqfunc, guess, optimset('TolX', 1e-10, 'TolFun', 1e-7, 'MaxFunEvals', 100000, 'MaxIter', 100000)); 
		if (succ_code == 0) % central failed
			disp(strcat('Failed on ', num2str(t1)));
			continue;
		end
		
		temp_guess(2) = abs(temp_guess(2));
		temp_guess(4) = abs(temp_guess(4));
		
	   % debug
	   % temp_guess
		
		if (output_wconst(i-tmin+1,19) > 1e-20 && succ_code == 1)
			output_wconst(i-tmin+1, 1) = t1;
			output_wconst(i-tmin+1, 2) = t2;
			output_wconst(i-tmin+1, 3) = temp_guess(1);
			output_wconst(i-tmin+1, 5) = abs(temp_guess(2));
			output_wconst(i-tmin+1, 11) = temp_guess(3);
			output_wconst(i-tmin+1, 13) = abs(temp_guess(4));
			output_wconst(i-tmin+1, 7) = temp_guess(5);
			output_wconst(i-tmin+1, 22) = temp_guess(1)*cosh(temp_guess(2)*nt/2);
			
		end
		
		output_wconst(i-tmin+1, 20) = 1.0-chi2cdf(output_wconst(i-tmin+1,19)*(numel(yval)-numel(temp_guess)), numel(yval)-numel(temp_guess));
		
		% Save the coefficients for a jackknife.
		coefficient_center = temp_guess;
		% New 01-26-2015: Keep amp*cosh(mass(Nt/2))
		coefficient_center(6) = temp_guess(1)*cosh(temp_guess(2)*nt/2);
		coefficient_blocks = zeros(num_blocks, numel(temp_guess)+1);
		guess = temp_guess; % for next pass. 
		
		flag = 1;
		tic
		for b=1:num_blocks
			t1 = i;
			if (fold == 1)
				t2 = parse_Nt/2;
			else
				t2 = parse_Nt-i;
			end
			
			xval = t1:t2;
			yval = correlator_jack((t1+1):(t2+1), b)';
			ycorr = correlator_cov_mat((t1+1):(t2+1),(t1+1):(t2+1));
			
			chisqfunc = @(x)(1.0/(size(yval,2)-5)*(yval(:)-x(1)*cosh(x(2)*(nt/2-xval(:)))-x(3)*(1-2*mod(xval(:),2)).*cosh(x(4)*(nt/2-xval(:)))-x(5))'*(ycorr\(yval(:)-x(1)*cosh(x(2)*(nt/2-xval(:)))-x(3)*(1-2*mod(xval(:),2)).*cosh(x(4)*(nt/2-xval(:)))-x(5))));
			
			% Fit the center.
			% #19 is the chisqperdof.
			[temp_guess, temp, succ_code] = fminsearch(chisqfunc, coefficient_center, optimset('TolX', 1e-15, 'TolFun', 1e-7, 'MaxFunEvals', 100000, 'MaxIter', 100000)); 
			if (succ_code == 0)
				flag = 0;
				continue;
			end
			
			temp_guess(2) = abs(temp_guess(2));
			temp_guess(4) = abs(temp_guess(4));
			coefficient_blocks(b,1:5) = temp_guess(1:5); 
			coefficient_blocks(b,6) = temp_guess(1)*cosh(temp_guess(2)*nt/2);
			
		end % b
		toc
		
		if flag == 0
			disp(strcat('Failed on ', num2str(t1)));
			continue;
		end
		coefficient_rep = repmat(coefficient_center, [num_blocks, 1]);
		coefficient_err = sqrt(sum((coefficient_blocks-coefficient_rep).^2, 1).*(num_blocks-1)./num_blocks);
		
		output_wconst(i-tmin+1, 4) = coefficient_err(1);
		output_wconst(i-tmin+1, 6) = coefficient_err(2);
		output_wconst(i-tmin+1, 12) = coefficient_err(3);
		output_wconst(i-tmin+1, 14) = coefficient_err(4);
		output_wconst(i-tmin+1, 8) = coefficient_err(5);
		output_wconst(i-tmin+1, 23) = coefficient_err(6);
		
	end
	
	guess = [output_wconst(tguess-tmin+1, 3), output_wconst(tguess-tmin+1, 5), ...
		output_wconst(tguess-tmin+1, 11), output_wconst(tguess-tmin+1, 13), ...
		output_wconst(tguess-tmin+1, 7)];
	
	% Alrighty! Now start at tguess+1 and work to tmax.
	for i=(tguess+1):tmax;
		
		disp(num2str(i))
		
		t1 = i;
		if (fold == 1)
			t2 = parse_Nt/2;
		else
			t2 = parse_Nt-i;
		end
		
		xval = t1:t2;
		yval = correlator_sum((t1+1):(t2+1))';
		ycorr = correlator_cov_mat((t1+1):(t2+1),(t1+1):(t2+1));
		output_wconst(i-tmin+1, 21) = cond(ycorr);
		
		chisqfunc = @(x)(1.0/(size(yval,2)-5)*(yval(:)-x(1)*cosh(x(2)*(nt/2-xval(:)))-x(3)*(1-2*mod(xval(:),2)).*cosh(x(4)*(nt/2-xval(:)))-x(5))'*(ycorr\(yval(:)-x(1)*cosh(x(2)*(nt/2-xval(:)))-x(3)*(1-2*mod(xval(:),2)).*cosh(x(4)*(nt/2-xval(:)))-x(5))));
		
		% Fit the center.
		% #19 is the chisqperdof.
		[temp_guess, output_wconst(i-tmin+1, 19), succ_code] = fminsearch(chisqfunc, guess, optimset('TolX', 1e-15, 'TolFun', 1e-7, 'MaxFunEvals', 100000, 'MaxIter', 100000)); 
		if (succ_code == 0) % central failed
			disp(strcat('Failed on ', num2str(t1)));
			continue;
		end
		
		temp_guess(2) = abs(temp_guess(2));
		temp_guess(4) = abs(temp_guess(4));
		
		if (output_wconst(i-tmin+1,19) > 1e-20 && succ_code == 1)
			output_wconst(i-tmin+1, 1) = t1;
			output_wconst(i-tmin+1, 2) = t2;
			output_wconst(i-tmin+1, 3) = temp_guess(1);
			output_wconst(i-tmin+1, 5) = abs(temp_guess(2));
			output_wconst(i-tmin+1, 11) = temp_guess(3);
			output_wconst(i-tmin+1, 13) = abs(temp_guess(4));
			output_wconst(i-tmin+1, 7) = temp_guess(5);
			output_wconst(i-tmin+1, 22) = temp_guess(1)*cosh(temp_guess(2)*nt/2);
			
		end
		
		output_wconst(i-tmin+1, 20) = 1.0-chi2cdf(output_wconst(i-tmin+1,19)*(numel(yval)-numel(temp_guess)), numel(yval)-numel(temp_guess));
		
		% Save the coefficients for a jackknife.
		coefficient_center = temp_guess;
		% New 01-26-2015: Keep amp*cosh(mass(Nt/2))
		coefficient_center(6) = temp_guess(1)*cosh(temp_guess(2)*nt/2);
		coefficient_blocks = zeros(num_blocks, numel(temp_guess)+1);
		guess = temp_guess; % for next pass. 
		
		flag = 1;
		tic
		for b=1:num_blocks
			t1 = i;
			if (fold == 1)
				t2 = parse_Nt/2;
			else
				t2 = parse_Nt-i;
			end
			
			xval = t1:t2;
			yval = correlator_jack((t1+1):(t2+1), b)';
			ycorr = correlator_cov_mat((t1+1):(t2+1),(t1+1):(t2+1));
			%output_wconst(i-tmin+1, 21) = cond(ycorr);
			
			chisqfunc = @(x)(1.0/(size(yval,2)-5)*(yval(:)-x(1)*cosh(x(2)*(nt/2-xval(:)))-x(3)*(1-2*mod(xval(:),2)).*cosh(x(4)*(nt/2-xval(:)))-x(5))'*(ycorr\(yval(:)-x(1)*cosh(x(2)*(nt/2-xval(:)))-x(3)*(1-2*mod(xval(:),2)).*cosh(x(4)*(nt/2-xval(:)))-x(5))));
			
			% Fit the center.
			% #19 is the chisqperdof.
			[temp_guess, temp, succ_code] = fminsearch(chisqfunc, coefficient_center, optimset('TolX', 1e-15, 'TolFun', 1e-7, 'MaxFunEvals', 100000, 'MaxIter', 100000)); 
			if (succ_code == 0)
				flag = 0;
				continue;
			end
			
			temp_guess(2) = abs(temp_guess(2));
			temp_guess(4) = abs(temp_guess(4));
			coefficient_blocks(b,1:5) = temp_guess(1:5); 
			coefficient_blocks(b,6) = temp_guess(1)*cosh(temp_guess(2)*nt/2);
			
		end % b
		toc
		
		if flag == 0
			disp(strcat('Failed on ', num2str(t1)));
			continue;
		end
		coefficient_rep = repmat(coefficient_center, [num_blocks, 1]);
		coefficient_err = sqrt(sum((coefficient_blocks-coefficient_rep).^2, 1).*(num_blocks-1)./num_blocks);
		
		output_wconst(i-tmin+1, 4) = coefficient_err(1);
		output_wconst(i-tmin+1, 6) = coefficient_err(2);
		output_wconst(i-tmin+1, 12) = coefficient_err(3);
		output_wconst(i-tmin+1, 14) = coefficient_err(4);
		output_wconst(i-tmin+1, 8) = coefficient_err(5);
		output_wconst(i-tmin+1, 23) = coefficient_err(6);
		
	end
	
	
	% aaaaaand hopefully we're done?
	output_wconst2 = zeros(tmax, 23);
	output_wconst2(tmin:tmax, :) = output_wconst;
	
	full_fname = strcat(fname, '/spectrum2/fits/fit_new.', state, '_zcen');
    save(full_fname, 'output_wconst2', '-ascii');
    
    % Also just save fit masses (and fpi, in the case of fpi. Man fpi is
    % annoying.
    
    % Cosh state.
    raw_output = zeros(size(output_wconst2,1), 3);
    for j=1:(size(output_wconst2,1))
        raw_output(j,1) = real(output_wconst2(j,1));
        raw_output(j,2) = real(output_wconst2(j, 5));
        raw_output(j,3) = real(output_wconst2(j, 6));
    end
    full_fname = strcat(fname, '/spectrum2/fits/fit.', state, '_zcen', num2str(1));
    save(full_fname, 'raw_output', '-ascii');
	
	% Oscil state.
    raw_output = zeros(size(output_wconst2,1), 3);
    for j=1:(size(output_wconst2,1))
        raw_output(j,1) = real(output_wconst2(j,1));
        raw_output(j,2) = real(output_wconst2(j, 13));
        raw_output(j,3) = real(output_wconst2(j, 14));
    end
    full_fname = strcat(fname, '/spectrum2/fits/fit.', state, '_zcen', num2str(2));
    save(full_fname, 'raw_output', '-ascii');
	
	% Constant.
    raw_output = zeros(size(output_wconst2,1), 3);
    for j=1:(size(output_wconst2,1))
        raw_output(j,1) = real(output_wconst2(j,1));
        raw_output(j,2) = real(output_wconst2(j, 7));
        raw_output(j,3) = real(output_wconst2(j, 8));
    end
    full_fname = strcat(fname, '/spectrum2/fits/fit.', state, '_zcen', num2str(3));
    save(full_fname, 'raw_output', '-ascii');
end
