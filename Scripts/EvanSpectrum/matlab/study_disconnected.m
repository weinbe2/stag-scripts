% Run disconnected analysis. blocksize, tmin_min and tmin_max are optional arguments!
% If the optional arguments aren't supplied, it looks for a fitparams file.
% 2015-01-13: New dc_three, sg_three to pass off to test_three_fits function.
function outform = study_disconnected(fname, state, number_bl, noerrors, blockval, tmin_min, tmin_max, tmax)
	% Load ensemble info.
	%[fl_l, fl_h, parse_Ns, parse_Nt, beta, m_l, m_h] = load_info_file(fname);
    [fl_l, fl_h, parse_Ns, parse_Nt, beta, m_l, m_h] = load_info_file(fname);
	fl_flavor = fl_l/4; % Staggered counting.
    
	check_errs = 0;
	
	if (exist('noerrors', 'var'))
		check_errs = noerrors;
	end
	
	% State:
	% dc_stoch -- just disconnected
	% dc_stoch_oscil -- disconnected, but fit cosh+oscil.
	% dc_stoch_ppp -- disconnected, positive parity project.
	% sg_stoch -- Nf/4 D - C, with oscil.
	
	% New things:
	% dc_stoch_three -- do a cosh+oscil, cosh+oscil+const, finite diff cosh+oscil fit on disconnected.
	% sg_stoch_three -- do a cosh+oscil, cosh+oscil+const, finite diff cosh+oscil fit on N_f/4 - C
	
	min_search = 0; max_search = 0; tmax_search = 0; blocksize = 0;
	temp_state = state;
	
	if (strcmp(state, 'dc_stoch_oscil') || strcmp(state, 'dc_stoch_ppp') || strcmp(state, 'dc_stoch_three'))
		state = 'dc_stoch';
	end
	
	if (strcmp(state, 'sg_stoch_three'))
		state = 'sg_stoch';
	end
	
	% If we didn't get range arguments, load disconnected...
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
	state = temp_state;
	
	% Offload to study_disconnected_three
	if (strcmp(state, 'sg_stoch_three') || strcmp(state, 'dc_stoch_three'))
		% We don't make study_disconnected_three load state parameters---it just gets it all.
		outform = study_disconnected_three(fname, state, number_bl, check_errs, blocksize, tmax_search);
	else % proceed as is.
	
		% For now, default to fitting dc_stoch.
		
		if (strcmp(state, 'dc_stoch_oscil') || strcmp(state, 'sg_stoch'))
			is_oscil = 1;
		else
			is_oscil = 0;
		end
		
		
		% Load up dc_stoch, while we're here.
		[disc_sum, disc_jack, disc_cov_mat, disc_err] = load_correlator_scalar(fname, fl_flavor, parse_Nt, parse_Ns, number_bl, blocksize);
		
		% Positive parity project if we need to.
		if (strcmp(state, 'dc_stoch_ppp'))
		
			disc_sum = ppp_data(disc_sum);
			disc_jack = ppp_data(disc_jack);
			
			% Recompute the cov mat, error.
			disc_rep = repmat(disc_sum, [1 size(disc_jack, 2)]);
			disc_cov_mat = zeros(parse_Nt);
			disc_err = zeros(parse_Nt, 1);
			for t1 = 1:parse_Nt
				for t2 = 1:parse_Nt
					disc_cov_mat(t1,t2) = sum((disc_rep(t1,:)-disc_jack(t1,:)).*(disc_rep(t2,:)-disc_jack(t2,:)),2).*(size(disc_jack,2)-1)./size(disc_jack,2);
					if t1 == t2
						disc_err(t1) = sqrt(disc_cov_mat(t1,t1));
					end
				end
			end
			clear('disc_rep');
			
		end
		
		% Subtract sigma, as well.
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
		
		% Count gets computed a little bit differently now.
		% Start at tmin = 1, go until tmax - dof-1.
		
		%       tmax          [dof]              [tmin + 1]
		count = tmax_search - (2 + 2*is_oscil);
		
		output_ps = zeros(count, 21);
		
		for i=1:count
			% Call the quick_results given blocks!
			% The first 0 refers to a non-oscillating fit.
			output_ps(i,:) = quick_results_blocks(disc_sum, disc_jack, disc_cov_mat, is_oscil, parse_Nt, i, tmax_search, 0, 1, check_errs);
			
		end
		
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
		
		outform = output_ps;
	end

end
