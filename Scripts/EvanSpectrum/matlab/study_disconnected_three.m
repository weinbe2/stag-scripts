% 2015-01-05 Test cosh+oscil, cosh+oscil+const, cosh+oscil diff fits. 
%function study_disconnected_three(fname, blockval, diag, fold)
function outform = study_disconnected_three(fname, state, number_bl, check_errs, blocksize, tmax_search)
	% Check if we're actually looking at a disconnected state.
	if (~(strcmp(state, 'dc_stoch_three') || strcmp(state, 'sg_stoch_three')))
		disp('Wrong type of state!')
		return;
    end
	
	if (strcmp(state, 'sg_stoch_three'))
		csub_flag = 1; % subtract the connected piece,
	else
		csub_flag = 0; % or not, for just D.
	end
	
	% Load ensemble info.
	%[fl_l, fl_h, parse_Ns, parse_Nt, beta, m_l, m_h] = load_info_file(fname);
    [fl_l, fl_h, parse_Ns, parse_Nt, beta, m_l, m_h] = load_info_file(fname);
	fl_flavor = fl_l/4; % Staggered counting.
    
	fold = 1; % fold it
	diag = 0; % don't use a diagonal cov matrix.
	mass_l = m_l;
	
	
	
	% Load both D and C.
	[disc_sum, disc_jack, disc_cov_mat, disc_err] = load_correlator_scalar(fname, fl_flavor, parse_Nt, parse_Ns, number_bl, blocksize);
	
	% Take derivatives!
	disc_deriv_sum = zeros(size(disc_sum));
	disc_deriv_jack = zeros(size(disc_jack));
	for t=1:(parse_Nt)
		disc_deriv_sum(t,:) = disc_sum(mod(t,parse_Nt)+1,:)-disc_sum(t,:);
		disc_deriv_jack(t,:) = disc_jack(mod(t,parse_Nt)+1,:)-disc_jack(t,:);
	end
	
	% Fake it if we don't need it.
	connected = load_correlator(fname, 'sc_stoch', parse_Nt);
	if (csub_flag == 0)
		connected(:,:) = 0;
	end

	% NEW: Take the derivative a la Claudio's method.
	% Times are now interpretted to be halfway between original
	% timeslices.
    
    % Do we fold?
    if fold == 1
		% We fold it.
		connected = fold_data(connected, 0); % 0 b/c not a baryon.
    end
	
	% We hold onto some extra data at the start.
	connected_deriv = zeros(size(connected));
	
	for t=1:(parse_Nt)
		connected_deriv(t,:) = connected(mod(t,parse_Nt)+1,:)-connected(t,:);
	end
	
	% Carry on with blocking and whatnot.
	
	% Block data!
	[connected_blocks, num_blocks] = block_data(connected, 2, blocksize);
	[connected_deriv_blocks, num_blocks] = block_data(connected_deriv, 2, blocksize);
	
	% Get jackknife blocks and covariance matrix.
	connected_sum = mean(connected_blocks, 2);
	connected_deriv_sum = mean(connected_deriv_blocks, 2);
	
	
	[connected_jack, connected_cov_mat, connected_err] = jackknife_from_blocks(connected_blocks);
	[connected_deriv_jack, connected_deriv_cov_mat, connected_deriv_err] = jackknife_from_blocks(connected_deriv_blocks);
	
	% We still need to form sigma, as well as the covariance matrix.
    sigma_sum = disc_sum - connected_sum;
    sigma_rep = repmat(sigma_sum, [1 num_blocks]);
    sigma_jack = disc_jack - connected_jack;
    
    sigma_deriv_sum = disc_deriv_sum - connected_deriv_sum;
    sigma_deriv_rep = repmat(sigma_deriv_sum, [1 num_blocks]);
    sigma_deriv_jack = disc_deriv_jack - connected_deriv_jack;
    
    sigma_cov_mat = zeros(parse_Nt);
    sigma_err = zeros(parse_Nt, 1);
    sigma_deriv_cov_mat = zeros(parse_Nt);
    sigma_deriv_err = zeros(parse_Nt, 1);
    for t1 = 1:parse_Nt
        for t2 = 1:parse_Nt
            sigma_cov_mat(t1,t2) = sum((sigma_rep(t1,:)-sigma_jack(t1,:)).*(sigma_rep(t2,:)-sigma_jack(t2,:)),2).*(num_blocks-1)./num_blocks;
            sigma_deriv_cov_mat(t1,t2) = sum((sigma_deriv_rep(t1,:)-sigma_deriv_jack(t1,:)).*(sigma_deriv_rep(t2,:)-sigma_deriv_jack(t2,:)),2).*(num_blocks-1)./num_blocks;
            if t1 == t2
                sigma_err(t1) = sqrt(sigma_cov_mat(t1,t1));
                sigma_deriv_err(t1) = sqrt(sigma_deriv_cov_mat(t1,t1));
            end
        end
    end
	
    % Save the derivative.
	data = zeros((parse_Nt/2)+1,3);
        for i=1:((parse_Nt/2)+1)
            data(i,1) = i-1;
            data(i,2) = sigma_deriv_sum(i);
            data(i,3) = sigma_deriv_err(i);
        end

        full_fname = strcat(fname, '/spectrum2/sum/sum.', state);

        save(full_fname, 'data', '-ascii');
    return; 

    % Alright, we now have everything.
    
    % Count gets computed a little bit differently now.
    % Start at tmin = 1, go until tmax - dof-1.
    
    %       extra *1 is because there's an oscillating term, +1 is for the constant.
	
    count = tmax_search - (2 + 2*1 +1);
    
	output_orig = zeros(count, 21);
	output_deriv = zeros(count, 25);
	output_wconst = zeros(count, 21);
    
    % We need to save some values at some point!
	
    for i=1:count
		% Call the quick_results given blocks!
		% The first 0 refers to a non-oscillating fit.
		output_deriv(i,:) = quick_results_blocks_deriv(sigma_deriv_sum, sigma_deriv_jack, sigma_deriv_cov_mat, 13, parse_Nt, i, tmax_search, diag, fold, check_errs); 
        output_wconst(i,:) = quick_results_blocks(sigma_sum, sigma_jack, sigma_cov_mat, 2, parse_Nt, i, tmax_search, diag, fold, check_errs); 
        output_orig(i,:) = quick_results_blocks(sigma_sum, sigma_jack, sigma_cov_mat, 1, parse_Nt, i, tmax_search, diag, fold, check_errs); 
    end
	
	% Fix the convention on deriv.
	output_deriv = output_deriv(:, [1 2 22 23 5 6 7 8 9 10 24 25 13 14 15 16 17 18 19 20 21 3 4 11 12]);
    
    % Save things!
    if (csub_flag == 0)
		full_fname = strcat(fname, '/spectrum2/fits/fit_new.dc_stoch_oscil');
		save(full_fname, 'output_orig', '-ascii');
		
		full_fname = strcat(fname, '/spectrum2/fits/fit_new.dc_stoch_deriv');
		save(full_fname, 'output_deriv', '-ascii');
		
		full_fname = strcat(fname, '/spectrum2/fits/fit_new.dc_stoch_wconst');
		save(full_fname, 'output_wconst', '-ascii');
		
	else
		full_fname = strcat(fname, '/spectrum2/fits/fit_new.sg_stoch');
		save(full_fname, 'output_orig', '-ascii');
		
		full_fname = strcat(fname, '/spectrum2/fits/fit_new.sg_stoch_deriv');
		save(full_fname, 'output_deriv', '-ascii');
		
		full_fname = strcat(fname, '/spectrum2/fits/fit_new.sg_stoch_wconst');
		save(full_fname, 'output_wconst', '-ascii');
	end
	
	% Also just save fit masses.
		
	% First, the original fit.
	
	% Cosh state.
	raw_output = zeros(size(output_orig,1), 3);
	for j=1:(size(output_orig,1))
		raw_output(j,1) = real(output_orig(j,1));
		raw_output(j,2) = real(output_orig(j, 5));
		raw_output(j,3) = real(output_orig(j, 6));
	end
	if (csub_flag == 0)
		full_fname = strcat(fname, '/spectrum2/fits/fit.dc_stoch_oscil1');
	else
		full_fname = strcat(fname, '/spectrum2/fits/fit.sg_stoch1');
	end
	save(full_fname, 'raw_output', '-ascii');

	% Oscil state, if it exists.
	raw_output = zeros(size(output_orig,1), 3);
	for j=1:(size(output_orig,1))
		raw_output(j,1) = real(output_orig(j,1));
		raw_output(j,2) = real(output_orig(j, 13));
		raw_output(j,3) = real(output_orig(j, 14));
	end
	if (csub_flag == 0)
		full_fname = strcat(fname, '/spectrum2/fits/fit.dc_stoch_oscil2');
	else
		full_fname = strcat(fname, '/spectrum2/fits/fit.sg_stoch2');
	end
	save(full_fname, 'raw_output', '-ascii');
     
	% Next, the derivative.
	
	% Cosh state.
	raw_output = zeros(size(output_deriv,1), 3);
	for j=1:(size(output_deriv,1))
		raw_output(j,1) = real(output_deriv(j,1));
		raw_output(j,2) = real(output_deriv(j, 5));
		raw_output(j,3) = real(output_deriv(j, 6));
	end
	if (csub_flag == 0)
		full_fname = strcat(fname, '/spectrum2/fits/fit.dc_stoch_deriv1');
	else
		full_fname = strcat(fname, '/spectrum2/fits/fit.sg_stoch_deriv1');
	end
	save(full_fname, 'raw_output', '-ascii');

	% Oscil state, if it exists.
	raw_output = zeros(size(output_deriv,1), 3);
	for j=1:(size(output_deriv,1))
		raw_output(j,1) = real(output_deriv(j,1));
		raw_output(j,2) = real(output_deriv(j, 13));
		raw_output(j,3) = real(output_deriv(j, 14));
	end
	if (csub_flag == 0)
		full_fname = strcat(fname, '/spectrum2/fits/fit.dc_stoch_deriv2');
	else
		full_fname = strcat(fname, '/spectrum2/fits/fit.sg_stoch_deriv2');
	end
	save(full_fname, 'raw_output', '-ascii');
	
	% Finally, the extra constant.
	
	% Cosh state.
	raw_output = zeros(size(output_wconst,1), 3);
	for j=1:(size(output_wconst,1))
		raw_output(j,1) = real(output_wconst(j,1));
		raw_output(j,2) = real(output_wconst(j, 5));
		raw_output(j,3) = real(output_wconst(j, 6));
	end
	if (csub_flag == 0)
		full_fname = strcat(fname, '/spectrum2/fits/fit.dc_stoch_wconst1');
	else
		full_fname = strcat(fname, '/spectrum2/fits/fit.sg_stoch_wconst1');
	end
	save(full_fname, 'raw_output', '-ascii');

	% Oscil state, if it exists.
	raw_output = zeros(size(output_wconst,1), 3);
	for j=1:(size(output_wconst,1))
		raw_output(j,1) = real(output_wconst(j,1));
		raw_output(j,2) = real(output_wconst(j, 13));
		raw_output(j,3) = real(output_wconst(j, 14));
	end
	if (csub_flag == 0)
		full_fname = strcat(fname, '/spectrum2/fits/fit.dc_stoch_wconst2');
	else
		full_fname = strcat(fname, '/spectrum2/fits/fit.sg_stoch_wconst2');
	end
	save(full_fname, 'raw_output', '-ascii');
	
	% Dat constant!
	raw_output = zeros(size(output_wconst,1), 3);
	for j=1:(size(output_wconst,1))
		raw_output(j,1) = real(output_wconst(j,1));
		raw_output(j,2) = real(output_wconst(j, 7));
		raw_output(j,3) = real(output_wconst(j, 8));
	end
	if (csub_flag == 0)
		full_fname = strcat(fname, '/spectrum2/fits/fit.dc_stoch_wconst3');
	else
		full_fname = strcat(fname, '/spectrum2/fits/fit.sg_stoch_wconst3');
	end
	save(full_fname, 'raw_output', '-ascii');
    
	outform = output_wconst;
	
end
