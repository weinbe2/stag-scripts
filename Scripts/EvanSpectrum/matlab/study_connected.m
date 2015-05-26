% Run analysis on a state. blocksize, tmin_min and tmin_max are optional arguments!
% If the optional arguments aren't supplied, it looks for a fitparams file.
% 08-20-2014: Load information from info.dat.
% 01-12-2015: Added 'noerrors' flag, which much be set to 0 or 1 if one wants to specify maxes and such.
%             Default behavior of noerrors is 0.
function outform = study_connected(fname, state, noerrors, blockval, tmin_min, tmin_max, tmax)
	% Load ensemble info.
	%[fl_l, fl_h, parse_Ns, parse_Nt, beta, m_l, m_h] = load_info_file(fname);
    [~, ~, parse_Ns, parse_Nt, ~, m_l, ~] = load_info_file(fname);
	
	min_search = 0; max_search = 0; tmax_search = 0; blocksize = 0;
	check_errs = 0;
	
	if (exist('noerrors', 'var'))
		check_errs = noerrors;
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
  
	% Learn about our measurement from the state name.
	is_baryon = 0;
	if (strcmp(state, 'nu') || strcmp(state, 'de'))
		is_baryon = 1;
	end
	
	is_fpi = 0;
	if (strcmp(state, 'ps2'))
		is_fpi = 1;
    end
	
	is_oscil = 1;
	% Hack in wflow states
	if (strcmp(state, 'ps') || strcmp(state, 'ps2') || strcmp(state, 'i5') || strcmp(state, 'ps_wt0_10') || strcmp(state, 'ps2_wt0_10') || strcmp(state, 'i5_wt0_10') || strcmp(state, 'ps_wt0_30') || strcmp(state, 'ps2_wt0_30') || strcmp(state, 'i5_wt0_30'))
		is_oscil = 0;
	end
	
	% Let's save a sum corr value while we're here, shall we?
	if (~strcmp(state, 'sc_stoch')) % sc_stoch gets saved elsewhere.
		connected = load_correlator(fname, state, parse_Nt);
	
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

	end
	
	% Spit out some effective mass data.
	effective_mass_err(fname, state, parse_Nt, 1+is_oscil, 2+2*is_oscil, 0, 1, 1, blocksize, 1);
	
	% Spit out autocorrelation info. As of 2015-01-12, we now save this to file.
	wall_corr = load_correlator(fname, state, parse_Nt);
	
		% Standard things we output.
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
    
    
	
    for i=1:count
		if (is_baryon == 0 && is_fpi == 0)
			output_ps(i,:) = quick_results(fname, state, parse_Ns, parse_Nt, blocksize, i, tmax_search, 0, 1, check_errs);
		elseif (is_baryon == 1)
			output_ps(i,:) = quick_results_baryon(fname, state, parse_Ns, parse_Nt, blocksize, i, tmax_search, 0, 1, check_errs); % We can fold.
		else % fpi
			output_ps(i,:) = quick_results_fpi(fname, state, parse_Ns, parse_Nt, blocksize, i, tmax_search, 0, 1, mass_l, check_errs);
		end
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
