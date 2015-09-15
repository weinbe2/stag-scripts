% 2015-09-15 Just save a correlator!
function build_correlator(fname, state, blockval)
	% Load the path with routines to load, bin, etc data.
	addpath('./process', '-end');
	

	% Load ensemble info.
	[fl_l, fl_h, parse_Ns, parse_Nt, beta, m_l, m_h] = load_info_file(fname);

    volume = parse_Ns^3;
    fl_flavor = fl_l/4; % Because staggered.
	
	blocksize = 0;
	fl_flavor = fl_l/4;
	
	% If we didn't get range arguments...
	if (~exist('blockval', 'var'))
		% See if a fitparams file exists.
		if (exist(fullfile(cd, strcat([fname '/spectrum2/fitparams/fitparam.' state])), 'file'))
			fitparam = importdata(strcat([fname '/spectrum2/fitparams/fitparam.' state]));
			blocksize = fitparam(5);
			
		else
			error('A fitparams file is expected but does not exist!');
		end
	else
		blocksize = blockval;
	end
	
	%%%%%%%%%%%%%%%%%%%%%%%%%
	% Learn based on state. %
	%%%%%%%%%%%%%%%%%%%%%%%%%
	
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
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Load and save correlator. %
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
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
    save(full_fname, 'output_autocorr_info', '-ascii', '-double');
		
	% End outputting autocorrelation info.
		
	% Fold the correlator, appropriately handling if it's a baryon.
    connected = fold_data(connected, is_baryon);

	% Block data!
	[connected_blocks, num_blocks] = block_data(connected, 2, blocksize);
		
	% Get jackknife blocks and covariance matrix.
	connected_sum = mean(connected_blocks, 2);
	connected_jack = jackknife_bins(connected_blocks, 2);
	[connected_cov_mat, connected_err] = errors_jackknife(connected_sum, connected_jack);
		
	% Save it!
    data = zeros((parse_Nt/2)+1,3);
	for i=1:((parse_Nt/2)+1)
		data(i,1) = i-1;
		data(i,2) = connected_sum(i);
		data(i,3) = connected_err(i);
	end
	full_fname = strcat(fname, '/spectrum2/sum/sum.', state);	
	save(full_fname, 'data', '-ascii', '-double');
	
	% Spit out some effective mass data.
	% This function computes and saves it.
	effective_mass_err(fname, state, parse_Nt, 1+is_oscil, 2+2*is_oscil, 0, 1, 1, blocksize, 1);
	

end