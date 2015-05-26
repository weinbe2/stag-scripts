function [ roots, roots_err, masses, masses_err, amps, amps_err ] = effective_mass_err(fname, state, parse_Nt, K, N, C, fold, do_err, blocksize, to_file)
    % Try this multi-state effective mass idea.
	
    % Load the data!
	connected = load_correlator(fname, state, parse_Nt);
    
    % Check if it's a Baryon. This only effects folding.
    is_baryon = 0;
    if (strcmp(state,'nu') || strcmp(state,'de'))
        is_baryon = 1;
    end
	
	% Fold it!
	if (fold == 1)
		connected = fold_data(connected, is_baryon);
	end

	% Block data!
	[connected_blocks, num_blocks] = block_data(connected, 2, blocksize);
	
	% Get jackknife blocks and covariance matrix.
	connected_sum = mean(connected_blocks, 2);
	[connected_jack, connected_cov_mat, connected_err] = jackknife_from_blocks(connected_blocks);
	
	% And we're off!

	% Get a central value!
	[ masses_center, roots_center, amps_center ] = effective_mass_utility(connected_sum, parse_Nt, K, N, C);
	
	if (do_err == 1)
	
		% Prepare for blocks!
		masses_jack = zeros([size(masses_center) num_blocks]);
		roots_jack = zeros([size(roots_center) num_blocks]);
		amps_jack = zeros([size(amps_center) num_blocks]);
		
		for b=1:num_blocks
			[ masses_jack(:,:,b), roots_jack(:,:,b), amps_jack(:,:,b)] = effective_mass_utility(connected_jack(:,b), parse_Nt, K, N, C);
		end
		
		% Get a jackknife error!
		masses_rep = repmat(masses_center, [1 1 num_blocks]);
		roots_rep = repmat(roots_center, [1 1 num_blocks]);
		amps_rep = repmat(amps_center, [1 1 num_blocks]);
		
		masses_err = sqrt(sum((masses_rep-masses_jack).^2, 3).*(num_blocks-1)/(num_blocks));
		roots_err = sqrt(sum((roots_rep-roots_jack).^2, 3).*(num_blocks-1)/(num_blocks));
		amps_err = sqrt(sum((amps_rep-amps_jack).^2, 3).*(num_blocks-1)/(num_blocks));
		
		if (to_file == 1)
		   for i=1:K
			  raw_output = zeros(size(masses_err,1), 3);
			  for j=1:(size(masses_err,1))
				 raw_output(j,1) = j-1+((N-1)/2);
				 raw_output(j,2) = real(masses_center(j, i));
				 raw_output(j,3) = real(masses_err(j, i));
			  end
			  full_fname = strcat(fname, '/spectrum2/effmass/effmass.', state, num2str(i));
        
			  save(full_fname, 'raw_output', '-ascii');
			   
		   end
			
		end
    else
        roots_err = [];
        masses_err = [];
        amps_err = [];
    end
	
	roots = roots_center;
	masses = masses_center;
	amps = amps_center; 

end
