function [ root, masses, amps ] = effective_mass(fname, state, parse_Nt, K, N, C)

    % Try this multi-state effective mass idea.
	connected = load_correlator(fname, state, parse_Nt);
	
	% Check if it's a Baryon. This only effects folding.
    is_baryon = 0;
    if (strcmp(state,'nu') || strcmp(state,'de'))
        is_baryon = 1;
    end

    % Default: We fold it.
	connected = fold_data(connected, is_baryon);

    connected_sum = mean(connected, 2);

	[ masses, root, amps ] = effective_mass_utility(connected_sum, parse_Nt, K, N, C);
	
	
end
