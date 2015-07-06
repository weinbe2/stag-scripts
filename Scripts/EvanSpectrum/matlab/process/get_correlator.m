function [connected_sum, connected_jack, connected_cov_mat, connected_err, num_blocks, connected_jack_single] = get_correlator(fname, state, parse_Nt, parse_Ns, fl_flavor, binsize, num_elim, subset_beg, subset_end)
	% Loads the central value and jackknife blocks in one go, independent of connected or disconnected. 
	% fname is relative to the current path.
	% state is 'ps', 'ps2', etc.
	% parse_Nt is the temporal dimension.
	% parse_Ns is the spatial dimension.
	% binsize is the size to bin the data before doing
	%   single elimination jackknife. 
	% num_elim is optional. If it's not set, it does a single-elimination.
	%   Otherwise, it does a 'num_elim' elimination.
	% subset_beg and subset_end are optional. If not set, we use all of the loaded data.
	%   Otherwise, the percent of the data between subset_beg (default 0) and subset_end
	%   (default 1) is kept.
	% Note: var-covar always comes from single-elimination jack.
	
	% The output 'connected_jack_single' is optional.
	% It'll give the jackknife blocks from a single elimination
	% if so desired. 
	
	if (~exist('num_elim', 'var'))
		num_elim = 1; % single eliminate
	end
	
	want_single = 0;
	if (nargout == 6) % want single elim!
		want_single = 1;
	end
	
	if (~exist('subset_beg', 'var') && ~exist('subset_end', 'var'))
        subset_beg = 0;
        subset_end = 1;
    end
    if (subset_beg < 0)
        subset_beg = 0;
    end
    if (subset_end > 1)
        subset_end = 1;
    end
	
	
	if (strcmp(state, 'dc_stoch') || strcmp(state, 'sg_stoch'))
        % Get that disconnected state!
        
		number_bl = 6; 
		
		% Load the pbp values. Recall pbp is (parse_Nt, number_bl_in, num_data);
		[pbp config_nums_pbp] = load_pbppart(fname, parse_Nt, parse_Ns, number_bl);
		num_data = size(pbp, 3);
		
		% Remove part of the data if we need to.
		if (~(subset_beg == 0 && subset_end == 1))
			num_data = size(pbp, 3);
			num_beg = floor(num_data*subset_beg)+1;
			num_end = floor(num_data*subset_end);
			pbp = pbp(:,:,num_beg:num_end);
		end
		
        % Rescale the ops so we don't need to multiply D by anything later.
		pbp = pbp.*sqrt(fl_flavor);
		
		% Build all the non-vev subtracted correlators.
		% (parse_Nt, number_bl_in, number_bl_in, num_data);
		disc = build_vev_correlator(pbp, 1); % fold it.
		
		% Now bin it!
		[pbp_blocks num_blocks] = block_data(pbp, 3, binsize);
		[disc_blocks num_blocks] = block_data(disc, 4, binsize);
		
		% Get central values.
		pbp_sum = mean(pbp_blocks, 3);
		disc_sum = mean(disc_blocks, 4);
		
		% Multi- elim jackknife.
		pbp_jack = jackknife_bins(pbp_blocks, 3, num_elim);
		disc_jack = jackknife_bins(disc_blocks, 4, num_elim);
		num_blocks = size(pbp_jack, 3);
		
		% Build the disconnected correlators.
		vev_center = mean(pbp_sum, 1);
		vev_sub = repmat(reshape(parse_Nt.*((vev_center')*vev_center),1,number_bl,number_bl), [parse_Nt, 1, 1]);
		disc_sum = disc_sum - vev_sub;
		for i=1:number_bl
			disc_sum(:,i,i,:) = 0;
		end
		disc_sum = sum(sum(disc_sum,2),3)/(number_bl*(number_bl-1));
		
		% Next on the jackknife!
		vev_center = mean(pbp_jack, 1);
		vev_1 = repmat(reshape(vev_center,1, number_bl, 1, num_blocks), [parse_Nt, 1, number_bl, 1]);
		vev_2 = repmat(reshape(vev_center,1, 1, number_bl, num_blocks), [parse_Nt, number_bl, 1, 1]);
		disc_jack = disc_jack - vev_1.*vev_2.*parse_Nt;
		for i=1:number_bl
			disc_jack(:,i,i,:) = 0;
		end
		disc_jack2 = sum(sum(disc_jack,2),3)/(number_bl*(number_bl-1));
		
		disc_jack = zeros(parse_Nt, num_blocks);
		disc_jack(:,:) = disc_jack2(:,1,1,:);
		
		% Grab the single elim if wanted.
		if (want_single == 1)
			% Single elim
			pbp_jack_sing = jackknife_bins(pbp_blocks, 3, 1);
			disc_jack_sing = jackknife_bins(disc_blocks, 4, 1);
			num_blocks_sing = size(pbp_jack_sing, 3);
			
			% Next on the jackknife!
			vev_center = mean(pbp_jack_sing, 1);
			vev_1 = repmat(reshape(vev_center,1, number_bl, 1, num_blocks_sing), [parse_Nt, 1, number_bl, 1]);
			vev_2 = repmat(reshape(vev_center,1, 1, number_bl, num_blocks_sing), [parse_Nt, number_bl, 1, 1]);
			disc_jack_sing = disc_jack_sing - vev_1.*vev_2.*parse_Nt;
			for i=1:number_bl
				disc_jack_sing(:,i,i,:) = 0;
			end
			disc_jack2 = sum(sum(disc_jack_sing,2),3)/(number_bl*(number_bl-1));
			
			disc_jack_sing = zeros(parse_Nt, num_blocks_sing);
			disc_jack_sing(:,:) = disc_jack2(:,1,1,:);
		end
		
		% Grab the connected piece if we need to!
		if (strcmp(state, 'sg_stoch'))
			
			% Load it, fold it...
			connected = load_correlator(fname, 'sc_stoch', parse_Nt);
			%connected = fold_data(connected, is_baryon);
			
			% Remove part of the data if we need to.
			if (~(subset_beg == 0 && subset_end == 1))
				num_data = size(connected, 2);
				num_beg = floor(num_data*subset_beg)+1;
				num_end = floor(num_data*subset_end);
				connected = connected(:,:,num_beg:num_end);
			end
			
			% Bin it!
			connected_blocks = block_data(connected, 2, binsize);
			
			% Sum it, jack it.
			connected_sum = mean(connected_blocks, 2);
			connected_jack = jackknife_bins(connected_blocks, 2, num_elim);
			
			% Grab the single elim if wanted.
			if (want_single == 1)
				connected_jack_single = jackknife_bins(connected_blocks, 2, 1);
			end
			
			% And modify disc.
			disc_sum = disc_sum - connected_sum;
			disc_jack = disc_jack - connected_jack;
			
			if (want_single == 1)
				disc_jack_sing = disc_jack_sing - connected_jack_single;
			end
		end
		
		% Rename!
		connected_sum = disc_sum;
		connected_jack = disc_jack; 
		
		if (want_single == 1)
			connected_jack_single = disc_jack_sing;
		end
		
    else

        % Load the correlator, run autocorrelation, block it.
        connected = load_correlator(fname, state, parse_Nt);

		% Remove part of the data if we need to.
		if (~(subset_beg == 0 && subset_end == 1))
			num_data = size(connected, 2);
			num_beg = floor(num_data*subset_beg)+1;
			num_end = floor(num_data*subset_end);
			connected = connected(:,:,num_beg:num_end);
		end
		
        % Fold after the fact. 
		%{
        if (strcmp(state, 'sc_stoch'))
            connected = -1*fold_data(connected, is_baryon);
        else
            connected = fold_data(connected, is_baryon);
        end
		%}
		
		if (strcmp(state, 'sc_stoch'))
			connected = -1*connected;
		end

        % Block it.
        connected_blocks = block_data(connected, 2, binsize);

        % Sum, blocks, errors, etc.
        connected_sum = mean(connected_blocks, 2);
        connected_jack = jackknife_bins(connected_blocks, 2, num_elim);
        %[connected_P_cov_mat, connected_P_err] = errors_jackknife(connected_P_sum, connected_P_jack);
		
		num_blocks = size(connected_jack, 2);
		
		% Single elim.
		if (want_single == 1)
			connected_jack_single = jackknife_bins(connected_blocks, 2, 1);
		end
    end

	[connected_cov_mat, connected_err] = errors_jackknife(connected_sum, connected_jack);
	
end