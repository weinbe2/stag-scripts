% Do single elimination jackknife along the dim specified. Num_elim is optional, and makes it a multi-elim jackknife.
function [out_jack] = jackknife_bins(in_data, dim, num_elim)
	
	if (~exist('num_elim', 'var'))
		num_elim = 1; % single eliminate
	end
	
	num_dimensions = numel(size(in_data));
	dims = 1:num_dimensions;
	dims(dim) = [];
	rearrange = [dim dims];
	
	in_data_perm = permute(in_data, rearrange);
	
	% Get the size of the data in that direction.
	num_data = size(in_data_perm, 1);
	other_sizes = size(in_data_perm);
	other_sizes(1) = [];
	
	% Get the number of blocks.
	num_blocks = floor(num_data/num_elim);
	
	% Get the overflow!
	overflow = num_data - num_elim*num_blocks;
	front_overflow = floor(overflow/2);
	end_overflow = overflow-front_overflow;
	
	% Iterate over all blocks.
	other_data_size = numel(in_data_perm)/num_data;
	
	% Prepare space!
	%in_data_perm = reshape(in_data_perm, [num_data, other_data_size]);
	
	the_blocks = zeros(num_blocks, other_data_size);
	
	% First block!
	file_conn_copy = in_data_perm;
	file_conn_copy(1:(num_elim+front_overflow),:) = [];
	the_blocks(1,:) = mean(file_conn_copy,1);
	
	% Last block!
	file_conn_copy = in_data_perm;
	file_conn_copy((num_data-num_elim-end_overflow+1):num_data,:) = [];
	the_blocks(num_blocks,:) = mean(file_conn_copy,1);
    
    % Get sum of first and last blocks, as well as all.
    %in_data_perm_orig = in_data_perm;
    in_data_perm = reshape(in_data_perm, [num_blocks other_data_size]);
    sum_all = sum(in_data_perm,1);
    %sum_first = sum(in_data_perm(1:(num_elim+front_overflow),:),1);
    %sum_last = sum(in_data_perm((num_data-num_elim-end_overflow+1):num_data,:),1);
    
    % Trim start and end.
	in_data_perm((num_data-num_elim-end_overflow+1):num_data,:) = [];
	in_data_perm(1:(num_elim+front_overflow),:) = [];
    
    % Get sum of all other bins.
	sum_other_bins = reshape(sum(reshape(in_data_perm, [num_elim (num_blocks-2) other_data_size]),1), [(num_blocks-2) other_data_size]);
    block_vals_tmp = (repmat(sum_all, [(num_blocks-2) 1]) - sum_other_bins)/(num_data-num_elim);
    the_blocks(2:(num_blocks-1),:) = block_vals_tmp;
    
    %meh = the_blocks;
    
	% And the rest!
	%for b=2:(num_blocks-1)
	%	file_conn_copy = in_data_perm_orig;
	%	file_conn_copy(((b-1)*num_elim+front_overflow+1):(b*num_elim+front_overflow),:) = [];
	%	the_blocks(b,:) = mean(file_conn_copy,1);
	%end
	
	% And we've jackknifed! Rearrange the data.
	unfolded_blocks = reshape(the_blocks, [num_blocks other_sizes]);
	
	% And fix the permutation.
	blocks = permute(unfolded_blocks, [2:dim 1 (dim+1):num_dimensions]);
	
	out_jack = blocks;
	
end