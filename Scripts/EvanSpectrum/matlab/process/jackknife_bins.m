% Do single elimination jackknife along the dim specified.
function [out_jack] = jackknife_bins(in_data, dim)
	
	num_dimensions = numel(size(in_data));
	dims = 1:num_dimensions;
	dims(dim) = [];
	rearrange = [dim dims];
	
	in_data_perm = permute(in_data, rearrange);
	
	% Get the size of the data in that direction.
	num_data = size(in_data_perm, 1);
	other_sizes = size(in_data_perm);
	other_sizes(1) = [];
	
	% Iterate over all blocks.
	other_data_size = numel(in_data_perm)/num_data;
	
	% Prepare space!
	%in_data_perm = reshape(in_data_perm, [num_data, other_data_size]);
	
	the_blocks = zeros(num_data, other_data_size);
	
	for b=1:num_data
		file_conn_copy = in_data_perm;
		file_conn_copy(b,:) = [];
		the_blocks(b,:) = mean(file_conn_copy,1);
	end
	
	% And we've jackknifed! Rearrange the data.
	unfolded_blocks = reshape(the_blocks, [num_data other_sizes]);
	
	% And fix the permutation.
	blocks = permute(unfolded_blocks, [2:dim 1 (dim+1):num_dimensions]);
	
	out_jack = blocks;
	
end