% Block data along the dimension specified.
function [blocks, num_blocks] = block_data(in_data, dim, blocksize)
	
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
	num_blocks = floor(num_data/blocksize);
	
	% Get the overflow!
	overflow = num_data - blocksize*num_blocks;
	front_overflow = floor(overflow/2);
	end_overflow = overflow-front_overflow;
	
	% Iterate over all blocks.
	other_data_size = numel(in_data_perm)/num_data;
	
	% Prepare space!
	the_blocks = zeros(num_blocks, other_data_size);
	
	% First block!
	first_set = in_data_perm(1:(blocksize+front_overflow), :);
	the_blocks(1,:) = mean(first_set, 1);
	
	% Last block!
	last_set = in_data_perm((num_data-blocksize-end_overflow+1):num_data,:);
	the_blocks(num_blocks,:) = mean(last_set, 1);
	
	% Trim start and end.
	%in_data_perm_orig = in_data_perm;
	in_data_perm((num_data-blocksize-end_overflow+1):num_data,:) = [];
	in_data_perm(1:(blocksize+front_overflow),:) = [];
	
	% Average over all other bins.
	the_blocks(2:(num_blocks-1),:) = reshape(mean(reshape(in_data_perm, [blocksize (num_blocks-2) other_data_size]), 1), [(num_blocks-2) other_data_size]);
	
	%meh = the_blocks;
	
	% And the rest!
	%for b=2:(num_blocks-1)
	%	a_range = ((b-1)*blocksize+1+front_overflow):(b*blocksize+front_overflow);
	%	middle_set = in_data_perm_orig(a_range,:);
	%	the_blocks(b,:) = mean(middle_set,1);
	%end
	
	% And we've blocked! Rearrange the data.
	unfolded_blocks = reshape(the_blocks, [num_blocks other_sizes]);
	
	% And fix the permutation.
	blocks = permute(unfolded_blocks, [2:dim 1 (dim+1):num_dimensions]);
	%blocks = permute(unfolded_blocks, [dim 1:(dim-1) (dim+1):num_dimensions]);
	
end