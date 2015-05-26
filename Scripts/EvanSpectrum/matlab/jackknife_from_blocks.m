function [out_jack, out_cov_mat, out_err] = jackknife_from_blocks(in_blocks)
	wall_blocks = in_blocks; % Eh, too lazy to rename.
	parse_Nt = size(wall_blocks, 1);
	num_blocks = size(wall_blocks, 2);

	% Now that we have all of this in line, form our jackknife blocks and cov matrix.
    wall_sum = mean(wall_blocks, 2); % Use blocks!
	%wall_unblock_err = std(wall_corr, 0, 2)/sqrt(size(wall_corr, 2));
	wall_rep = repmat(wall_sum, [1 num_blocks]);
	wall_jack = zeros(parse_Nt, num_blocks);
    for b=1:num_blocks
       file_conn_copy = wall_blocks;
       file_conn_copy(:,b) = [];
       wall_jack(:,b) = mean(file_conn_copy,2);
    end 
	clear('file_conn_copy');
    wall_cov_mat = zeros(parse_Nt);
    wall_err = zeros(parse_Nt, 1);
    for t1 = 1:parse_Nt
        for t2 = 1:parse_Nt
            wall_cov_mat(t1,t2) = sum((wall_rep(t1,:)-wall_jack(t1,:)).*(wall_rep(t2,:)-wall_jack(t2,:)),2).*(num_blocks-1)./num_blocks;
            if t1 == t2
                wall_err(t1) = sqrt(wall_cov_mat(t1,t1));
            end
        end
    end
    
	out_jack = wall_jack; out_cov_mat = wall_cov_mat; out_err = wall_err;
	
end