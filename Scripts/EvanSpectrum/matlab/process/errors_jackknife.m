% Get the covariance and errors from jackknife values.
function [out_cov_mat, out_err] = errors_jackknife(in_sum, in_jack)
	
	parse_Nt = size(in_sum,1);
	num_blocks = size(in_jack,2);
	
	in_rep = repmat(in_sum, [1 num_blocks]);
	
	in_cov_mat = zeros(parse_Nt);
	in_err = zeros(parse_Nt,1);
	
	for t1 = 1:parse_Nt
        for t2 = 1:parse_Nt
            in_cov_mat(t1,t2) = sum((in_rep(t1,:)-in_jack(t1,:)).*(in_rep(t2,:)-in_jack(t2,:)),2).*(num_blocks-1)./num_blocks;
            if t1 == t2
                in_err(t1) = sqrt(in_cov_mat(t1,t1));
            end
        end
    end
	
	out_err = in_err;
    out_cov_mat = in_cov_mat;
	
end