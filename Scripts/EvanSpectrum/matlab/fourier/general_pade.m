function yval = general_pade(coeff, xval, up_degree, down_degree)
    
    % Get top half.
    top = polyval(cat(1, coeff(up_degree:-1:1), [1]), xval);
    
    % Get bottom half.
    bottom = polyval(coeff((up_degree+down_degree+1):-1:(up_degree+1)), xval);
    
    yval = top./bottom;

end