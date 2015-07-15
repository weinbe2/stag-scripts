function correlator = load_correlator(fname, state, parse_Nt)
    full_fname = strcat(fname, '/spectrum2/corr/corr.', state);
	raw_data = importdata(full_fname);
    num_data = floor(size(raw_data, 1)/(parse_Nt));
    %correlator = zeros(parse_Nt, num_data);
    %for i=1:num_data
    %    for j=1:parse_Nt
    %        % normalize it.
    %        correlator(j,i) = raw_data((i-1)*parse_Nt+j,3);
    %    end
    %end
    
    correlator = reshape(raw_data(:,3), [parse_Nt, num_data]);
end