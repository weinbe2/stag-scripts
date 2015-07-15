% Return <pbp>. If requested, also return the configuration.
function [correlator varargout] = load_pbppart(fname, parse_Nt, parse_Ns, number_bl)
	% Number of extra arguments.
	nout = max(nargout,1) - 1;

	% Business as usual.
	full_fname = strcat(fname, '/spectrum2/stoch/PBPPART.dat');
	fd = fopen(full_fname,'rt');
	raw_data = textscan(fd, '%d%d%d%f%f', 'Whitespace', ' \t');
	fclose(fd);    
    
	raw_data{1} = double(raw_data{1});
	raw_data{2} = double(raw_data{2});
	raw_data{3} = double(raw_data{3});
	raw_data = cell2mat(raw_data);
    %raw_data = importdata(full_fname);
	num_data = floor(size(raw_data, 1)/(parse_Nt*number_bl));
    %correlator = zeros(parse_Nt, number_bl, num_data);
	
	if (nout == 1)
		varargout{1} = zeros(1, num_data);
    end
	
    % Old! The reshape method is faster.
    %{
	for i=1:num_data
        for j=1:number_bl
            for k=1:parse_Nt
                correlator(k,j,i) = raw_data((i-1)*number_bl*parse_Nt+(j-1)*parse_Nt+k, 5)*sqrt((parse_Ns)^3*parse_Nt);
				if (nout == 1)
					varargout{1}(1,i) = raw_data((i-1)*number_bl*parse_Nt+(j-1)*parse_Nt+k, 1);
				end
            end
        end    
    end
    %}
	

	correlator = reshape(raw_data(:,5), [parse_Nt, number_bl, num_data])*sqrt((parse_Ns)^3*parse_Nt);
    if (nout == 1)
        tmpinfo = reshape(raw_data(:,1), [parse_Nt, number_bl, num_data]);
        varargout{1}(1,:) = tmpinfo(1,1,:);
        clear('tmpinfo');
    end
        
    
end

