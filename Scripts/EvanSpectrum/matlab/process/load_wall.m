function [pbp_w, conn_ww, pion_ww] = load_wall(fname, parse_Nt, parse_Ns)
    
	volume = (parse_Ns.^3)*parse_Nt;
	
    % PBPPART_WALL.dat: disconnected piece of wall-wall.
    % Business as usual.
	full_fname = strcat(fname, '/spectrum2/wall/PBPPART_WALL.dat');
	fd = fopen(full_fname,'rt');
	raw_data = textscan(fd, '%d%d%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f', 'Whitespace', ' \t');
	fclose(fd);
    for i=1:18
        raw_data{i} = double(raw_data{i});
    end
	raw_data = cell2mat(raw_data);
    %raw_data = importdata(full_fname);
	num_data = floor(size(raw_data, 1)/(parse_Nt));
	
	raw_data = sum(raw_data(:,3:2:17),2); % sum over all 8 hypercube ops.
	pbp_w = reshape(raw_data(:,1), [parse_Nt 1 num_data]);
    
    % SPECTRUM_WALL.dat: connected piece of wall-wall.
	full_fname = strcat(fname, '/spectrum2/wall/SPECTRUM_WALL.dat');
	fd = fopen(full_fname,'rt');
	raw_data = textscan(fd, '%d%d%f%f', 'Whitespace', ' \t');
	fclose(fd);
	raw_data{1} = double(raw_data{1});
	raw_data{2} = double(raw_data{2});
   	raw_data{3} = double(raw_data{3});
	raw_data{4} = double(raw_data{4});
    conn_values = cell2mat(raw_data);

	num_data = size(conn_values,1)/parse_Nt;

	conn_ww = -reshape(conn_values(:, 3), [parse_Nt num_data]);
	
	% Due to a bug in the output.
	conn_ww = circshift(conn_ww, [1 0]);
	
	% SPECTRUM_WALL_PION.dat: connected piece of wall-wall.
	full_fname = strcat(fname, '/spectrum2/wall/SPECTRUM_WALL_PION.dat');
	fd = fopen(full_fname,'rt');
	raw_data = textscan(fd, '%d%d%f%f', 'Whitespace', ' \t');
	fclose(fd);
	raw_data{1} = double(raw_data{1});
	raw_data{2} = double(raw_data{2});
   	raw_data{3} = double(raw_data{3});
	raw_data{4} = double(raw_data{4});
    conn_values = cell2mat(raw_data);

	num_data = size(conn_values,1)/parse_Nt;

	pion_ww = reshape(conn_values(:,3), [parse_Nt num_data]);
	
    pion_ww = circshift(pion_ww, [1 0]);
	
    num_data_in = num_data; 
    
	% Do some rescaling by 3-volume.
	pbp_w = pbp_w/(volume/parse_Nt);
	conn_ww = conn_ww/((volume/parse_Nt)^2);
	pion_ww = pion_ww/((volume/parse_Nt)^2);
	
	% Send 'em back!
end