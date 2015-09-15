% Update 08-17-2014, Fix vev issue sort of!
% Update 08-20-2014, Read in flavor count, Ns, Nt from file.
% Update 09-18-2014, Updated for SCC analysis. Do jackknifing right!
% Update 01-17-2015, Save the dc and sg correlator to the corr directory.
% Update 09-15-2014, Use new loading interface.
function build_correlator_scalar(fname, stoch_src, blocksize)

	% Load the path with routines to load, bin, etc data.
	addpath('./process', '-end');
	

	% Load ensemble info.
	[fl_l, fl_h, parse_Ns, parse_Nt, beta, m_l, m_h] = load_info_file(fname);

	number_src = stoch_src;
    number_bl = stoch_src;
    volume = parse_Ns^3;
    fl_flavor = fl_l/4; % Because staggered.
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Load and save the vev! %
	%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	% First, get some info on the vev.
	
	[pbp_stats, vev_stats] = load_vev_scalar(fname, fl_flavor, parse_Nt, parse_Ns, number_bl, blocksize);
	
	% Save the vev info!
	
	raw_output = zeros(size(pbp_stats,1), 3);
	for j=1:size(pbp_stats,1)
		raw_output(j,1) = pbp_stats(j,1);
		raw_output(j,2) = pbp_stats(j,2);
		raw_output(j,3) = pbp_stats(j,4);
	end
	full_fname = strcat(fname, '/spectrum2/vev/history.vev1');
	save(full_fname, 'raw_output', '-ascii', '-double');
	
	raw_output = zeros(size(vev_stats,1), 3);
	for j=1:size(pbp_stats,1)
		raw_output(j,1) = vev_stats(j,1);
		raw_output(j,2) = vev_stats(j,2);
		raw_output(j,3) = vev_stats(j,4);
	end
	full_fname = strcat(fname, '/spectrum2/vev/history.vev2');
	save(full_fname, 'raw_output', '-ascii', '-double');

	% This has a new structure now!
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Load and save sum for D! %
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	% First, load the disconnected correlator.
	% This doesn't depend on the connected piece, which doesn't
	% exist yet.
	
	[disc_sum, disc_jack, disc_cov_mat, disc_err, num_blocks] = ...
		get_correlator(fname, 'dc_stoch', parse_Nt, parse_Ns, ...
		fl_flavor, blocksize);
		
		
	% Save D.
	data = zeros(parse_Nt/2+1,3);
    for i=1:(parse_Nt/2+1)
        data(i,1) = i-1;
        data(i,2) = disc_sum(i);
        data(i,3) = disc_err(i);
    end
	full_fname = strcat(fname, '/spectrum2/sum/sum.dc_stoch');
    save(full_fname,'data','-ascii', '-double');
	
	
	% Deal with some old variables kicking around. 
	parse_Nt_in = numel(disc_sum);
	%num_data = size(disc_jack, 2)*blocksize;
	%num_data_in = num_data;
	number_files_in = 1;
	number_bl_in = stoch_src;
	
	%%%%%%%%%%%%%%%%%%%%
	% Load and save C! %
	%%%%%%%%%%%%%%%%%%%%
	
	% Next, load the connected correlator.

	full_fname = strcat(fname, '/spectrum2/stoch/SPECTRUM.dat');
	fd = fopen(full_fname,'rt');
	raw_data = textscan(fd, '%d%d%f', 'Whitespace', ' \t');
	fclose(fd);
	raw_data{1} = double(raw_data{1});
	raw_data{2} = double(raw_data{2});
	raw_data{3} = double(raw_data{3});
    conn_values = cell2mat(raw_data);

	num_data = size(conn_values,1)/parse_Nt;
	num_data_in = num_data;

	% Old! Use reshape.
	%{
    conn_corr = zeros(parse_Nt, num_data);
	config_nums_conn = zeros(1, num_data);
    for i=1:num_data
        for j=1:parse_Nt
            % normalize it.
            conn_corr(j,i) = conn_values((i-1)*parse_Nt+j,3)/(parse_Nt*volume*stoch_src*(stoch_src-1)*0.5);
			config_nums_conn(1,i) = conn_values((i-1)*parse_Nt+j,1);
        end
    end
	
	connected = conn_corr;
	%}
	
	% Get corr.
	connected = reshape(conn_values(:,3), [parse_Nt, num_data])/(parse_Nt*volume*stoch_src*(stoch_src-1)*0.5);
	
	% Get cfg numbers, too.
	tmpinfo = reshape(conn_values(:,1), [parse_Nt, num_data]);
	config_nums_conn = tmpinfo(1,:);
	clear('tmpinfo');
	
	connected_unfolded = connected;

    connected = fold_data(connected, 0); % It's not a baryon!

	% Block data!
	[connected_blocks, num_blocks] = block_data(connected, 2, blocksize);
	
	% Sum, blocks, errors, etc.
	connected_sum = mean(connected_blocks, 2);
	connected_jack = jackknife_bins(connected_blocks, 2);
	[connected_cov_mat, connected_err] = errors_jackknife(connected_sum, connected_jack);
	
	% Spit the connected corr out now so we can get sg_stoch.
	
	% Save corr.sc_stoch
	conn_output = zeros(num_data_in*parse_Nt_in, 3);
    for i=1:parse_Nt_in
        for j=1:num_data_in
			conn_output((j-1)*parse_Nt_in+i, 1) = config_nums_conn(1,j);
            conn_output((j-1)*parse_Nt_in+i, 2) = i-1;
            conn_output((j-1)*parse_Nt_in+i, 3) = connected_unfolded(i,j);
		end
	end
	full_fname = strcat(fname, '/spectrum2/corr/corr.sc_stoch');
    save(full_fname, 'conn_output', '-ascii', '-double');
	
	% Save sum.sc_stoch
	data = zeros(parse_Nt/2+1,3);
    for i=1:(parse_Nt/2+1)
        data(i,1) = i-1;
        data(i,2) = connected_sum(i);
        data(i,3) = connected_err(i);
    end
	full_fname = strcat(fname, '/spectrum2/sum/sum.sc_stoch');
    save(full_fname,'data','-ascii','-double');
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Load and save sum for #D-C! %
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
    % While we're at it, also build sigma.
	
	% We put in the factor of N_f/4 earlier.
	sigma_sum = disc_sum - connected_sum;
	sigma_jack = disc_jack - connected_jack;	
	[sigma_cov_mat, sigma_err] = errors_jackknife(sigma_sum, sigma_jack);
	
	% Save!
	data = zeros(parse_Nt/2+1,3);
    for i=1:(parse_Nt/2+1)
        data(i,1) = i-1;
        data(i,2) = sigma_sum(i);
        data(i,3) = sigma_err(i);
    end

	full_fname = strcat(fname, '/spectrum2/sum/sum.sg_stoch');
    save(full_fname,'data','-ascii', '-double');
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Save #D, #D-C with global vev subtracted %
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	% New 01-17-2015: Save the dc, sg corr where the total vev is subbed from all of them.
	% This code is mostly borrowed from load_correlator_scalar.
	[realops config_nums_pbp] = load_pbppart(fname, parse_Nt, parse_Ns, number_bl);
	% Rescale the ops so we don't need to multiply D by anything later.
	realops = realops.*sqrt(fl_flavor);
	% Build all the non-vev subtracted correlators.
	% (parse_Nt, number_bl_in, number_bl_in, num_data);
	disc_corr_output = build_vev_correlator(realops, 1); % fold it.	
	% Recall realops is (parse_Nt, number_bl_in, num_data);
	% Okay, now we're ready to get the correlators! Vev subtract!
	vev_center = mean(mean(realops, 3), 1);
	vev_sub = repmat(reshape(parse_Nt.*((vev_center')*vev_center),1,number_src,number_src), [parse_Nt, 1, 1, size(disc_corr_output,4)]);
	% Subtract the vev.
	disc_corr_output = disc_corr_output - vev_sub; 
	for i=1:number_src
		disc_corr_output(:,i,i,:) = 0;
    end
  	disc_corr_reduced = sum(sum(disc_corr_output, 2),3)/(number_src*(number_src-1));
	disc_corr_fixed = zeros(parse_Nt, size(disc_corr_reduced,4));
	disc_corr_fixed(:,:) = disc_corr_reduced(:,1,1,:);
	
	sigma_corr_fixed = disc_corr_fixed - connected;
	
	% We save disc_corr_fixed, sigma_corr_fixed. 
	
	conn_output = zeros(num_data_in*parse_Nt_in, 3);
    for i=1:parse_Nt_in
        for j=1:num_data_in
			conn_output((j-1)*parse_Nt_in+i, 1) = config_nums_conn(1,j);
            conn_output((j-1)*parse_Nt_in+i, 2) = i-1;
            conn_output((j-1)*parse_Nt_in+i, 3) = disc_corr_fixed(i,j);
		end
	end
	full_fname = strcat(fname, '/spectrum2/corr/corr.dc_stoch');
    save(full_fname, 'conn_output', '-ascii', '-double');
	
	conn_output = zeros(num_data_in*parse_Nt_in, 3);
    for i=1:parse_Nt_in
        for j=1:num_data_in
			conn_output((j-1)*parse_Nt_in+i, 1) = config_nums_conn(1,j);
            conn_output((j-1)*parse_Nt_in+i, 2) = i-1;
            conn_output((j-1)*parse_Nt_in+i, 3) = sigma_corr_fixed(i,j);
		end
	end
	full_fname = strcat(fname, '/spectrum2/corr/corr.sg_stoch');
    save(full_fname, 'conn_output', '-ascii', '-double');
	
	%%%%%%%%%%%%%%
	% ppp stuff. %
	%%%%%%%%%%%%%%
	
	% NEW 09-29-2014
	% Save ppp'd correlators, too.
	disc_ppp_sum = ppp_data(disc_sum);
	disc_ppp_jack = ppp_data(disc_jack);
	[disc_ppp_cov_mat, disc_ppp_err] = errors_jackknife(disc_ppp_sum, disc_ppp_jack);
	
	
	sigma_ppp_sum = ppp_data(sigma_sum);
	sigma_ppp_jack = ppp_data(sigma_jack);
	[sigma_ppp_cov_mat, sigma_ppp_err] = errors_jackknife(sigma_ppp_sum, sigma_ppp_jack);
	
	data = zeros(parse_Nt/2+1,3);
    for i=1:(parse_Nt/2+1)
        data(i,1) = i-1;
        data(i,2) = disc_ppp_sum(i);
        data(i,3) = disc_ppp_err(i);
    end
	full_fname = strcat(fname, '/spectrum2/sum/sum.dc_stoch_ppp', '-double');
    save(full_fname,'data','-ascii');
    
	
	data = zeros(parse_Nt/2+1,3);
    for i=1:(parse_Nt/2+1)
        data(i,1) = i-1;
        data(i,2) = sigma_ppp_sum(i);
        data(i,3) = sigma_ppp_err(i);
    end
	full_fname = strcat(fname, '/spectrum2/sum/sum.sg_stoch_ppp', '-double');
    save(full_fname,'data','-ascii');
	
	%%%%%%%%%%%%%%%%%%%%%
	% Effective masses. %
	%%%%%%%%%%%%%%%%%%%%%
	
	% Get some effective masses of D, too!
	[disc_mass_sum, disc_roots_sum, disc_amps_sum] = effective_mass_utility(disc_sum, parse_Nt, 1, 2, 0);
	disc_mass_jack = zeros([size(disc_mass_sum) num_blocks]);
	for b=1:num_blocks
		[ disc_mass_jack(:,:, b), ~, ~] = effective_mass_utility(disc_jack(:,b), parse_Nt, 1, 2, 0);
	end
	[disc_mass_cov_mat, disc_mass_err] = errors_jackknife(disc_mass_sum, disc_mass_jack);
	clear('disc_mass_cov_mat');
	clear('disc_roots_sum');
	clear('disc_amps_sum');
	
	% Two state effective masses, too.
	[disc_mass2_sum, disc_roots2_sum, disc_amps2_sum] = effective_mass_utility(disc_sum, parse_Nt, 2, 4, 0);
	disc_mass2_jack = zeros([size(disc_mass2_sum) num_blocks]);
	for b=1:num_blocks
		[ disc_mass2_jack(:,:, b), ~, ~] = effective_mass_utility(disc_jack(:,b), parse_Nt, 2, 4, 0);
	end
	disc_mass2_rep = repmat(disc_mass2_sum, [1 1 num_blocks]);
	disc_mass2_err = sqrt(sum((disc_mass2_rep-disc_mass2_jack).^2, 3).*(num_blocks-1)/(num_blocks));
	clear('disc_roots2_sum');
	clear('disc_amps2_sum');
	clear('disc_mass2_rep');

	
	% NEW: 09-22-2014, get some new effective masses.
	% First, the Kuti-esque effective mass.
	

	% Kuti no zero.
	disc_kuti2_mass = zeros(parse_Nt-3,1);
	disc_kuti2_mass_jack = zeros(parse_Nt-3,size(disc_jack,2));
	for i=2:(parse_Nt-2)
		disc_kuti2_mass(i-1,1) = abs(log((disc_sum(i)+2*disc_sum(i+1)+disc_sum(i+2))./(disc_sum(i-1)+2*disc_sum(i)+disc_sum(i+1))));
		disc_kuti2_mass_jack(i-1,:) = abs(log((disc_jack(i,:)+2*disc_jack(i+1,:)+disc_jack(i+2,:))./(disc_jack(i-1,:)+2*disc_jack(i,:)+disc_jack(i+1,:))));
	end
	disc_kuti2_mass_err = sqrt((num_blocks-1)/(num_blocks)*sum((disc_kuti2_mass_jack - repmat(disc_kuti2_mass, [1 num_blocks])).^2, 2));

	% Zero the center of the correlator.
	% NEW 01-25-2015 subtract a constant center.
	disc_twid_sum = disc_sum - disc_sum(parse_Nt/2+1); 
	disc_twid_jack = disc_jack - disc_sum(parse_Nt/2+1); % changed here.
	disc_kuti_mass = zeros(parse_Nt-3,1);
	disc_kuti_mass_jack = zeros(parse_Nt-3,size(disc_jack,2));
	% Kuti's weird form.
	for i=2:(parse_Nt-2)
		disc_kuti_mass(i-1, 1) = abs(log((disc_twid_sum(i) + 2*disc_twid_sum(i+1) + disc_twid_sum(i+2))./(disc_twid_sum(i-1)+2*disc_twid_sum(i)+disc_twid_sum(i+1))));
		disc_kuti_mass_jack(i-1, :) = abs(log((disc_twid_jack(i,:) + 2*disc_twid_jack(i+1,:) + disc_twid_jack(i+2,:))./(disc_twid_jack(i-1,:)+2*disc_twid_jack(i,:)+disc_twid_jack(i+1,:))));
	end
	disc_kuti_mass_err = sqrt((num_blocks-1)/(num_blocks)*sum((disc_kuti_mass_jack - repmat(disc_kuti_mass, [1 num_blocks])).^2, 2));
	
	% And while we're here, let's get a sane KMI effective mass.
	disc_kmi_mass = zeros(parse_Nt-1,1);
	disc_kmi_mass_jack = zeros(parse_Nt-1,size(disc_jack,2));
	% Kuti's weird form.
	for i=1:(parse_Nt-1)
		disc_kmi_mass(i, 1) = abs(log(disc_sum(i)./disc_sum(i+1)));
		disc_kmi_mass_jack(i, :) = abs(log(disc_jack(i,:)./disc_jack(i+1,:)));
	end
	disc_kmi_mass_err = sqrt((num_blocks-1)/(num_blocks)*sum((disc_kmi_mass_jack - repmat(disc_kmi_mass, [1 num_blocks])).^2, 2));
	
	% NEW: 09-29-2014
	% Repeat everything on D for the sigma.
	[sigma_mass_sum, sigma_roots_sum, sigma_amps_sum] = effective_mass_utility(sigma_sum, parse_Nt, 1, 2, 0);
	sigma_mass_jack = zeros([size(sigma_mass_sum) num_blocks]);
	for b=1:num_blocks
		[ sigma_mass_jack(:,:, b), ~, ~] = effective_mass_utility(sigma_jack(:,b), parse_Nt, 1, 2, 0);
	end
	sigma_mass_rep = repmat(sigma_mass_sum, [1 1 num_blocks]);
	sigma_mass_err = sqrt(sum((sigma_mass_rep-sigma_mass_jack).^2, 3).*(num_blocks-1)/(num_blocks));
	
	% Two state effective masses, too.
	[sigma_mass2_sum, sigma_roots2_sum, sigma_amps2_sum] = effective_mass_utility(sigma_sum, parse_Nt, 2, 4, 0);
	sigma_mass2_jack = zeros([size(sigma_mass2_sum) num_blocks]);
	for b=1:num_blocks
		[ sigma_mass2_jack(:,:, b), ~, ~] = effective_mass_utility(sigma_jack(:,b), parse_Nt, 2, 4, 0);
	end
	sigma_mass2_rep = repmat(sigma_mass2_sum, [1 1 num_blocks]);
	sigma_mass2_err = sqrt(sum((sigma_mass2_rep-sigma_mass2_jack).^2, 3).*(num_blocks-1)/(num_blocks));

	% Kuti no zero.
	sigma_kuti2_mass = zeros(parse_Nt-3,1);
	sigma_kuti2_mass_jack = zeros(parse_Nt-3,size(sigma_jack,2));
	for i=2:(parse_Nt-2)
		sigma_kuti2_mass(i-1,1) = abs(log((sigma_sum(i)+2*sigma_sum(i+1)+sigma_sum(i+2))./(sigma_sum(i-1)+2*sigma_sum(i)+sigma_sum(i+1))));
		sigma_kuti2_mass_jack(i-1,:) = abs(log((sigma_jack(i,:)+2*sigma_jack(i+1,:)+sigma_jack(i+2,:))./(sigma_jack(i-1,:)+2*sigma_jack(i,:)+sigma_jack(i+1,:))));
	end
	sigma_kuti2_mass_err = sqrt((num_blocks-1)/(num_blocks)*sum((sigma_kuti2_mass_jack - repmat(sigma_kuti2_mass, [1 num_blocks])).^2, 2));

	% Zero the center of the correlator.
	sigma_twid_sum = sigma_sum - sigma_sum(parse_Nt/2+1); 
	sigma_twid_jack = sigma_jack - sigma_sum(parse_Nt/2+1);
	sigma_kuti_mass = zeros(parse_Nt-3,1);
	sigma_kuti_mass_jack = zeros(parse_Nt-3,size(sigma_jack,2));
	% Kuti's weird form.
	for i=2:(parse_Nt-2)
		sigma_kuti_mass(i-1, 1) = abs(log((sigma_twid_sum(i) + 2*sigma_twid_sum(i+1) + sigma_twid_sum(i+2))./(sigma_twid_sum(i-1)+2*sigma_twid_sum(i)+sigma_twid_sum(i+1))));
		sigma_kuti_mass_jack(i-1, :) = abs(log((sigma_twid_jack(i,:) + 2*sigma_twid_jack(i+1,:) + sigma_twid_jack(i+2,:))./(sigma_twid_jack(i-1,:)+2*sigma_twid_jack(i,:)+sigma_twid_jack(i+1,:))));
	end
	sigma_kuti_mass_err = sqrt((num_blocks-1)/(num_blocks)*sum((sigma_kuti_mass_jack - repmat(sigma_kuti_mass, [1 num_blocks])).^2, 2));
	
	% And while we're here, let's get a sane KMI effective mass.
	sigma_kmi_mass = zeros(parse_Nt-1,1);
	sigma_kmi_mass_jack = zeros(parse_Nt-1,size(sigma_jack,2));
	% Kuti's weird form.
	for i=1:(parse_Nt-1)
		sigma_kmi_mass(i, 1) = abs(log(sigma_sum(i)./sigma_sum(i+1)));
		sigma_kmi_mass_jack(i, :) = abs(log(sigma_jack(i,:)./sigma_jack(i+1,:)));
	end
	sigma_kmi_mass_err = sqrt((num_blocks-1)/(num_blocks)*sum((sigma_kmi_mass_jack - repmat(sigma_kmi_mass, [1 num_blocks])).^2, 2));
	

	% Spit out effective masses.
	
	raw_output = zeros(size(disc_mass_err,1), 3);
	for j=1:(size(disc_mass_err,1))
		raw_output(j,1) = j-1+((2-1)/2);
		raw_output(j,2) = real(disc_mass_sum(j, 1));
		raw_output(j,3) = real(disc_mass_err(j, 1));
	end
	full_fname = strcat(fname, '/spectrum2/effmass/effmass.dc_stoch', num2str(1));

	save(full_fname, 'raw_output', '-ascii', '-double');


	raw_output = zeros(size(disc_mass2_err,1), 3);
	for j=1:(size(disc_mass2_err,1))
		raw_output(j,1) = j-1+((4-1)/2);
		raw_output(j,2) = real(disc_mass2_sum(j, 1));
		raw_output(j,3) = real(disc_mass2_err(j, 1));
	end
	full_fname = strcat(fname, '/spectrum2/effmass/effmass.dc_stoch_oscil', num2str(1));
	save(full_fname, 'raw_output', '-ascii', '-double');
	
	raw_output = zeros(size(disc_mass2_err,1), 3);
	for j=1:(size(disc_mass2_err,1))
		raw_output(j,1) = j-1+((4-1)/2);
		raw_output(j,2) = real(disc_mass2_sum(j, 2));
		raw_output(j,3) = real(disc_mass2_err(j, 2));
	end
	full_fname = strcat(fname, '/spectrum2/effmass/effmass.dc_stoch_oscil', num2str(2));
	save(full_fname, 'raw_output', '-ascii', '-double');

	
    raw_output = zeros(size(disc_kuti_mass_err,1), 3);
	for j=1:(size(disc_kuti_mass_err,1))
		raw_output(j,1) = j;
		raw_output(j,2) = real(disc_kuti_mass(j, 1));
		raw_output(j,3) = real(disc_kuti_mass_err(j, 1));
	end
	full_fname = strcat(fname, '/spectrum2/effmass/effmass.dc_stoch_kuti', num2str(1));

	save(full_fname, 'raw_output', '-ascii', '-double');

	raw_output = zeros(size(disc_kuti2_mass_err,1), 3);
	for j=1:(size(disc_kuti2_mass_err,1))
		raw_output(j,1) = j;
		raw_output(j,2) = real(disc_kuti2_mass(j,1));
		raw_output(j,3) = real(disc_kuti2_mass_err(j,1));
	end

	full_fname = strcat(fname, '/spectrum2/effmass/effmass.dc_stoch_ppp', num2str(1));
	save(full_fname, 'raw_output', '-ascii', '-double');
    
	raw_output = zeros(size(disc_kmi_mass_err,1), 3);
	for j=1:(size(disc_kmi_mass_err,1))
		raw_output(j,1) = j-1;
		raw_output(j,2) = real(disc_kmi_mass(j, 1));
		raw_output(j,3) = real(disc_kmi_mass_err(j, 1));
	end
	full_fname = strcat(fname, '/spectrum2/effmass/effmass.dc_stoch_kmi', num2str(1));

	save(full_fname, 'raw_output', '-ascii', '-double');
    
	% Sigma effective masses, too.
	raw_output = zeros(size(sigma_mass_err,1), 3);
	for j=1:(size(sigma_mass_err,1))
		raw_output(j,1) = j-1+((2-1)/2);
		raw_output(j,2) = real(sigma_mass_sum(j, 1));
		raw_output(j,3) = real(sigma_mass_err(j, 1));
	end
	full_fname = strcat(fname, '/spectrum2/effmass/effmass.sg_stoch', num2str(1));

	save(full_fname, 'raw_output', '-ascii', '-double');


	raw_output = zeros(size(sigma_mass2_err,1), 3);
	for j=1:(size(sigma_mass2_err,1))
		raw_output(j,1) = j-1+((4-1)/2);
		raw_output(j,2) = real(sigma_mass2_sum(j, 1));
		raw_output(j,3) = real(sigma_mass2_err(j, 1));
	end
	full_fname = strcat(fname, '/spectrum2/effmass/effmass.sg_stoch_oscil', num2str(1));
	save(full_fname, 'raw_output', '-ascii', '-double');
	
	raw_output = zeros(size(sigma_mass2_err,1), 3);
	for j=1:(size(sigma_mass2_err,1))
		raw_output(j,1) = j-1+((4-1)/2);
		raw_output(j,2) = real(sigma_mass2_sum(j, 2));
		raw_output(j,3) = real(sigma_mass2_err(j, 2));
	end
	full_fname = strcat(fname, '/spectrum2/effmass/effmass.sg_stoch_oscil', num2str(2));
	save(full_fname, 'raw_output', '-ascii', '-double');

	
    raw_output = zeros(size(sigma_kuti_mass_err,1), 3);
	for j=1:(size(sigma_kuti_mass_err,1))
		raw_output(j,1) = j;
		raw_output(j,2) = real(sigma_kuti_mass(j, 1));
		raw_output(j,3) = real(sigma_kuti_mass_err(j, 1));
	end
	full_fname = strcat(fname, '/spectrum2/effmass/effmass.sg_stoch_kuti', num2str(1));

	save(full_fname, 'raw_output', '-ascii', '-double');

	raw_output = zeros(size(sigma_kuti2_mass_err,1), 3);
	for j=1:(size(sigma_kuti2_mass_err,1))
		raw_output(j,1) = j;
		raw_output(j,2) = real(sigma_kuti2_mass(j,1));
		raw_output(j,3) = real(sigma_kuti2_mass_err(j,1));
	end

	full_fname = strcat(fname, '/spectrum2/effmass/effmass.sg_stoch_ppp', num2str(1));
	save(full_fname, 'raw_output', '-ascii', '-double');
    
	raw_output = zeros(size(sigma_kmi_mass_err,1), 3);
	for j=1:(size(sigma_kmi_mass_err,1))
		raw_output(j,1) = j-1;
		raw_output(j,2) = real(sigma_kmi_mass(j, 1));
		raw_output(j,3) = real(sigma_kmi_mass_err(j, 1));
	end
	full_fname = strcat(fname, '/spectrum2/effmass/effmass.sg_stoch_kmi', num2str(1));

	save(full_fname, 'raw_output', '-ascii', '-double');
	
end

