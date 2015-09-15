% Do a study to learn about the vev!
function [pbp_stats, vev_stats] = load_vev_scalar(fname, fl_flavor, parse_Nt, parse_Ns, number_bl, blocksize)
	

    [realops config_nums_pbp] = load_pbppart(fname, parse_Nt, parse_Ns, number_bl);
	num_data = size(realops, 3);
	num_blocks = floor(num_data/blocksize);
	parse_Nt_in = size(realops, 1);
    num_data_in = size(realops, 3);
    number_bl_in = size(realops, 2);
    number_files_in = size(realops, 4);
	number_src = number_bl;
	
	% NEW: 2015-01-13 Get proper autocorrelation info.
	% Reduce pbp in all the right directions.
		output_autocorr_info = zeros(2, 6);
		
		reduced_pbp = mean(mean(realops.*sqrt(fl_flavor), 2), 1);
		reduced_pbp_resize = zeros(1, num_data);
		reduced_pbp_resize(1,:) = reduced_pbp(1,1,:);
		[the_mean the_error the_error_error the_tau the_tau_err] = UWerr(transpose(reduced_pbp_resize), 1.5, size(reduced_pbp_resize,2), 'Name', 1);
		disp(strcat([num2str(the_mean) ' ' num2str(the_error) ' ' num2str(the_error_error) ' ' num2str(the_tau) ' ' num2str(the_tau_err)]))
		output_autocorr_info(1,1) = 0;
		output_autocorr_info(1,2) = the_mean;
		output_autocorr_info(1,3) = the_error;
		output_autocorr_info(1,4) = the_error_error;
		output_autocorr_info(1,5) = the_tau;
		output_autocorr_info(1,6) = the_tau_err;
		
		% Also, get proper autocorrelation of the properly scaled vacuum?
		vevfunc = @(x)(sum(sum(parse_Nt*((x')*x-diag(x.^2)),2),1)/(number_bl*(number_bl-1)));
		reduced_pbp = mean(realops.*sqrt(fl_flavor),1);
		reduced_pbp_resize = zeros(number_bl, num_data);
		reduced_pbp_resize(:,:) = reduced_pbp(1,:,:);
		[the_mean the_error the_error_error the_tau the_tau_err] = UWerr(transpose(reduced_pbp_resize), 1.5, size(reduced_pbp_resize,2), 'Name', vevfunc);
		disp(strcat([num2str(the_mean) ' ' num2str(the_error) ' ' num2str(the_error_error) ' ' num2str(the_tau) ' ' num2str(the_tau_err)]))
		output_autocorr_info(2,1) = 0;
		output_autocorr_info(2,2) = the_mean;
		output_autocorr_info(2,3) = the_error;
		output_autocorr_info(2,4) = the_error_error;
		output_autocorr_info(2,5) = the_tau;
		output_autocorr_info(2,6) = the_tau_err;
	
		full_fname = strcat(fname, '/spectrum2/uwerr/uwerr.vev');
		save(full_fname, 'output_autocorr_info', '-ascii');
	
	% Bin the operators!
	[realops_blocks num_blocks] = block_data(realops, 3, blocksize);
	
	% Rescale the ops so we don't need to multiply D by anything later.
	realops_blocks = realops_blocks.*sqrt(fl_flavor);
	
	% For simplicity.
	realops = realops_blocks;
	num_data = num_blocks;
	num_data_in = num_blocks;
	
	% Get a history plot type deal of the vev.
	tmp = mean(mean(realops, 2), 1);
	realops_reduced = zeros(1,num_data);
	realops_reduced(1,:) = tmp(1,1,:);
	
	% Get essentially what we subtract, too.
	tmp1 = mean(mean(realops(:,1:(number_bl/2),:), 2), 1);
	tmp2 = mean(mean(realops(:,(number_bl/2+1):(number_bl),:), 2), 1);
	%tmp = mean(mean(realops(:,1:(number_bl/2),:).*realops(:,(number_bl/2+1):(number_bl),:)*parse_Nt, 2), 1);
	vev_reduced = zeros(2, num_data);
	vev_reduced(1,:) = tmp1(1,1,:);
	vev_reduced(2,:) = tmp2(1,1,:);
	
        
	num_steps = floor(2*(log(num_data)/log(2)-1))*0.5;
    
    count = 1;
    pbp_stats = zeros(2*num_steps, 4); % data, mean, stddev, stderr
	
	vev_stats = zeros(2*num_steps, 4);
    
    
    for p=1:0.5:num_steps
		grab_num = ceil(2.^(p+1));
		grab_num = grab_num - mod(grab_num, 2);

		% Reduce the data.

		pbp_sub = realops_reduced(1,1:grab_num);

		% Get some basic data about pbp.

		pbp_stats(count, 1) = grab_num;
		pbp_stats(count, 2) = mean(pbp_sub);
		pbp_stats(count, 3) = std(pbp_sub);
		pbp_stats(count, 4) = std(pbp_sub)/sqrt(grab_num);

		vev_sub = vev_reduced(:,1:grab_num);
		vev_mean = mean(vev_sub,2);
		
		% Get some basic data about pbp.

		vev_stats(count, 1) = grab_num;
		vev_stats(count, 2) = vev_mean(1)*vev_mean(2)*parse_Nt;
		vev_stats(count, 3) = 2*std(pbp_sub)*parse_Nt*mean(vev_mean,1);
		vev_stats(count, 4) = 2*std(pbp_sub)*parse_Nt*mean(vev_mean,1)/sqrt(grab_num);
		
		count = count + 1;
	end
	
	pbp_stats(count, 1) = numel(realops_reduced);
	pbp_stats(count, 2) = mean(realops_reduced);
	pbp_stats(count, 3) = std(realops_reduced);
	pbp_stats(count, 4) = std(realops_reduced)/sqrt(numel(realops_reduced));
	
	vev_sub = vev_reduced;
	vev_mean = mean(vev_sub,2);
	
	% Get some basic data about pbp.

	vev_stats(count, 1) = size(vev_sub,2);
	vev_stats(count, 2) = vev_mean(1)*vev_mean(2)*parse_Nt;
	vev_stats(count, 3) = 2*std(pbp_sub)*parse_Nt*mean(vev_mean,1);
	vev_stats(count, 4) = 2*std(pbp_sub)*parse_Nt*mean(vev_mean,1)/sqrt(size(vev_sub,2));
	
end
