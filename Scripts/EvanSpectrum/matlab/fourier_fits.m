function outform = fourier_fits_new(fname, state, num_roots, noerrors, blockval, tmin_min, tmin_max, tmax)
	% Load the path with routines to load, bin, etc data.
	addpath('.\process', '-end');
	% Get routines for fourier fits.
	addpath('.\fourier', '-end');

	% Load ensemble info.
	[fl_l, fl_h, parse_Ns, parse_Nt, beta, m_l, m_h] = load_info_file(fname);
	
	min_search = 0; max_search = 0; tmax_search = 0; blocksize = 0;
	check_errs = 0; diag_only = 0; number_bl = 6;
	fl_flavor = fl_l/4;
	
	if (exist('noerrors', 'var'))
		if (mod(noerrors, 2) == 1)
			check_errs = 1;
			noerrors = noerrors - 1;
		end
		noerrors = noerrors / 2;
		if (noerrors == 1)
			diag_only = 1;
        end
    end
	
	% If we didn't get range arguments...
	if (~exist('blockval', 'var') && ~exist('tmin_min', 'var') && ~exist('tmin_max', 'var') && ~exist('tmax', 'var'))
		% See if a fitparams file exists.
		if (exist(fullfile(cd, strcat([fname '/spectrum2/fitparams/fitparam.' state])), 'file'))
			fitparam = importdata(strcat([fname '/spectrum2/fitparams/fitparam.' state]));
			min_search = fitparam(1);
			max_search = fitparam(3);
			tmax_search = fitparam(4);
			blocksize = fitparam(5);
			
		else
			error('A fitparams file is expected but does not exist!');
		end
	else
		min_search = tmin_min;
		max_search = tmin_max;
		tmax_search = tmax;
		blocksize = blockval;
	end
		
	mass_l = m_l;
  
	% Figure out what fit type we need.
	fit_type = 3; % default cosh+oscil
  
	% Learn about our measurement from the state name.
	is_baryon = 0;
	if (strcmp(state, 'nu') || strcmp(state, 'de'))
		is_baryon = 1;
		fit_type = 12; % baryon
	end
	
	is_fpi = 0;
	if (strcmp(state, 'ps2'))
		is_fpi = 1;
    end
	
	is_oscil = 1;
	if (strcmp(state, 'ps') || strcmp(state, 'ps2') || strcmp(state, 'i5'))
		is_oscil = 0;
		fit_type = 1; % single cosh
    end
	
	% Load everything.
	[connected_sum, connected_jack, connected_cov_mat, connected_err, num_blocks] = get_correlator(fname, state, parse_Nt, parse_Ns, fl_flavor, blocksize);
	clear('connected_cov_mat');
	clear('connected_err');
	
	% Fix sign of sc_stoch.
	if (strcmp(state, 'sc_stoch'))
		connected_sum = -1*connected_sum;
		connected_jack = -1*connected_jack;
	end
	
	% Fourier transform, take the real part to fold it, do the analysis.
	% /sqrt(parse_Nt) is Claudio's convention.
	connected_P_sum = real(fft(connected_sum, [], 1))/sqrt(parse_Nt);
	connected_P_jack = real(fft(connected_jack, [], 1))/sqrt(parse_Nt);
	[connected_P_cov_mat, connected_P_err] = errors_jackknife(connected_P_sum, connected_P_jack);
	
	%{
	
    % Extend to also handle disconnected.
    
    if (strcmp(state, 'dc_stoch') || strcmp(state, 'sg_stoch'))
        % Get that disconnected state!
        
		% Load the pbp values. Recall pbp is (parse_Nt, number_bl_in, num_data);
		[pbp config_nums_pbp] = load_pbppart(fname, parse_Nt, parse_Ns, number_bl);
		num_data = size(pbp, 3);
		
        % Rescale the ops so we don't need to multiply D by anything later.
		pbp = pbp.*sqrt(fl_flavor);
		
		% Build all the non-vev subtracted correlators.
		% (parse_Nt, number_bl_in, number_bl_in, num_data);
		disc = build_vev_correlator(pbp, 1); % fold it.
		
		% Now bin it!
		[pbp_blocks num_blocks] = block_data(pbp, 3, blocksize);
		[disc_blocks num_blocks] = block_data(disc, 4, blocksize);
		
		% Get central values.
		pbp_sum = mean(pbp_blocks, 3);
		disc_sum = mean(disc_blocks, 4);
		
		% Single elim jackknife.
		pbp_jack = jackknife_bins(pbp_blocks, 3);
		disc_jack = jackknife_bins(disc_blocks, 4);
		
		% Build the disconnected correlators.
		vev_center = mean(pbp_sum, 1);
		vev_sub = repmat(reshape(parse_Nt.*((vev_center')*vev_center),1,number_bl,number_bl), [parse_Nt, 1, 1]);
		disc_sum = disc_sum - vev_sub;
		for i=1:number_bl
			disc_sum(:,i,i,:) = 0;
		end
		disc_sum = sum(sum(disc_sum,2),3)/(number_bl*(number_bl-1));
		
		% Next on the jackknife!
		vev_center = mean(pbp_jack, 1);
		vev_1 = repmat(reshape(vev_center,1, number_bl, 1, num_blocks), [parse_Nt, 1, number_bl, 1]);
		vev_2 = repmat(reshape(vev_center,1, 1, number_bl, num_blocks), [parse_Nt, number_bl, 1, 1]);
		disc_jack = disc_jack - vev_1.*vev_2.*parse_Nt;
		for i=1:number_bl
			disc_jack(:,i,i,:) = 0;
		end
		disc_jack2 = sum(sum(disc_jack,2),3)/(number_bl*(number_bl-1));
		
		disc_jack = zeros(parse_Nt, num_blocks);
		disc_jack(:,:) = disc_jack2(:,1,1,:);
		
		% Grab the connected piece if we need to!
		if (strcmp(state, 'sg_stoch'))
			
			% Load it, fold it...
			connected = load_correlator(fname, 'sc_stoch', parse_Nt);
			connected = fold_data(connected, is_baryon);
			
			% Bin it!
			connected_blocks = block_data(connected, 2, blocksize);
			
			% Sum it, jack it.
			connected_sum = mean(connected_blocks, 2);
			connected_jack = jackknife_bins(connected_blocks, 2);
			
			% And modify disc.
			disc_sum = disc_sum - connected_sum;
			disc_jack = disc_jack - connected_jack;
			
		end
		
		% Now give the rest of the code what it wants!
		connected_P_sum = real(fft(disc_sum, [], 1))/sqrt(parse_Nt); % match Claudio's conventions.
		connected_P_jack = real(fft(disc_jack, [], 1))/sqrt(parse_Nt);
		[connected_P_cov_mat, connected_P_err] = errors_jackknife(connected_P_sum, connected_P_jack);
		
    else

        % Load the correlator, run autocorrelation, block it.
        connected = load_correlator(fname, state, parse_Nt);

        % Fold the correlator, appropriately handling if it's a baryon.

        if (strcmp(state, 'sc_stoch'))
            connected = -1*fold_data(connected, is_baryon);
        else
            connected = fold_data(connected, is_baryon);
        end

        % Take the fourier transform!

        connected_P = real(fft(connected, [], 1));
        connected_P = connected_P/sqrt(parse_Nt); % match claudio's convention

        % Block it.
        connected_P_blocks = block_data(connected_P, 2, blocksize);
        num_blocks = size(connected_P_blocks, 2);

        % Sum, blocks, errors, etc.
        connected_P_sum = mean(connected_P_blocks, 2);
        connected_P_jack = jackknife_bins(connected_P_blocks, 2);
        [connected_P_cov_mat, connected_P_err] = errors_jackknife(connected_P_sum, connected_P_jack);
    end
	
	%}
    
    % Plot it.
    hold on; axis([-pi/parse_Nt pi 0 max(connected_P_sum)*1.25]); errorbar((0:(parse_Nt-1))/parse_Nt*2*pi, connected_P_sum, connected_P_err, '.k'); hold off;

    %pause
    
    % Get some fits going?
	
    xval = ((0:(parse_Nt/2))/parse_Nt*2*pi);
	yval = connected_P_sum(1:(parse_Nt/2+1))';
	ycorr = connected_P_cov_mat(1:(parse_Nt/2+1),1:(parse_Nt/2+1));
    
    %{
    
    % We can ramp this up to an arbitrary order intelligently.
    
    %ycorr = diag(diag(ycorr));
    
    % Get initial guesses for a one pole fit. Picked various locations
    % heuristically.
    guess = inv([-cos(xval(2)), yval(2), cos(xval(2))*yval(2); ...
            -cos(xval(parse_Nt/8+1)), yval(parse_Nt/8+1), cos(xval(parse_Nt/8+1))*yval(parse_Nt/8+1); ...
            -cos(xval(3*parse_Nt/8+1)), yval(3*parse_Nt/8+1), cos(xval(3*parse_Nt/8+1))*yval(3*parse_Nt/8+1);]) * ...
            [1;1;1];
        
    % Do the 1 pole fit. 1,1 means 1st order in top, 1st order in bottom.
    padefunc = @(x,xval)(general_pade(x,cos(xval),1,1));
    ycorrtmp = diag(diag(ycorr));
    chisqfunc = @(x)(1.0/(size(yval,2)-2*1-1).*(yval-padefunc(x,xval))*(ycorrtmp\((yval-padefunc(x,xval))')));
    [temp_guess, chisq, succ_code] = fminsearch(chisqfunc, guess, optimset('TolX', 1e-15,'TolFun', 1e-7,'MaxFunEvals',100000,'MaxIter',100000));
    
    errorbar(xval,yval,sqrt(diag(ycorr)),'.k'); hold on;
    title('1 pole fit.');
    axis([-pi/parse_Nt pi 0 1]);
    tmp=general_pade(temp_guess,cos(xval),1,1);  plot (xval, tmp); hold off;
    
    pause;
    
    % Now we ramp this up, adding parameters, up to num_roots.
    if (num_roots>1)
        for i=2:num_roots
            
            % Add some extra zeros into the guess.
            guess = [temp_guess(1:(i-1)); 0; temp_guess(i:(2*i-1)); 0];
            
            % Do the i pole fit. i,i means i'th order in top, i'th order in
            % bottom. We just do these as uncorrelated fits for stability.
            padefunc = @(x,xval)(general_pade(x,cos(xval),i,i));
            ycorrtmp = diag(diag(ycorr));
            chisqfunc = @(x)(1.0/(size(yval,2)-2*i-1).*(yval-padefunc(x,xval))*(ycorrtmp\((yval-padefunc(x,xval))')));
            [temp_guess, chisq, succ_code] = fminsearch(chisqfunc, guess, optimset('TolX', 1e-15,'TolFun', 1e-7,'MaxFunEvals',100000,'MaxIter',100000));

            % Plot it!
            errorbar(xval,yval,sqrt(diag(ycorr)),'.k');
            hold on;
            title(strcat(num2str(i), ' pole fit.'));
            axis([-pi/parse_Nt pi 0 1]);
            tmp=general_pade(temp_guess,cos(xval),i,i);  plot (xval, tmp); hold off;

            pause;
            
        end % for
    end % if
    
    %}
    
    % Use 'ratpolyfit' to get initial guesses.
	[num, denom] = ratpolyfit(cos(xval), yval, num_roots-1, num_roots);
    %[num, denom] = ratpolyfit(cos(xval(1:(parse_Nt/4-1))), yval(1:(parse_Nt/4-1)), num_roots-1, num_roots);
    
    % Convert to my convention.
    denom = denom./num(size(num,2));
    num = num./num(size(num,2));
    num(size(num,2)) = [];
    
    temp_guess = [num(size(num,2):-1:1), denom(size(denom,2):-1:1)];
        
    % Now redo the fit as a correlated fit.
    if (diag_only == 1)
        ycorr = diag(diag(ycorr));
    end
    
    padefunc = @(x,xval)(general_pade(x,cos(xval),num_roots-1,num_roots));
    chisqfunc = @(x)(1.0/(size(yval,2)-2*num_roots).*(yval-padefunc(x,xval))*(ycorr\((yval-padefunc(x,xval))')));
    [temp_guess, chisq, succ_code] = fminsearch(chisqfunc, temp_guess', optimset('TolX', 1e-15,'TolFun', 1e-7,'MaxFunEvals',100000,'MaxIter',100000));
    [temp_guess, chisq, succ_code] = fminsearch(chisqfunc, temp_guess, optimset('TolX', 1e-15,'TolFun', 1e-7,'MaxFunEvals',100000,'MaxIter',100000));
    [temp_guess, chisq, succ_code] = fminsearch(chisqfunc, temp_guess, optimset('TolX', 1e-15,'TolFun', 1e-7,'MaxFunEvals',100000,'MaxIter',100000));
    
    errorbar(xval,yval,sqrt(diag(ycorr)),'.k');
    hold on;
    title(strcat(num2str(num_roots), ' pole correlated fit.'));
    axis([-pi/parse_Nt pi 0 max(connected_P_sum)*1.25]);
    tmp=general_pade(temp_guess,cos(xval),num_roots-1,num_roots);  plot (xval, tmp); hold off;

    pause;
    
    results_vals = temp_guess; 
    
    % Get roots of denominator.
    
    temp_guess(1:(num_roots-1))
    temp_guess((num_roots):(2*num_roots))
    
    rvals = roots(temp_guess((2*num_roots):-1:(num_roots)));
    
    masses_sum = acosh(rvals);
    
    masses_sum
    
    chisq
    
    pval = 1.0-chi2cdf(chisq*(size(yval,2)-2*num_roots), size(yval,2)-2*num_roots)
    
    %return
    
    if (check_errs == 1)
        % I guess do it all under a jackknife?
        masses_jack = zeros(num_roots, num_blocks);

        %tic
        
        for b=1:num_blocks   
            yval = connected_P_jack(1:(parse_Nt/2+1),b)';

            padefunc = @(x,xval)(general_pade(x,cos(xval),num_roots-1,num_roots));
            chisqfunc = @(x)(1.0/(size(yval,2)-2*num_roots).*(yval-padefunc(x,xval))*(ycorr\((yval-padefunc(x,xval))')));
            [temp_guess, chisq, succ_code] = fminsearch(chisqfunc, results_vals, optimset('TolX', 1e-15,'TolFun', 1e-7,'MaxFunEvals',100000,'MaxIter',100000));
            [temp_guess, chisq, succ_code] = fminsearch(chisqfunc, temp_guess, optimset('TolX', 1e-15,'TolFun', 1e-7,'MaxFunEvals',100000,'MaxIter',100000));
            [temp_guess, chisq, succ_code] = fminsearch(chisqfunc, temp_guess, optimset('TolX', 1e-15,'TolFun', 1e-7,'MaxFunEvals',100000,'MaxIter',100000));

            rvals = roots(temp_guess((2*num_roots):-1:(num_roots)));
            masses_jack(:,b) = acosh(rvals);
            
            %{
            if (mod(b,50) == 0)
                toc
                tic
            end
            %}
            
        end
        %toc

        [masses_cov_mat, masses_err] = errors_jackknife(masses_sum, masses_jack);
		
		masses_err

    end
    
    
    
    meh=1;
    
    % claudio's values
    %temp_guess = [-0.869425, 0.0294295, 4.11378, -5.70907, 1.82169];
    
    % Get the roots a different way...
    %{
    b0 = temp_guess(3); %const
    b1 = temp_guess(4); %lin
    b2 = temp_guess(5); %quad
    
    r1 = (-b1 - sqrt(b1*b1-4*b0*b2))/(2*b2);
    r2 = (-b1 + sqrt(b1*b1-4*b0*b2))/(2*b2);
    
    m1 = log(r1+sqrt(r1*r1-1));
    m2 = log(r2+sqrt(r2*r2-1));
    %}
    
    
    
 end
    