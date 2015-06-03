% Prepare data for a fit depending on flags.

% Throw the central value and the jackknife blocks together.
rescale_corr = cat(2, scalar_sum, scalar_jack);

% First---check if we fold!
if (fit_fold == 1) % fold!
    for tmp_i = 1:(parse_Nt/2)
        rescale_corr(tmp_i+1,:) = 0.5*(rescale_corr(tmp_i+1,:)+rescale_corr((parse_Nt-tmp_i)+1,:));
    end
    for tmp_i = (parse_Nt/2+1):(parse_Nt-1)
        rescale_corr(tmp_i+1,:) = rescale_corr(parse_Nt-tmp_i+1,:);
    end
end

% Next---check if we positive parity project it!
if (fit_ppp == 1) % positive parity project!
    tmp_meh = rescale_corr(1,:); % last term gets special treatment
    for tmp_i = 0:(parse_Nt-3)
        rescale_corr(tmp_i+1,:) = rescale_corr(tmp_i+1,:)+2*rescale_corr(tmp_i+2,:)+rescale_corr(tmp_i+3,:);
    end
    rescale_corr(parse_Nt-1,:) = rescale_corr(parse_Nt-1,:)+2*rescale_corr(parse_Nt,:)+tmp_meh;
    clear('tmp_meh');
    % Shift things forward.
    for tmp_i = (parse_Nt):-1:2
        rescale_corr(tmp_i,:) = 0.25*rescale_corr(tmp_i-1,:);
    end
    rescale_corr(1,:) = rescale_corr(2,:); % we lose the very edge.
elseif (fit_ppp == -1) % negative parity project!
    tmp_meh = rescale_corr(1,:); % last term gets special treatment
    for tmp_i = 0:(parse_Nt-3)
        rescale_corr(tmp_i+1,:) = (1-2*mod(tmp_i,2))*(-rescale_corr(tmp_i+1,:)+2*rescale_corr(tmp_i+2,:)-rescale_corr(tmp_i+3,:));
    end
    rescale_corr(parse_Nt-1,:) = rescale_corr(parse_Nt-1,:)-2*rescale_corr(parse_Nt,:)+tmp_meh;
    clear('tmp_meh');
    % Shift things forward.
    for tmp_i = (parse_Nt):-1:2
        rescale_corr(tmp_i,:) = 0.25*rescale_corr(tmp_i-1,:);
    end
    rescale_corr(1,:) = rescale_corr(2,:); % we lose the very edge.
end

% Next---check if we zero shift it!
% Zero shifting just subtracts the value at the
% center of the correlator from the mean and
% from the jackknife blocks, so the value and error
% at the center is exactly zero.

% Gap: subtract the value at the center of the mean
% from the mean and all jackknife blocks, so there's
% still errors on the center (errors are, in fact,
% unchanged. If there's a constant in
% the fit, it'll exactly absorb that.
if (fit_zero == 1) % shift to zero
	shift_middle = rescale_corr((parse_Nt/2+1),:); 
    rep_shift_middle = repmat(shift_middle, [parse_Nt 1]);
    %shift_middle = rescale_corr((parse_Nt/2+1),1); 
    %rep_shift_middle = repmat(shift_middle, [parse_Nt size(rescale_corr, 2)]);
    rescale_corr = rescale_corr - rep_shift_middle;
    clear('shift_middle'); clear('rep_shift_middle');
elseif (fit_zero == -1) % 'gap' it.
	shift_middle = rescale_corr((parse_Nt/2+1),1); 
    rep_shift_middle = repmat(shift_middle, [parse_Nt size(rescale_corr, 2)]);
    rescale_corr = rescale_corr - rep_shift_middle;
    clear('shift_middle'); clear('rep_shift_middle');
end

% Once all these things are done, get the error.

rescale_sum = rescale_corr(:,1);
rescale_jack = rescale_corr(:,2:end);
[rescale_cov_mat, rescale_err] = ...
	errors_jackknife(rescale_sum, rescale_jack);






