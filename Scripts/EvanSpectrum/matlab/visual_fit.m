function visual_fit(blockval)
	% A windowed-ish fitting program!

	% Load some useful directories!
	addpath('.\process', '-end');
	addpath('.\multifit', '-end');
	
	% First, we pick a directory.


	directory_answer = inputdlg('Enter an input directory.', 'Input', 1, {'Eight\f8l24t48b48m00889'});
	if (size(directory_answer,1) == 0)
		return;
	end
		
	while (exist(strcat('..\..\..\', directory_answer{1}),'dir') == 0)

		directory_answer = inputdlg('Directory does not exist! Enter an input diretory.', 'Error!', 1, directory_answer);
		
		if (size(directory_answer,1) == 0)
			return;
		end
		
	end
	
	% Load ensemble info.
	%[fl_l, fl_h, parse_Ns, parse_Nt, beta, m_l, m_h] = load_info_file(fname);
    [fl_l, fl_h, parse_Ns, parse_Nt, beta, m_l, m_h] = load_info_file(strcat(['..\..\..\' directory_answer{1}]));

	% Get all files in said directory.
	files_in_dir = dir(strcat('..\..\..\', directory_answer{1}, '\spectrum2\corr\'));
	spectrum_list = {};
	for i=1:size(files_in_dir,1)
		[tmppath, tmpname, tmpext] = fileparts(files_in_dir(i).name);
		if (strcmp(tmpname, 'corr')) % if the file has the right name
			spectrum_list{size(spectrum_list,2)+1} = tmpext(2:end);
		end
    end
    
    % If a spectrum.dat and a pbppart.dat file exist, add "dc_stoch" and
    % "sg_stoch" to the list of states.
    if (exist(strcat(['..\..\..\' directory_answer{1} '\spectrum2\stoch\PBPPART.dat']), 'file') && exist(strcat(['..\..\..\' directory_answer{1} '\spectrum2\stoch\SPECTRUM.dat']), 'file'))
        % We can build dc_stoch and sg_stoch
        spectrum_list{size(spectrum_list,2)+1} = 'dc_stoch';
        spectrum_list{size(spectrum_list,2)+1} = 'sg_stoch';
    end

	clear('files_in_dir'); clear('tmppath'); clear('tmpname'); clear('tmpext');




	% Once the directory exists, find which input code to look at.
	%spectrum_list = {'ps', 'sc', 'i5', 'ij', 'r0', 'ris', 'rij', 'ri5', 'ps2', 'nu', 'de', 'sg', 'dc', 'sgv', 'dcv'};
	spectrum_answer = listdlg('PromptString','Enter a spectrum to look at.','SelectionMode', 'single', 'ListString', spectrum_list);

	if (isempty(spectrum_answer))
		return;
	end

	spectrum_text = spectrum_list{spectrum_answer};

	% Start loading the data!
	blocksize = 0;
	
	% If we didn't get range arguments...
	if (~exist('blockval', 'var'))
		% See if a fitparams file exists.
		if (exist(fullfile(cd, strcat(['..\..\..\' directory_answer{1} '\spectrum2\fitparams\fitparam.' spectrum_text])), 'file'))
			fitparam = importdata(strcat(['..\..\..\' directory_answer{1} '\spectrum2\fitparams\fitparam.' spectrum_text]));
			blocksize = fitparam(5);
			
		else
			error('A fitparams file is expected but does not exist!');
		end
	else
		blocksize = blockval;
	end
		
	mass_l = m_l;

	number_files = 1;
	%num_data = 600;
	%number_bl = 4;
	%parse_Nt = 48;
	%blocksize = 10;
	volume = parse_Ns^3; %*48;
	fl_flavor = fl_l/4;
	
	% Import data.
    
    [scalar_sum, scalar_jack, scalar_cov_mat, ...
		scalar_err, num_blocks] = ...
		get_correlator(strcat(['../../../' directory_answer{1}]), ...
		spectrum_text, parse_Nt, parse_Ns, fl_flavor, blocksize);
	
	% Get rolling!

	xrange = (0:(parse_Nt-1))';

    %amp_1, mass_1, amp_2, mass_2, amp_3, mass_3, amp_4, mass_4;
    coefficients = zeros(12, 1);
	
	% New: add constraints.
	constraints = zeros(12, 1);
	
	% Errors and things!
	errors = zeros(12, 1);
	
	% Hold onto saved results.
	saving = 0; 
	results_save = zeros(1, 40);
    
    % Are we blocking?
    is_blocking = 0;
	
	%amp_1 = 0.0;
	%mass_1 = 0.0;
	%amp_2 = 0.0;
	%mass_2 = 0.0;
	chisq_dof = 0.0;
	p_val = 0.0;
	cond_num = 0.0;

	form_list = {'One cosh', ...
					'Two cosh', ...
					'Three cosh', ...
					'One Oscil', ...
					'One cosh + One Oscil', ...
					'Two cosh + One Oscil', ...
					'Three cosh + One Oscil', ...
					'Two Oscil', ...
					'One cosh + Two Oscil', ...
					'Two cosh + Two Oscil', ...
					'Three cosh + Two Oscil', ...
					'Three Oscil', ...
					'One cosh + Three Oscil', ...
					'Two cosh + Three Oscil', ...
					'Three cosh + Three Oscil'};
					
	rescale_y = scalar_sum;
	rescale_yerr = scalar_err;
	mass_answer = 0; % Effectively rescaling by a zero mass, which gives 1.
	vev_answer = 0; % subtract off this vev
	fit_form = 1; % One cosh. (2 = Two cosh, 3 = One cosh, one osc, 4 = Baryon, 5 = Cosh + const, 6 = cosh with zero shift).

	% Options: tmin, tmax
	% Fold? y/n
	% Positive Parity Project? y/n
	% Fit to only even time data? y/n
	% Full or diagonal correlator matrix? 

	fit_minimum = 5;
	fit_maximum = parse_Nt-fit_minimum;
	fit_fold = 0; % don't fold
	fit_ppp = 0; % don't positive parity project
	fit_zero = 0; % don't normalize center to zero.
	fit_even = 0; % fit all
	fit_diag = 0; % fit full correlation matrix instead of diagonal.
	fit_x_prec = 1e-10; % set the x precision of the fit.

	h = figure();
	movegui(h, 'northwest');
	set(h, 'Name', 'Correlator');
	errorbar(xrange, scalar_sum, scalar_err, '.k');  hold on;
	[funcdata] = get_function(xrange, parse_Nt, fit_form, coefficients);
	plot(xrange, funcdata); hold off;
	axis([-1 parse_Nt -max(abs(scalar_sum))*1.1 max(abs(scalar_sum))*1.1]);

	h2 = figure('menu', 'none', 'toolbar', 'none', 'units', 'normalized', 'position', [0 0 .3 .8]);
	movegui(h2, 'northeast');
	set (h2, 'Name', 'Info');
	ph = uipanel(h2, 'Units', 'normalized', 'position', [0.05 0.05 0.9 0.9], 'title', 'Display Window');
	th = uicontrol(ph, 'style','text','Units','normalized','position',[0 0 1 1],'FontSize', 9, 'string', {'Rescale mass: 0.0', 'Fit Form: One cosh', 'Amplitude 1: 0.0', 'Mass 1: 0.0', 'Amplitude 2: 0.0', 'Mass 2: 0.0', 'Chisq per DoF: 0.0'}, 'horizontalalignment', 'left', 'FontName', 'FixedWidth');


	run render_update;

	figure(h);


	% Enter the main loop of options.

	flag = 1;

	while (flag == 1)

		option_list = {'Rescale+Shift (Visual Only)', ...
						'Set Fit Form', ...
						'Set Fit Options', ...
						'Set Initial Guess', ...
						'Set Constraints', ...
						'Perform Fit', ...
						'Perform Jackknife', ...
						'Begin Saving Results', ...
						'End Saving Results', ...
						'Change Directory', ...
						'Change State', ...
						'Reset Fit Options', ...
						'Reset Initial Guess', ...
						'Reset Constraints', ...
                        'Exit'};
		option_answer = listdlg('PromptString','Choose an option.','SelectionMode', 'single', 'ListString', option_list);
		
		if (isempty(option_answer))
			option_answer = 11;
		end
		
		switch option_answer
			case 1 % Rescale (Purely Visual Effect)
				
				new_mass_answer = inputdlg({'Enter a mass to rescale by.','Enter a vev to subtract.'}, 'Input', 1, {num2str(mass_answer),num2str(vev_answer)});
				
				if (~(size(new_mass_answer, 1) == 0))
					tmpmass = str2num(new_mass_answer{1});
					if (~(size(tmpmass, 1) == 0))
						% Rescale!
						mass_answer = tmpmass;
						
						
					end
					
					tmpvev = str2num(new_mass_answer{2});
					if (~(size(tmpvev, 1) == 0))
						% Rescale!
						vev_answer = tmpvev;
						
					end
					
					
				end
				
				run render_update;
				
				clear('new_mass_answer');
				clear('tmpmass');
				clear('rescale_vals');
			case 2 % Set Fit Form
				form_answer = listdlg('PromptString','Choose a fit form.','SelectionMode', 'single', 'ListString', form_list);
				
				if (~isempty(form_answer))
					fit_form = form_answer;
					
					run render_update;
				end
				
				clear('form_answer'); 

			case 3 % Set Fit Options
				
				% Options: tmin, tmax
				% Fold? y/n
				% Positive Parity Project? y/n
				% Fit to only even time data? y/n
				% Full or diagonal correlator matrix? 
				
			
				new_fit_min = inputdlg({'Minimum Fit t:', 'Maximum Fit t (-1 for symmetric fit):', 'Fit X Precion:', ...
										'Fold (0 no, 1 yes):', 'Parity Project (0 no, 1 positive, -1 negative):', 'Zero Center (0 no, 1 yes)', ...
										'Fit Every Other (0 no, 1 yes)', 'Diagonal Correlator (0 no, 1 yes)'}, ...
										'Input', 1, {num2str(fit_minimum), num2str(fit_maximum), num2str(fit_x_prec), num2str(fit_fold), num2str(fit_ppp), num2str(fit_zero), num2str(fit_even), num2str(fit_diag)});
				
				%fit_minimum = 5;
				%fit_maximum = parse_Nt-fit_minimum;
				%fit_fold = 0; % don't fold
				%fit_ppp = 0; % don't positive parity project
				%fit_zero = 0;
				%fit_even = 0; % fit all
				%fit_diag = 0; % fit full correlation matrix instead of diagonal.
				
				
				if (~(size(new_fit_min, 1) == 0)) % Make sure we didn't get a cancel!
					tmpfitmin = str2num(new_fit_min{1});
					if (~(size(tmpfitmin, 1) == 0))
						if (tmpfitmin >= 0 && tmpfitmin < (parse_Nt/2))
							fit_minimum = tmpfitmin;
						end
					end
					
					tmpfitmax = str2num(new_fit_min{2});
					if (~(size(tmpfitmax, 1) == 0))
						if (tmpfitmax > fit_minimum && tmpfitmax < (parse_Nt))
							fit_maximum = tmpfitmax;
						elseif (tmpfitmax == -1)
							fit_maximum = parse_Nt - fit_minimum;
						end
					end
					
					
					tmpxprec = str2num(new_fit_min{3});
					if (~(size(tmpxprec, 1) == 0))
						if (tmpxprec > 0)
							fit_x_prec = tmpxprec;
						end
					end
					
					tmpfitfold = str2num(new_fit_min{4});
					if (~(size(tmpfitfold, 1) == 0))
						if (tmpfitfold == 0 || tmpfitfold == 1)
							fit_fold = tmpfitfold;
						end
					end
					
					tmpfitppp = str2num(new_fit_min{5});
					if (~(size(tmpfitppp, 1) == 0))
						if (tmpfitppp == 0 || tmpfitppp == 1 || tmpfitppp == -1)
							fit_ppp = tmpfitppp;
						end
					end
					
					tmpfitzero = str2num(new_fit_min{6});
					if (~(size(tmpfitzero, 1) == 0))
						if (tmpfitzero == 0 || tmpfitzero == 1)
							fit_zero = tmpfitzero;
						end
					end
					
					tmpfiteven = str2num(new_fit_min{7});
					if (~(size(tmpfiteven, 1) == 0))
						if (tmpfiteven == 0 || tmpfiteven == 1)
							fit_even = tmpfiteven;
						end
					end
					
					tmpfitdiag = str2num(new_fit_min{8});
					if (~(size(tmpfitdiag, 1) == 0))
						if (tmpfitdiag == 0 || tmpfitdiag == 1)
							fit_diag = tmpfitdiag;
						end
					end
					
					
				end
				
				clear('tmpfitmin');
				clear('tmpfitmax');
				clear('tmpfitfold');
				clear('tmpfitppp');
				clear('tmpfitzero');
				clear('tmpfiteven');
				clear('tmpfitdiag');
				clear('tmpxprec');
				
				run render_update;
				
			case 12 % Reset Fit Options
			
				fit_minimum = 5;
				fit_maximum = parse_Nt-fit_minimum;
				fit_fold = 0; % don't fold
				fit_ppp = 0; % don't positive parity project
				fit_zero = 0; % don't zero the center of the correlator.
				fit_even = 0; % fit all
				fit_diag = 0; % fit full correlation matrix instead of diagonal.

				run render_update;
				
			case 4 % Set Initial Params
				new_guess_answers = inputdlg({'Cosh: Amplitude 1:', 'Cosh: Mass 1:', 'Cosh: Amplitude 2:', 'Cosh: Mass 2:', 'Cosh: Amplitude 3:', 'Cosh: Mass 3:',...
                    'Osc: Amplitude 4:', 'Osc: Mass 4:', 'Osc: Amplitude 5:', 'Osc: Mass 5:', 'Osc: Amplitude 6:', 'Osc: Mass 6:'}, ...
                    'Input', 1, {num2str(coefficients(1)), num2str(coefficients(2)), num2str(coefficients(3)), num2str(coefficients(4)), ...
                    num2str(coefficients(5)), num2str(coefficients(6)), num2str(coefficients(7)), num2str(coefficients(8)), num2str(coefficients(9)), num2str(coefficients(10)), num2str(coefficients(11)), num2str(coefficients(12))});
				
				if (~(size(new_guess_answers, 1) == 0))
					tmpamp1 = str2num(new_guess_answers{1});
					tmpmass1 = str2num(new_guess_answers{2});
					tmpamp2 = str2num(new_guess_answers{3});
					tmpmass2 = str2num(new_guess_answers{4});
                    tmpamp3 = str2num(new_guess_answers{5});
					tmpmass3 = str2num(new_guess_answers{6});
					tmpamp4 = str2num(new_guess_answers{7});
					tmpmass4 = str2num(new_guess_answers{8});
					tmpamp5 = str2num(new_guess_answers{9});
					tmpmass5 = str2num(new_guess_answers{10});
					tmpamp6 = str2num(new_guess_answers{11});
					tmpmass6 = str2num(new_guess_answers{12});
					
					if (~(size(tmpamp1, 1) == 0))
						coefficients(1) = tmpamp1;
					end
					
					if (~(size(tmpmass1, 1) == 0))
						coefficients(2) = tmpmass1;
					end
					
					if (~(size(tmpamp2, 1) == 0))
						coefficients(3) = tmpamp2;
					end
					
					if (~(size(tmpmass2, 1) == 0))
						coefficients(4) = tmpmass2;
                    end
                    
                    if (~(size(tmpamp3, 1) == 0))
						coefficients(5) = tmpamp3;
					end
					
					if (~(size(tmpmass3, 1) == 0))
						coefficients(6) = tmpmass3;
					end
					
					if (~(size(tmpamp4, 1) == 0))
						coefficients(7) = tmpamp4;
					end
					
					if (~(size(tmpmass4, 1) == 0))
						coefficients(8) = tmpmass4;
					end
					
					if (~(size(tmpamp5, 1) == 0))
						coefficients(9) = tmpamp5;
					end
					
					if (~(size(tmpmass5, 1) == 0))
						coefficients(10) = tmpmass5;
					end
					
					if (~(size(tmpamp6, 1) == 0))
						coefficients(11) = tmpamp6;
					end
					
					if (~(size(tmpmass6, 1) == 0))
						coefficients(12) = tmpmass6;
					end
					
					run render_update;
					
					clear('tmpamp1');
					clear('tmpamp2');
                    clear('tmpamp3');
					clear('tmpamp4');
					clear('tmpamp5');
					clear('tmpamp6');
					clear('tmpmass1');
					clear('tmpmass2');
					clear('tmpmass3');
					clear('tmpmass4');
					clear('tmpmass5');
					clear('tmpmass6');

                end
				
				clear('new_guess_answers'); 
			
			case 5 % Set Constraints
				new_guess_answers = inputdlg({'Cosh: Amplitude 1 constraint:', 'Cosh: Mass 1 constraint:', 'Cosh: Amplitude 2 constraint:', 'Cosh: Mass 2 constraint:', 'Cosh: Amplitude 3 constraint:', 'Cosh: Mass 3 constraint:',...
                    'Osc: Amplitude 4 constraint:', 'Osc: Mass 4 constraint:', 'Osc: Amplitude 5 constraint:', 'Osc: Mass 5 constraint:', 'Osc: Amplitude 6 constraint:', 'Osc: Mass 6 constraint:'}, ...
                    'Input', 1, {num2str(constraints(1)), num2str(constraints(2)), num2str(constraints(3)), num2str(constraints(4)), ...
                    num2str(constraints(5)), num2str(constraints(6)), num2str(constraints(7)), num2str(constraints(8)), num2str(constraints(9)), num2str(constraints(10)), num2str(constraints(11)), num2str(constraints(12))});
				
				if (~(size(new_guess_answers, 1) == 0))
					tmpamp1 = str2num(new_guess_answers{1});
					tmpmass1 = str2num(new_guess_answers{2});
					tmpamp2 = str2num(new_guess_answers{3});
					tmpmass2 = str2num(new_guess_answers{4});
                    tmpamp3 = str2num(new_guess_answers{5});
					tmpmass3 = str2num(new_guess_answers{6});
					tmpamp4 = str2num(new_guess_answers{7});
					tmpmass4 = str2num(new_guess_answers{8});
					tmpamp5 = str2num(new_guess_answers{9});
					tmpmass5 = str2num(new_guess_answers{10});
					tmpamp6 = str2num(new_guess_answers{11});
					tmpmass6 = str2num(new_guess_answers{12});
					
					if (~(size(tmpamp1, 1) == 0))
						constraints(1) = tmpamp1;
					end
					
					if (~(size(tmpmass1, 1) == 0))
						constraints(2) = tmpmass1;
					end
					
					if (~(size(tmpamp2, 1) == 0))
						constraints(3) = tmpamp2;
					end
					
					if (~(size(tmpmass2, 1) == 0))
						constraints(4) = tmpmass2;
                    end
                    
                    if (~(size(tmpamp3, 1) == 0))
						constraints(5) = tmpamp3;
					end
					
					if (~(size(tmpmass3, 1) == 0))
						constraints(6) = tmpmass3;
					end
					
					if (~(size(tmpamp4, 1) == 0))
						constraints(7) = tmpamp4;
					end
					
					if (~(size(tmpmass4, 1) == 0))
						constraints(8) = tmpmass4;
					end
					
					if (~(size(tmpamp5, 1) == 0))
						constraints(9) = tmpamp5;
					end
					
					if (~(size(tmpmass5, 1) == 0))
						constraints(10) = tmpmass5;
					end
					
					if (~(size(tmpamp6, 1) == 0))
						constraints(11) = tmpamp6;
					end
					
					if (~(size(tmpmass6, 1) == 0))
						constraints(12) = tmpmass6;
					end
					
					run render_update;
					
					clear('tmpamp1');
					clear('tmpamp2');
                    clear('tmpamp3');
					clear('tmpamp4');
					clear('tmpamp5');
					clear('tmpamp6');
					clear('tmpmass1');
					clear('tmpmass2');
					clear('tmpmass3');
					clear('tmpmass4');
					clear('tmpmass5');
					clear('tmpmass6');

                end
				
				clear('new_guess_answers'); 
			
				
			case 13 % Reset Params
                coefficients = zeros(12,1);
				errors = zeros(12, 1);
                chisq_dof = 0.0; p_val = 0.0; cond_num = 0;
				
				run render_update;
				
			case 14 % Reset Constraints
                constraints = zeros(12,1);
                chisq_dof = 0.0; p_val = 0.0; cond_num = 0;
				
				run render_update;
				
			case 6 % Perform Fit
			
				run prepare_data;
			
				% Whelp, here we go!
				the_fit_output = get_all_nlfit_multi(rescale_sum, rescale_cov_mat, fit_minimum, fit_maximum, parse_Nt, fit_even, mod(fit_form,4), (fit_form-mod(fit_form,4))/4, fit_diag, fit_x_prec, coefficients, constraints);
				
				clear('rescale_sum'); clear('rescale_err'); clear('rescale_cov_mat');
				
				if (~(size(the_fit_output, 1) == 0))
				
                    coefficients(:) = the_fit_output(1, 3:14);
                    errors = zeros(1, 12);

					chisq_dof = the_fit_output(1, 15);
					p_val = the_fit_output(1, 16);
					cond_num = the_fit_output(1, 17);
					
					% Push it to the save stack, overwriting if it's the same t-min.
					if (saving == 1)
						if (results_save(end, 1) == fit_minimum)
							% Override.
							results_save(end, 1) = fit_minimum;
							results_save(end, 2) = fit_maximum;
							results_save(end, 3:2:25) = coefficients(:);
							results_save(end, 4:2:26) = errors(:);
							results_save(end, 27) = chisq_dof;
							results_save(end, 28) = p_val;
							results_save(end, 29:40) = constraints(:);
						else
							% Add to the end.
							results_save(end+1, 1) = fit_minimum;
							results_save(end, 2) = fit_maximum;
							results_save(end, 3:2:25) = coefficients(:);
							results_save(end, 4:2:26) = errors(:);
							results_save(end, 27) = chisq_dof;
							results_save(end, 28) = p_val;
							results_save(end, 29:40) = constraints(:);
						end
					end
					
					run render_update
				
				else
					run render_update
				end
				
			case 7 % Perform Jackknife, finally.
			
				run prepare_data;
			
				% Whelp, here we go!
				the_fit_output = get_all_nlfit_multi(rescale_sum, rescale_cov_mat, fit_minimum, fit_maximum, parse_Nt, fit_even, mod(fit_form,4), (fit_form-mod(fit_form,4))/4, fit_diag, fit_x_prec, coefficients, constraints);
				
				
				if (~(size(the_fit_output, 1) == 0))
				
                    %coefficients(:) = the_fit_output(1, 3:14);
					
					chisq_dof = the_fit_output(1, 15);
					p_val = the_fit_output(1, 16);
					cond_num = the_fit_output(1, 17);
					
					% Into the rabbit hole...
					coefficients_blocks = zeros(num_blocks, 12);
					flag = 1; is_blocking = 1;
					for b=1:num_blocks
                        if (flag == 1)
    						block_fit_output = get_all_nlfit_multi(rescale_jack(:,b), rescale_cov_mat, fit_minimum, fit_maximum, parse_Nt, fit_even, mod(fit_form,4), (fit_form-mod(fit_form,4))/4, fit_diag, fit_x_prec, coefficients, constraints);
						
        					if (numel(block_fit_output) == 0)
            					flag = 0;
                				continue;
                            end
						
                        	coefficients_blocks(b, :) = block_fit_output(3:14);
                            
                            if (mod(b,10) == 0)
                                run render_update;
                            end
                        end
					end
					if flag == 1 % it worked!
						coefficients_rep = repmat(coefficients', [num_blocks 1]);
						errors = sqrt(sum((coefficients_blocks - coefficients_rep).^2,1).*(num_blocks-1)./num_blocks);
					else % it failed
						errors = zeros(1,12);
                    end
                    
                    is_blocking  = 0;
					
					% Push it to the save stack, overwriting if it's the same t-min.
					if (saving == 1)
						if (results_save(end, 1) == fit_minimum)
							% Override.
							results_save(end, 1) = fit_minimum;
							results_save(end, 2) = fit_maximum;
							results_save(end, 3:2:25) = coefficients(:);
							results_save(end, 4:2:26) = errors(:);
							results_save(end, 27) = chisq_dof;
							results_save(end, 28) = p_val;
							results_save(end, 29:40) = constraints(:);
						else
							% Add to the end.
							results_save(end+1, 1) = fit_minimum;
							results_save(end, 2) = fit_maximum;
							results_save(end, 3:2:25) = coefficients(:);
							results_save(end, 4:2:26) = errors(:);
							results_save(end, 27) = chisq_dof;
							results_save(end, 28) = p_val;
							results_save(end, 29:40) = constraints(:);
						end
					end
					
					run render_update
				
				else
					errors = zeros(12, 1);
					run render_update
				end
			
			
				clear('rescale_sum'); clear('rescale_err'); clear('rescale_cov_mat');
				
			
			case 10 % Change Directory
				directory_tmp = inputdlg('Enter an input directory.', 'Input', 1, directory_answer);
				
                
				if (~(size(directory_tmp,1) == 0) && ~(exist(strcat('..\..\..\', directory_tmp{1}),'dir') == 0))
					directory_answer = directory_tmp;
					
					% Load ensemble info.
					%[fl_l, fl_h, parse_Ns, parse_Nt, beta, m_l, m_h] = load_info_file(fname);
					[fl_l, fl_h, parse_Ns, parse_Nt, beta, m_l, m_h] = load_info_file(strcat(['..\..\..\' directory_answer{1}]));
					
					% Get all files in said directory to see if we can match
					% the current answer.
					match = 0;
					cnt = 0;
					files_in_dir = dir(strcat('..\..\..\', directory_answer{1}, '\spectrum2\corr\'));
					spectrum_new_list = {};
					for i=1:size(files_in_dir,1)
						[tmppath, tmpname, tmpext] = fileparts(files_in_dir(i).name);
						if (strcmp(tmpname, 'corr')) % if the file has the right name
							cnt = cnt+1;
							spectrum_new_list{size(spectrum_new_list,2)+1} = tmpext(2:end);
							if (strcmp(tmpext(2:end), spectrum_text))
								match = cnt;
							end
						end
                    end
                    
                    % If a spectrum.dat and a pbppart.dat file exist, add "dc_stoch" and
                    % "sg_stoch" to the list of states.
                    if (exist(strcat(['..\..\..\' directory_answer{1} '\spectrum2\stoch\PBPPART.dat']), 'file') && exist(strcat(['..\..\..\' directory_answer{1} '\spectrum2\stoch\SPECTRUM.dat']), 'file'))
                        % We can build dc_stoch and sg_stoch
                        spectrum_new_list{size(spectrum_list,2)+1} = 'dc_stoch';
                        spectrum_new_list{size(spectrum_list,2)+1} = 'sg_stoch';
                    end


					clear('files_in_dir'); clear('tmppath'); clear('tmpname'); clear('tmpext');

					spectrum_list = spectrum_new_list;
					
					if (match == 0) % no match

						% Once the directory exists, find which input code to look at.
						%spectrum_list = {'ps', 'sc', 'i5', 'ij', 'r0', 'ris', 'rij', 'ri5', 'ps2', 'nu', 'de', 'sg', 'dc', 'sgv', 'dcv'};
						spectrum_answer = listdlg('PromptString','Original state does not exist. Enter a spectrum to look at.', ...
							'SelectionMode', 'single', 'ListString', spectrum_list);

						if (isempty(spectrum_answer))
							return;
						end

						spectrum_text = spectrum_list{spectrum_answer};
					else
						spectrum_answer = match;
                    end
                    
					
					% Start loading the data!
					blocksize = 0;
					
					% If we didn't get range arguments...
					if (~exist('blockval', 'var'))
						% See if a fitparams file exists.
						if (exist(fullfile(cd, strcat(['..\..\..\' directory_answer{1} '\spectrum2\fitparams\fitparam.' spectrum_text])), 'file'))
							fitparam = importdata(strcat(['..\..\..\' directory_answer{1} '\spectrum2\fitparams\fitparam.' spectrum_text]));
							blocksize = fitparam(5);
							
						else
							error('A fitparams file is expected but does not exist!');
						end
					else
						blocksize = blockval;
					end
						
					mass_l = m_l;

					number_files = 1;
					%num_data = 600;
					%number_bl = 4;
					%parse_Nt = 48;
					%blocksize = 10;
					volume = parse_Ns^3; %*48;
					
					% Import data.
					[scalar_sum, scalar_jack, scalar_cov_mat, ...
						scalar_err, num_blocks] = ...
						get_correlator(strcat(['../../../' directory_answer{1}]), ...
						spectrum_text, parse_Nt, parse_Ns, fl_flavor, blocksize);
						
					
                    run render_update;
					
					saving = 0;
					results_save = zeros(1, 40);
					
				end
			
			case 11 % Change State
				%spectrum_list = {'ps', 'sc', 'i5', 'ij', 'r0', 'ris', 'rij', 'ri5', 'ps2', 'nu', 'de'};
				
				spectrum_guess = listdlg('PromptString','Enter a spectrum to look at.','SelectionMode', 'single', 'ListString', spectrum_list, 'InitialValue', spectrum_answer);

				if (~isempty(spectrum_guess))
				
					spectrum_answer = spectrum_guess;
					spectrum_text = spectrum_list{spectrum_answer};
				
					
					% Start loading the data!
					blocksize = 0;
					
					% If we didn't get range arguments...
					if (~exist('blockval', 'var'))
						% See if a fitparams file exists.
						if (exist(fullfile(cd, strcat(['..\..\..\' directory_answer{1} '\spectrum2\fitparams\fitparam.' spectrum_text])), 'file'))
							fitparam = importdata(strcat(['..\..\..\' directory_answer{1} '\spectrum2\fitparams\fitparam.' spectrum_text]));
							blocksize = fitparam(5);
							
						else
							error('A fitparams file is expected but does not exist!');
						end
					else
						blocksize = blockval;
					end
						
					mass_l = m_l;

					number_files = 1;
					%num_data = 600;
					%number_bl = 4;
					%parse_Nt = 48;
					%blocksize = 10;
					volume = parse_Ns^3; %*48;
					
					% Import data.
					[scalar_sum, scalar_jack, scalar_cov_mat, ...
						scalar_err, num_blocks] = ...
						get_correlator(strcat(['../../../' directory_answer{1}]), ...
						spectrum_text, parse_Nt, parse_Ns, fl_flavor, blocksize);
		
					saving = 0;
					results_save = zeros(1, 40);
					
					run render_update;
					
				end
				
				clear('spectrum_guess');
				
			case 8 % Begin saving.
				if (saving == 0)
					saving = 1;
					results_save = zeros(1,40);
					run render_update;
				end
			case 9 % End saving.
				if (saving == 1)
					saving = 0;
					
					% Get a filename!
					[temp_fname, temp_directory] = uiputfile('*.*', 'Save file...', strcat(['../../../' directory_answer{1} '/spectrum2/multifits/']));
					
					if (temp_fname ~= 0)
                        save(strcat([temp_directory, temp_fname]), 'results_save', '-ascii');
						results_save = zeros(1,40);
                    end
					
					run render_update;
				end
			case 15
				flag = 0;
		end
		
	end

	close(h);
	close(h2);

	return; 

end