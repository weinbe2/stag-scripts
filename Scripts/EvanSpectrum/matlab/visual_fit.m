function visual_fit(blockval, num_elim)
	% A windowed-ish fitting program!

	% Load some useful directories!
	addpath('.\process', '-end');
	addpath('.\multifit', '-end');
	
	% Set num_elim to 1 if it's unset.
	if (~exist('num_elim', 'var'))
		num_elim = 1; % single eliminate
	end
	
	% First, we pick a directory.


	directory_answer = inputdlg('Enter an input directory.', 'Input', 1, {'FourPlusEight\f4plus8l36t64b40m003m080'});
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
		scalar_err, num_blocks, scalar_jack_single] = ...
		get_correlator(strcat(['../../../' directory_answer{1}]), ...
		spectrum_text, parse_Nt, parse_Ns, fl_flavor, blocksize, num_elim);
	
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
	log_answer = 0; % plot on a log scale. 
	fit_form = 1; % One cosh. (2 = Two cosh, 3 = One cosh, one osc, 4 = Baryon, 5 = Cosh + const, 6 = cosh with zero shift).

	run set_visual_defaults;
	
	% set values for effective masses.
	
	eff_K = 2; 
	eff_N = 5; % minimum 2*n, maximum... Nt.
	eff_C = 0; % number to cut. Must be less than K.
	
	% See if it's a baryon...
	if (strcmp(spectrum_text, 'nu') || strcmp(spectrum_text, 'de'))
		is_baryon = 1;
	else
		is_baryon = 0;
	end
	
	% Get correct cosh functions in place.
	run get_cosh; 
	
	% Build plots strings as needed.
	plot_str{1} = '+r';
	plot_str{2} = 'og';
	plot_str{3} = 'xb';
	plot_str{4} = 'sc';
	plot_str{5} = 'dm';
	plot_str{6} = 'pk';
	
	% Prepare figures.
	h = figure();
	movegui(h, 'northwest');
	set(h, 'Name', 'Correlator');
	errorbar(xrange, scalar_sum, scalar_err, '.k');  hold on;
	[funcdata] = get_function(xrange, parse_Nt, fit_form, coefficients);
	plot(xrange, funcdata); hold off;
	axis([-1 parse_Nt -max(abs(scalar_sum))*1.1 max(abs(scalar_sum))*1.1]);

	h2 = figure('menu', 'none', 'toolbar', 'none', 'units', 'normalized', 'position', [0 0 .3 .9]);
	movegui(h2, 'northeast');
	set (h2, 'Name', 'Info');
	ph = uipanel(h2, 'Units', 'normalized', 'position', [0.05 0.05 0.9 0.9], 'title', 'Display Window');
	th = uicontrol(ph, 'style','text','Units','normalized','position',[0 0 1 1],'FontSize', 9, 'string', {'Rescale mass: 0.0', 'Fit Form: One cosh', 'Amplitude 1: 0.0', 'Mass 1: 0.0', 'Amplitude 2: 0.0', 'Mass 2: 0.0', 'Chisq per DoF: 0.0'}, 'horizontalalignment', 'left', 'FontName', 'FixedWidth');

	% Get the correct cosh functions in place.
	run get_cosh; 
	
	run render_update;

	figure(h);


	% Enter the main loop of options.

	flag = 1;

	while (flag == 1)

		option_list = {'Rescale+Shift (Visual Only)', ...
						'Set Fit Form', ...
						'Set Fit Options', ...
						'Set Data Modifications', ...
						'Set Binning (Reload)', ...
						'Set Initial Guess', ...
						'Set Constraints', ...
						'Perform Fit', ...
						'Perform Jackknife', ...
						'Save Modified Correlator', ...
						'Begin Saving Results', ...
						'End Saving Results', ...
						'Change Directory', ...
						'Change State', ...
						'Reset Fit Options', ...
						'Reset Initial Guess', ...
						'Reset Constraints', ...
						'Visualize Singular Values', ...
						'Visualize Effective Masses', ...
						'Save Effective Masses', ...
                        'Exit'};
		option_answer = listdlg('PromptString','Choose an option.','SelectionMode', 'single', 'ListString', option_list);
		
		if (isempty(option_answer))
			option_answer = 21;
		end
		
		switch option_answer
			case 1 % Rescale (Purely Visual Effect)
				
				new_mass_answer = inputdlg({'Enter a mass to rescale by.','Enter a vev to subtract.', 'Plot on a log (0 no, 1 yes).'}, 'Input', 1, {num2str(mass_answer),num2str(vev_answer), num2str(log_answer)});
				
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
					
					tmplog = str2num(new_mass_answer{3});
					if (~(size(tmplog, 1) == 0))
						% Rescale!
						log_answer = tmplog;
						
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
				% Jackknife: number to eliminate on jackknife?
			
				new_fit_min = inputdlg({'Minimum Fit t:', ...
										'Maximum Fit t (-1 for symmetric fit):', ...
										'Fit X Precion:', ...
										'Fit Every Other (0 no, 1 yes)', ...
										'Cut Singular Values', ...
										'Diagonal Correlator (0 no, 1 yes)', ...
										'Zeroed Cosh (0 no, 1 yes)', ...
										'Finite Difference Cosh (0 no, 1 yes)'}, ...
										'Input', 1, {num2str(fit_minimum), num2str(fit_maximum), num2str(fit_x_prec),  num2str(fit_even), num2str(fit_cut), num2str(fit_diag), num2str(func_zshift),  num2str(func_diff)});
				
				
				
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
					
					tmpfiteven = str2num(new_fit_min{4});
					if (~(size(tmpfiteven, 1) == 0))
						if (tmpfiteven == 0 || tmpfiteven == 1)
							fit_even = tmpfiteven;
						end
					end
					
					tmpfitsvd = str2num(new_fit_min{5});
					if (~(size(tmpfitsvd, 1) == 0))
						if (tmpfitsvd >= 0 && tmpfitsvd <= parse_Nt)
							fit_cut = tmpfitsvd;
						end
					end
					
					tmpfitdiag = str2num(new_fit_min{6});
					if (~(size(tmpfitdiag, 1) == 0))
						if (tmpfitdiag == 0 || tmpfitdiag == 1)
							fit_diag = tmpfitdiag;
						end
					end
					
					tmpfunczshift = str2num(new_fit_min{7});
					if (~(size(tmpfunczshift, 1) == 0))
						if (tmpfunczshift == 0 || tmpfunczshift == 1)
							func_zshift = tmpfunczshift;
						end
					end
					
					tmpfuncdiff = str2num(new_fit_min{8});
					if (~(size(tmpfuncdiff, 1) == 0))
						if (tmpfuncdiff == 0 || tmpfuncdiff == 1)
							func_diff = tmpfuncdiff;
						end
					end
				end
				
				clear('tmpfitmin');
				clear('tmpfitmax');
				clear('tmpfiteven');
				clear('tmpfitdiag');
				clear('tmpfitsvd');
				clear('tmpxprec');
				clear('tmpfunczshift');
				clear('tmpfuncdiff');
				
				% Fix the cosh.
				run get_cosh;
				
				run render_update;
				
			case 4 % Set Data Manipulations
				
				% Options: tmin, tmax
				% Fold? y/n
				% Positive Parity Project? y/n
				% Fit to only even time data? y/n
				% Full or diagonal correlator matrix? 
				% Jackknife: number to eliminate on jackknife?
			
				new_fit_min = inputdlg({'Fold (0 no, 1 yes):', 'Parity Project (0 no, 1 positive, -1 negative):', 'Zero Center (0 no, 1 yes, -1 gap)', ...
										'Finite Difference Data (0 no, 1 yes)'}, ...
										'Input', 1, {num2str(func_fold), num2str(fit_ppp), num2str(fit_zero),   num2str(fit_diff)});
				
				
				
				if (~(size(new_fit_min, 1) == 0)) % Make sure we didn't get a cancel!
					
					tmpfuncfold = str2num(new_fit_min{1});
					if (~(size(tmpfuncfold, 1) == 0))
						if (tmpfuncfold == 0 || tmpfuncfold == 1)
							func_fold = tmpfuncfold;
						end
					end
					
					tmpfitppp = str2num(new_fit_min{2});
					if (~(size(tmpfitppp, 1) == 0))
						if (tmpfitppp == 0 || tmpfitppp == 1 || tmpfitppp == -1)
							fit_ppp = tmpfitppp;
						end
					end
					
					tmpfitzero = str2num(new_fit_min{3});
					if (~(size(tmpfitzero, 1) == 0))
						if (tmpfitzero == 0 || tmpfitzero == 1 || tmpfitzero == -1)
							fit_zero = tmpfitzero;
						end
					end
					
					tmpfitdiff = str2num(new_fit_min{4});
					if (~(size(tmpfitdiff, 1) == 0))
						if (tmpfitdiff == 0 || tmpfitdiff == 1)
							fit_diff = tmpfitdiff;
						end
					end
					
				end
				
				
				clear('tmpfitfold');
				clear('tmpfitppp');
				clear('tmpfitzero');
				clear('tmpfitdiff');
				
				% Fix the cosh.
				run get_cosh;
				
				run render_update;
				
			case 15 % Reset Fit Options
				
				% Ask first!
				choice = questdlg('Are you sure you want to reset the fit options?', 'Check!', ...
					'Yes', 'No', 'Yes');
				
				if (strcmp(choice, 'Yes'))
					
					run set_visual_defaults;
					
				end
				
				run get_cosh; 
				
				run render_update;
			
			case 5 % Set binning/jackknife. (requires reload).
				
				% Options: bin size, jackknife.
			
				new_fit_min = inputdlg({'Bin size:', 'Jackknife elim:'}, ...
										'Input', 1, {num2str(blocksize), num2str(num_elim)});
				
				
				
				if (~(size(new_fit_min, 1) == 0)) % Make sure we didn't get a cancel!
					
					tmpbinsize = str2num(new_fit_min{1});
					if (~(size(tmpbinsize, 1) == 0))
						if (tmpbinsize > 0)
							blocksize = tmpbinsize;
						end
					end
					
					tmpnumelim = str2num(new_fit_min{2});
					if (~(size(tmpnumelim, 1) == 0))
						if (tmpnumelim > 0)
							num_elim = tmpnumelim;
						end
					end
					
					% Import data.
					[scalar_sum, scalar_jack, scalar_cov_mat, ...
						scalar_err, num_blocks, scalar_jack_single] = ...
						get_correlator(strcat(['../../../' directory_answer{1}]), ...
						spectrum_text, parse_Nt, parse_Ns, fl_flavor, blocksize, num_elim);
					
					saving = 0;
					results_save = zeros(1, 40);
					
				end
				
				
				clear('tmpbinsize');
				clear('tmpnumelim');
				
				% Fix the cosh.
				run get_cosh;
				
				run render_update;
			
			
			case 6 % Set Initial Params
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
			
			case 7 % Set Constraints
			
				% Note! If the constraint is >1e-22, <1e-18,
				% it is set as an exact equality (basically...) 
				% and doesn't get counted in the dof!
			
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
			
				
			case 16 % Reset Params
			
				% Ask first!
				choice = questdlg('Are you sure you want to reset the fit params?', 'Check!', ...
					'Yes', 'No', 'Yes');
				
				if (strcmp(choice, 'Yes'))
					coefficients = zeros(12,1);
					errors = zeros(12, 1);
					chisq_dof = 0.0; p_val = 0.0; cond_num = 0;
				end
				
				run render_update;
				
			case 17 % Reset Constraints
                constraints = zeros(12,1);
                chisq_dof = 0.0; p_val = 0.0; cond_num = 0;
				
				run render_update;
				
			case 8 % Perform Fit
			
				run prepare_data;
			
				% Whelp, here we go!
				the_fit_output = get_all_nlfit_multi(rescale_sum, rescale_cov_mat, fit_minimum, fit_maximum, parse_Nt, fit_even, fit_cut, func_diff, mod(fit_form,4), (fit_form-mod(fit_form,4))/4, fit_diag, fit_x_prec, coefficients, constraints);
				
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
				
			case 9 % Perform Jackknife, finally.
			
				run prepare_data;
			
				% Whelp, here we go!
				the_fit_output = get_all_nlfit_multi(rescale_sum, rescale_cov_mat, fit_minimum, fit_maximum, parse_Nt, fit_even, fit_cut, func_diff, mod(fit_form,4), (fit_form-mod(fit_form,4))/4, fit_diag, fit_x_prec, coefficients, constraints);
				
				
				if (~(size(the_fit_output, 1) == 0))
				
                    %coefficients(:) = the_fit_output(1, 3:14);
					
					chisq_dof = the_fit_output(1, 15);
					p_val = the_fit_output(1, 16);
					cond_num = the_fit_output(1, 17);
					
					% Into the rabbit hole...
					coefficients_blocks = zeros(num_blocks, 12);
					jack_flag = 1; is_blocking = 1;
					for b=1:num_blocks
                        if (jack_flag == 1)
    						block_fit_output = get_all_nlfit_multi(rescale_jack(:,b), rescale_cov_mat, fit_minimum, fit_maximum, parse_Nt, fit_even, fit_cut, func_diff, mod(fit_form,4), (fit_form-mod(fit_form,4))/4, fit_diag, fit_x_prec, coefficients, constraints);
						
        					if (numel(block_fit_output) == 0)
            					jack_flag = 0;
                				continue;
                            end
						
                        	coefficients_blocks(b, :) = block_fit_output(3:14);
                            
                            if (mod(b,10) == 0)
                                run render_update;
                            end
                        end
					end
					if jack_flag == 1 % it worked!
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
				
			
			case 13 % Change Directory
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
						cnt = cnt+1;
                        spectrum_new_list{size(spectrum_new_list,2)+1} = 'dc_stoch';
						if (strcmp('dc_stoch', spectrum_text))
							match = cnt;
						end
						
						cnt = cnt+1;
                        spectrum_new_list{size(spectrum_new_list,2)+1} = 'sg_stoch';
						if (strcmp('sg_stoch', spectrum_text))
							match = cnt;
						end
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
					fl_flavor = fl_l/4;
					
					xrange = (0:(parse_Nt-1))';
					
					% Import data.
					[scalar_sum, scalar_jack, scalar_cov_mat, ...
						scalar_err, num_blocks, scalar_jack_single] = ...
						get_correlator(strcat(['../../../' directory_answer{1}]), ...
						spectrum_text, parse_Nt, parse_Ns, fl_flavor, blocksize, num_elim);
						
					% See if it's a baryon...
					if (strcmp(spectrum_text, 'nu') || strcmp(spectrum_text, 'de'))
						is_baryon = 1;
					else
						is_baryon = 0;
					end
					
					run get_cosh; 
					
                    run render_update;
					
					saving = 0;
					results_save = zeros(1, 40);
					
				end
			
			case 14 % Change State
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
					fl_flavor = fl_l/4;
					xrange = (0:(parse_Nt-1))';
					
					% Import data.
					[scalar_sum, scalar_jack, scalar_cov_mat, ...
						scalar_err, num_blocks, scalar_jack_single] = ...
						get_correlator(strcat(['../../../' directory_answer{1}]), ...
						spectrum_text, parse_Nt, parse_Ns, fl_flavor, blocksize, num_elim);
		
					saving = 0;
					results_save = zeros(1, 40);
					
					% See if it's a baryon...
					if (strcmp(spectrum_text, 'nu') || strcmp(spectrum_text, 'de'))
						is_baryon = 1;
					else
						is_baryon = 0;
					end
					
					run get_cosh; 
					
					run render_update;
					
				end
				
				clear('spectrum_guess');
				
			case 11 % Begin saving.
				if (saving == 0)
					saving = 1;
					results_save = zeros(1,40);
					run render_update;
				end
			case 12 % End saving.
				if (saving == 1)
					saving = 0;
					
					% Get a filename!
					[temp_fname, temp_directory] = uiputfile('*.*', 'Save file...', strcat(['../../../' directory_answer{1} '/spectrum2/multifits/']));
					
					if (temp_fname ~= 0)
                        save(strcat([temp_directory, temp_fname]), 'results_save', '-ascii', '-double');
						results_save = zeros(1,40);
                    end
					
					run render_update;
				end
			case 10 % Save modified correlator 
			
				% get the data in shape!
				run render_update;
				
				% Prepare data it!
				save_data = zeros(parse_Nt,3);
				save_data(:,1) = 0:(parse_Nt-1);
				if (fit_diff == 1) % shift forward by 0.5
					save_data(:,1) = save_data(:,1) + 0.5;
				end
				save_data(:,2) = rescale_sum(:);
				save_data(:,3) = rescale_err(:);
				
				% First, ask for directory to save the correlator.
				[temp_fname, temp_directory] = uiputfile('*.*', 'Save modified correlator...', strcat(['../../../' directory_answer{1} '/spectrum2/sum/']));
				
				if (temp_fname ~= 0)
					save(strcat([temp_directory, temp_fname]), 'save_data', '-ascii', '-double');
				end
				
				% Next, ask for a place to save the covariance.
				[temp_fname, temp_directory] = uiputfile('*.*', 'Save modified variance-covariance...', strcat(['../../../' directory_answer{1} '/spectrum2/sum/']));
				
				if (temp_fname ~= 0)
					save(strcat([temp_directory, temp_fname]), 'rescale_cov_mat', '-ascii', '-double');
				end
				
			case 18 % visualize singular values.
				% Get the appropriate covariance matrix.
				
				t1 = fit_minimum;
				t2 = fit_maximum;
				
				if (fit_even == 0) % fit to all
					if (func_diff == 0) % still all
						xval = t1:t2;
						yval = rescale_sum((t1+1):(t2+1))';
						ycorr = rescale_cov_mat((t1+1):(t2+1),(t1+1):(t2+1));
					else % finite diff! one less data.
						xval = t1:(t2-1);
						yval = rescale_sum((t1+1):t2)';
						ycorr = rescale_cov_mat((t1+1):t2,(t1+1):t2);
					end
				else % only fit to even...
					if (func_diff == 0) % no finite diff.
						if (mod(t1,2)==0)
							xval = t1:2:t2;
						else
							xval = (t1+1):2:t2;
						end
						yval = rescale_sum((xval+1))';
						ycorr = rescale_cov_mat((xval+1),(xval+1));
					else % finite diff! careful about end.
						if (mod(t1,2)==0)
							xval = t1:2:(t2-1);
						else
							xval = (t1+1):2:(t2-1);
						end
						yval = rescale_sum((xval+1))';
						ycorr = rescale_cov_mat((xval+1),(xval+1));
					end
                end
				
                diagvals = sort(diag(ycorr)); % get the diagonal values.
				svals = fliplr(svd(ycorr)'); % minimum first.
				
				
				h_svd = figure();
				%movegui(h, 'northwest');
				set(h_svd, 'Name', 'SVD visualization');
				h_svals = semilogy(svals,'xb'); hold on;
                h_dvals = semilogy(diagvals, 'or'); hold off;
                legend([h_svals, h_dvals],'Singular Values', 'Diagonal Elements', 'Location', 'southeast');
				waitfor(h_svd);
				
				% clear everything.
				clear('h_svd');
                clear('h_svals');
                clear('h_dvals');
				clear('svals');
				clear('xval');
				clear('yval');
				clear('ycorr');
				clear('t1');
				clear('t2');
			
			case 19 % Visualize Effective Masses
			
				% Ask what values of K, N, and C to use.
				% eff_K = 2; 
				% eff_N = 5; % minimum 2*n, maximum... Nt.
				% eff_C = 0; % number to cut. Must be less than K.
				
				% Ask for values then check sanity.
				
				are_valid = 0;
				while (are_valid == 0)
					are_valid = 1;
					new_effmass = inputdlg({'Number of states (1 to 6)', 'Values to use (2*N_states to N_t/2)', 'States to SVD Cut (0 to N_states)'}, ...
										'Input', 1, {num2str(eff_K), num2str(eff_N), num2str(eff_C)});
				
				
				
					if (~(size(new_effmass, 1) == 0)) % Make sure we didn't get a cancel!
						
						tmpeffK = str2num(new_effmass{1});
						if (~(size(tmpeffK, 1) == 0))
							if (tmpeffK > 0 && tmpeffK < 6)
								eff_K = tmpeffK;
							else
								are_valid = 0;
							end
						end
						
						tmpeffN = str2num(new_effmass{2});
						if (~(size(tmpeffN, 1) == 0))
							if (tmpeffN >= (2*eff_K) && tmpeffN <= (parse_Nt/2))
								eff_N = tmpeffN;
							else
								are_valid = 0;
							end
						end
						
						tmpeffC = str2num(new_effmass{3});
						if (~(size(tmpeffC, 1) == 0))
							if (tmpeffC >= 0 && tmpeffC < eff_K)
								eff_C = tmpeffC;
							else
								are_valid = 0;
							end
						end
						
						if (are_valid == 0)
							msgbox('Parameter values do not match constraints.', 'Error', 'error');
						end
						
						clear('tmpeffK');
						clear('tmpeffN');
						clear('tmpeffC');
					else
						are_valid = -1; % exit
					end
				end
				
				if (are_valid == 1)
				
					% Compute masses and get errors under jackknife!
					[masses, tmproot, tmpamps] = effective_mass_utility(rescale_sum, parse_Nt, eff_K, eff_N, eff_C);
					
					% parse_Nt, num_blocks, which mass.
					% The weird order is because 
					% errors_jackknife expects (timespan, jack).
					masses_jack = zeros(size(masses,1), num_blocks, eff_K);
                    
                    % Keep things updated!
                    is_blocking = 1;
                    
                    for (b=1:num_blocks)
						[masses_jack(:,b,:), tmproot, tmpamps] = effective_mass_utility(rescale_jack(:,b), parse_Nt, eff_K, eff_N, eff_C);
                        
                        if (mod(b,10) == 0)
                            run render_update;
                        end
                    end
                    
                    is_blocking = 0;
                    
					% Save some memory...
					clear('tmproot'); clear('tmpamps');
					
					masses_err = zeros(size(masses,1), eff_K);
					
					for b=1:eff_K
						[tmpthing, masses_err(:,b)] = errors_jackknife(masses(:,b), masses_jack(:,:,b));
					end
					
					clear('tmpthing');
					
                    % Clean up masses!
                    bool_flag = ((abs(imag(masses)) > 1e-8) | (real(masses_err) > 1));
                    masses(bool_flag) = real(masses(bool_flag)) - 1e10;
                    
                    % Plot it!                    
                    h_svd = figure();
                    %movegui(h, 'northwest');
                    set(h_svd, 'Name', 'Effective mass visualization');
                    x_eff = ((eff_N-1)/2):(size(masses,1)-1+(eff_N-1)/2);
                    errorbar(x_eff-0.2, masses(:,1), masses_err(:,1), plot_str{1}); hold on;
                    if (eff_K > 1)
                        for b=2:eff_K
                            errorbar(x_eff+(0.4/(eff_K-1))*(b-1)-0.2, masses(:,b), masses_err(:,b), plot_str{b}); 
                        end
                    end
                    axis([-inf, inf, 0.000001, 1]);
                    %h_dvals = semilogy(diagvals, 'or'); hold off;
                    %legend([h_svals, h_dvals],'Singular Values', 'Diagonal Elements', 'Location', 'southeast');
                    waitfor(h_svd);
                    hold off;
                    
				end
							
			case 20 % Save Effective Masses
                    % Ask what values of K, N, and C to use.
				% eff_K = 2; 
				% eff_N = 5; % minimum 2*n, maximum... Nt.
				% eff_C = 0; % number to cut. Must be less than K.
				
				% Ask for values then check sanity.
				
				are_valid = 0;
				while (are_valid == 0)
					are_valid = 1;
					new_effmass = inputdlg({'Number of states (1 to 6)', 'Values to use (2*N_states to N_t/2)', 'States to SVD Cut (0 to N_states)'}, ...
										'Input', 1, {num2str(eff_K), num2str(eff_N), num2str(eff_C)});
				
				
				
					if (~(size(new_effmass, 1) == 0)) % Make sure we didn't get a cancel!
						
						tmpeffK = str2num(new_effmass{1});
						if (~(size(tmpeffK, 1) == 0))
							if (tmpeffK > 0 && tmpeffK < 6)
								eff_K = tmpeffK;
							else
								are_valid = 0;
							end
						end
						
						tmpeffN = str2num(new_effmass{2});
						if (~(size(tmpeffN, 1) == 0))
							if (tmpeffN >= (2*eff_K) && tmpeffN <= (parse_Nt/2))
								eff_N = tmpeffN;
							else
								are_valid = 0;
							end
						end
						
						tmpeffC = str2num(new_effmass{3});
						if (~(size(tmpeffC, 1) == 0))
							if (tmpeffC >= 0 && tmpeffC < eff_K)
								eff_C = tmpeffC;
							else
								are_valid = 0;
							end
						end
						
						if (are_valid == 0)
							msgbox('Parameter values do not match constraints.', 'Error', 'error');
						end
						
						clear('tmpeffK');
						clear('tmpeffN');
						clear('tmpeffC');
					else
						are_valid = -1; % exit
					end
				end
				
				if (are_valid == 1)
				
					% Compute masses and get errors under jackknife!
					[masses, tmproot, tmpamps] = effective_mass_utility(rescale_sum, parse_Nt, eff_K, eff_N, eff_C);
					
					% parse_Nt, num_blocks, which mass.
					% The weird order is because 
					% errors_jackknife expects (timespan, jack).
					masses_jack = zeros(size(masses,1), num_blocks, eff_K);
                    
                    % Keep things updated!
                    is_blocking = 1;
                    
                    for (b=1:num_blocks)
						[masses_jack(:,b,:), tmproot, tmpamps] = effective_mass_utility(rescale_jack(:,b), parse_Nt, eff_K, eff_N, eff_C);
                        
                        if (mod(b,10) == 0)
                            run render_update;
                        end
                    end
                    
                    is_blocking = 0;
                    
					% Save some memory...
					clear('tmproot'); clear('tmpamps');
					
					masses_err = zeros(size(masses,1), eff_K);
					
					for b=1:eff_K
						[tmpthing, masses_err(:,b)] = errors_jackknife(masses(:,b), masses_jack(:,:,b));
					end
					
					clear('tmpthing');
					
                    % Clean up masses!
                    %bool_flag = ((abs(imag(masses)) > 1e-8) | (real(masses_err) > 1));
                    %masses(bool_flag) = real(masses(bool_flag)) - 1e10;
                    
					% Package data up!
                    x_eff = ((eff_N-1)/2):(size(masses,1)-1+(eff_N-1)/2);
					
					% time + (masses+errs for each).
					save_data = zeros(size(masses, 1), 1+eff_K*2);
					save_data(:,1) = x_eff(:);
					for b=1:eff_K
						save_data(:,2*b) = real(masses(:, b));
						save_data(:,2*b+1) = real(masses_err(:,b));
					end
					
				
					% Ask for directory to save the effective masses.
					[temp_fname, temp_directory] = uiputfile('*.*', 'Save effective masses...', strcat(['../../../' directory_answer{1} '/spectrum2/effmass/']));
					
					if (temp_fname ~= 0)
						save(strcat([temp_directory, temp_fname]), 'save_data', '-ascii', '-double');
					end
				
				end
				
                    
			case 21
				flag = 0;
		end
		
	end

	close(h);
	close(h2);

	return; 

end