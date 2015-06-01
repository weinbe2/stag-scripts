function [fit_output] = get_all_nlfit_multi(corr_fcn, corr_mat, tmin, tmax, nt, fiteven, fitcosh, fitoscil, fitdiag, fit_xprec, coeff, constraints)
    % Fit some variable number of cosh or oscillating terms, up to 3 of each!
	fitfunc = fitcosh+4*fitoscil;
	
    fit_output = zeros(2, 17);

    count = 1;

    t1 = tmin;
	t2 = tmax;
	
	if (fiteven == 0) % fit to all
		xval = t1:t2;
		yval = corr_fcn((t1+1):(t2+1))';
		ycorr = corr_mat((t1+1):(t2+1),(t1+1):(t2+1));
	else % only fit to even...
		if (mod(t1,2)==0)
			xval = t1:2:t2;
		else
			xval = (t1+1):2:t2;
		end
		yval = corr_fcn((xval+1))';
		ycorr = corr_mat((xval+1),(xval+1));
	end
	
	
	if (fitdiag == 1) % if we only fit to the diagonal corr
		ycorr = diag(diag(ycorr));
	end
	
	fit_output(count, 17) = cond(ycorr);
    
    
    choice = 'Yes';
    
    if (fit_output(count, 17) > 1e13)
        choice = 'No';
        %choice = questdlg(strcat(['Condition number = ' num2str(fit_output(count, 8)) '. Fit may crash MATLAB! Try anyway?']), 'Warning!');
    end
	
	% Convert the coefficients on the fly!
	coeff(6) = -log(coeff(6) - coeff(4));
	coeff(4) = -log(coeff(4) - coeff(2));
	coeff(2) = -log(coeff(2));
	
	coeff(12) = -log(coeff(8) - coeff(12));
	coeff(10) = -log(coeff(10) - coeff(8));
	coeff(8) = -log(coeff(8));
	
    
    if (strcmp(choice, 'Yes'))
	
        switch fitfunc
         case 1 % cosh 1 oscil 0
            guess = [coeff(1) coeff(2) ];
            chisqfunc = @(x)(1.0/(size(yval,2)-numel(guess))*(yval(:)-(+x(1)*cosh((+exp(-x(2)))*(nt/2-xval(:)))))'*(ycorr\(yval(:)-(+x(1)*cosh((+exp(-x(2)))*(nt/2-xval(:)))))));
         case 2 % cosh 2 oscil 0
            guess = [coeff(1) coeff(2) coeff(3) coeff(4) ];
            chisqfunc = @(x)(1.0/(size(yval,2)-numel(guess))*(yval(:)-(+x(1)*cosh((+exp(-x(2)))*(nt/2-xval(:)))+x(3)*cosh((+exp(-x(2))+exp(-x(4)))*(nt/2-xval(:)))))'*(ycorr\(yval(:)-(+x(1)*cosh((+exp(-x(2)))*(nt/2-xval(:)))+x(3)*cosh((+exp(-x(2))+exp(-x(4)))*(nt/2-xval(:)))))));
         case 3 % cosh 3 oscil 0
            guess = [coeff(1) coeff(2) coeff(3) coeff(4) coeff(5) coeff(6) ];
            chisqfunc = @(x)(1.0/(size(yval,2)-numel(guess))*(yval(:)-(+x(1)*cosh((+exp(-x(2)))*(nt/2-xval(:)))+x(3)*cosh((+exp(-x(2))+exp(-x(4)))*(nt/2-xval(:)))+x(5)*cosh((+exp(-x(2))+exp(-x(4))+exp(-x(6)))*(nt/2-xval(:)))))'*(ycorr\(yval(:)-(+x(1)*cosh((+exp(-x(2)))*(nt/2-xval(:)))+x(3)*cosh((+exp(-x(2))+exp(-x(4)))*(nt/2-xval(:)))+x(5)*cosh((+exp(-x(2))+exp(-x(4))+exp(-x(6)))*(nt/2-xval(:)))))));
         case 4 % cosh 0 oscil 1
            guess = [coeff(7) coeff(8) ];
            chisqfunc = @(x)(1.0/(size(yval,2)-numel(guess))*(yval(:)-(+x(1)*(1-2*mod(xval(:),2)).*cosh((+exp(-x(2)))*(nt/2-xval(:)))))'*(ycorr\(yval(:)-(+x(1)*(1-2*mod(xval(:),2)).*cosh((+exp(-x(2)))*(nt/2-xval(:)))))));
         case 5 % cosh 1 oscil 1
            guess = [coeff(1) coeff(2) coeff(7) coeff(8) ];
            chisqfunc = @(x)(1.0/(size(yval,2)-numel(guess))*(yval(:)-(+x(1)*cosh((+exp(-x(2)))*(nt/2-xval(:)))+x(3)*(1-2*mod(xval(:),2)).*cosh((+exp(-x(4)))*(nt/2-xval(:)))))'*(ycorr\(yval(:)-(+x(1)*cosh((+exp(-x(2)))*(nt/2-xval(:)))+x(3)*(1-2*mod(xval(:),2)).*cosh((+exp(-x(4)))*(nt/2-xval(:)))))));
         case 6 % cosh 2 oscil 1
            guess = [coeff(1) coeff(2) coeff(3) coeff(4) coeff(7) coeff(8) ];
            chisqfunc = @(x)(1.0/(size(yval,2)-numel(guess))*(yval(:)-(+x(1)*cosh((+exp(-x(2)))*(nt/2-xval(:)))+x(3)*cosh((+exp(-x(2))+exp(-x(4)))*(nt/2-xval(:)))+x(5)*(1-2*mod(xval(:),2)).*cosh((+exp(-x(6)))*(nt/2-xval(:)))))'*(ycorr\(yval(:)-(+x(1)*cosh((+exp(-x(2)))*(nt/2-xval(:)))+x(3)*cosh((+exp(-x(2))+exp(-x(4)))*(nt/2-xval(:)))+x(5)*(1-2*mod(xval(:),2)).*cosh((+exp(-x(6)))*(nt/2-xval(:)))))));
         case 7 % cosh 3 oscil 1
            guess = [coeff(1) coeff(2) coeff(3) coeff(4) coeff(5) coeff(6) coeff(7) coeff(8) ];
            chisqfunc = @(x)(1.0/(size(yval,2)-numel(guess))*(yval(:)-(+x(1)*cosh((+exp(-x(2)))*(nt/2-xval(:)))+x(3)*cosh((+exp(-x(2))+exp(-x(4)))*(nt/2-xval(:)))+x(5)*cosh((+exp(-x(2))+exp(-x(4))+exp(-x(6)))*(nt/2-xval(:)))+x(7)*(1-2*mod(xval(:),2)).*cosh((+exp(-x(8)))*(nt/2-xval(:)))))'*(ycorr\(yval(:)-(+x(1)*cosh((+exp(-x(2)))*(nt/2-xval(:)))+x(3)*cosh((+exp(-x(2))+exp(-x(4)))*(nt/2-xval(:)))+x(5)*cosh((+exp(-x(2))+exp(-x(4))+exp(-x(6)))*(nt/2-xval(:)))+x(7)*(1-2*mod(xval(:),2)).*cosh((+exp(-x(8)))*(nt/2-xval(:)))))));
         case 8 % cosh 0 oscil 2
            guess = [coeff(7) coeff(8) coeff(9) coeff(10) ];
            chisqfunc = @(x)(1.0/(size(yval,2)-numel(guess))*(yval(:)-(+x(1)*(1-2*mod(xval(:),2)).*cosh((+exp(-x(2)))*(nt/2-xval(:)))+x(3)*(1-2*mod(xval(:),2)).*cosh((+exp(-x(2))+exp(-x(4)))*(nt/2-xval(:)))))'*(ycorr\(yval(:)-(+x(1)*(1-2*mod(xval(:),2)).*cosh((+exp(-x(2)))*(nt/2-xval(:)))+x(3)*(1-2*mod(xval(:),2)).*cosh((+exp(-x(2))+exp(-x(4)))*(nt/2-xval(:)))))));
         case 9 % cosh 1 oscil 2
            guess = [coeff(1) coeff(2) coeff(7) coeff(8) coeff(9) coeff(10) ];
            chisqfunc = @(x)(1.0/(size(yval,2)-numel(guess))*(yval(:)-(+x(1)*cosh((+exp(-x(2)))*(nt/2-xval(:)))+x(3)*(1-2*mod(xval(:),2)).*cosh((+exp(-x(4)))*(nt/2-xval(:)))+x(5)*(1-2*mod(xval(:),2)).*cosh((+exp(-x(4))+exp(-x(6)))*(nt/2-xval(:)))))'*(ycorr\(yval(:)-(+x(1)*cosh((+exp(-x(2)))*(nt/2-xval(:)))+x(3)*(1-2*mod(xval(:),2)).*cosh((+exp(-x(4)))*(nt/2-xval(:)))+x(5)*(1-2*mod(xval(:),2)).*cosh((+exp(-x(4))+exp(-x(6)))*(nt/2-xval(:)))))));
         case 10 % cosh 2 oscil 2
            guess = [coeff(1) coeff(2) coeff(3) coeff(4) coeff(7) coeff(8) coeff(9) coeff(10) ];
            chisqfunc = @(x)(1.0/(size(yval,2)-numel(guess))*(yval(:)-(+x(1)*cosh((+exp(-x(2)))*(nt/2-xval(:)))+x(3)*cosh((+exp(-x(2))+exp(-x(4)))*(nt/2-xval(:)))+x(5)*(1-2*mod(xval(:),2)).*cosh((+exp(-x(6)))*(nt/2-xval(:)))+x(7)*(1-2*mod(xval(:),2)).*cosh((+exp(-x(6))+exp(-x(8)))*(nt/2-xval(:)))))'*(ycorr\(yval(:)-(+x(1)*cosh((+exp(-x(2)))*(nt/2-xval(:)))+x(3)*cosh((+exp(-x(2))+exp(-x(4)))*(nt/2-xval(:)))+x(5)*(1-2*mod(xval(:),2)).*cosh((+exp(-x(6)))*(nt/2-xval(:)))+x(7)*(1-2*mod(xval(:),2)).*cosh((+exp(-x(6))+exp(-x(8)))*(nt/2-xval(:)))))));
         case 11 % cosh 3 oscil 2
            guess = [coeff(1) coeff(2) coeff(3) coeff(4) coeff(5) coeff(6) coeff(7) coeff(8) coeff(9) coeff(10) ];
            chisqfunc = @(x)(1.0/(size(yval,2)-numel(guess))*(yval(:)-(+x(1)*cosh((+exp(-x(2)))*(nt/2-xval(:)))+x(3)*cosh((+exp(-x(2))+exp(-x(4)))*(nt/2-xval(:)))+x(5)*cosh((+exp(-x(2))+exp(-x(4))+exp(-x(6)))*(nt/2-xval(:)))+x(7)*(1-2*mod(xval(:),2)).*cosh((+exp(-x(8)))*(nt/2-xval(:)))+x(9)*(1-2*mod(xval(:),2)).*cosh((+exp(-x(8))+exp(-x(10)))*(nt/2-xval(:)))))'*(ycorr\(yval(:)-(+x(1)*cosh((+exp(-x(2)))*(nt/2-xval(:)))+x(3)*cosh((+exp(-x(2))+exp(-x(4)))*(nt/2-xval(:)))+x(5)*cosh((+exp(-x(2))+exp(-x(4))+exp(-x(6)))*(nt/2-xval(:)))+x(7)*(1-2*mod(xval(:),2)).*cosh((+exp(-x(8)))*(nt/2-xval(:)))+x(9)*(1-2*mod(xval(:),2)).*cosh((+exp(-x(8))+exp(-x(10)))*(nt/2-xval(:)))))));
         case 12 % cosh 0 oscil 3
            guess = [coeff(7) coeff(8) coeff(9) coeff(10) coeff(11) coeff(12) ];
            chisqfunc = @(x)(1.0/(size(yval,2)-numel(guess))*(yval(:)-(+x(1)*(1-2*mod(xval(:),2)).*cosh((+exp(-x(2)))*(nt/2-xval(:)))+x(3)*(1-2*mod(xval(:),2)).*cosh((+exp(-x(2))+exp(-x(4)))*(nt/2-xval(:)))+x(5)*(1-2*mod(xval(:),2)).*cosh((+exp(-x(2))+exp(-x(4))+exp(-x(6)))*(nt/2-xval(:)))))'*(ycorr\(yval(:)-(+x(1)*(1-2*mod(xval(:),2)).*cosh((+exp(-x(2)))*(nt/2-xval(:)))+x(3)*(1-2*mod(xval(:),2)).*cosh((+exp(-x(2))+exp(-x(4)))*(nt/2-xval(:)))+x(5)*(1-2*mod(xval(:),2)).*cosh((+exp(-x(2))+exp(-x(4))+exp(-x(6)))*(nt/2-xval(:)))))));
         case 13 % cosh 1 oscil 3
            guess = [coeff(1) coeff(2) coeff(7) coeff(8) coeff(9) coeff(10) coeff(11) coeff(12) ];
            chisqfunc = @(x)(1.0/(size(yval,2)-numel(guess))*(yval(:)-(+x(1)*cosh((+exp(-x(2)))*(nt/2-xval(:)))+x(3)*(1-2*mod(xval(:),2)).*cosh((+exp(-x(4)))*(nt/2-xval(:)))+x(5)*(1-2*mod(xval(:),2)).*cosh((+exp(-x(4))+exp(-x(6)))*(nt/2-xval(:)))+x(7)*(1-2*mod(xval(:),2)).*cosh((+exp(-x(4))+exp(-x(6))+exp(-x(8)))*(nt/2-xval(:)))))'*(ycorr\(yval(:)-(+x(1)*cosh((+exp(-x(2)))*(nt/2-xval(:)))+x(3)*(1-2*mod(xval(:),2)).*cosh((+exp(-x(4)))*(nt/2-xval(:)))+x(5)*(1-2*mod(xval(:),2)).*cosh((+exp(-x(4))+exp(-x(6)))*(nt/2-xval(:)))+x(7)*(1-2*mod(xval(:),2)).*cosh((+exp(-x(4))+exp(-x(6))+exp(-x(8)))*(nt/2-xval(:)))))));
         case 14 % cosh 2 oscil 3
            guess = [coeff(1) coeff(2) coeff(3) coeff(4) coeff(7) coeff(8) coeff(9) coeff(10) coeff(11) coeff(12) ];
            chisqfunc = @(x)(1.0/(size(yval,2)-numel(guess))*(yval(:)-(+x(1)*cosh((+exp(-x(2)))*(nt/2-xval(:)))+x(3)*cosh((+exp(-x(2))+exp(-x(4)))*(nt/2-xval(:)))+x(5)*(1-2*mod(xval(:),2)).*cosh((+exp(-x(6)))*(nt/2-xval(:)))+x(7)*(1-2*mod(xval(:),2)).*cosh((+exp(-x(6))+exp(-x(8)))*(nt/2-xval(:)))+x(9)*(1-2*mod(xval(:),2)).*cosh((+exp(-x(6))+exp(-x(8))+exp(-x(10)))*(nt/2-xval(:)))))'*(ycorr\(yval(:)-(+x(1)*cosh((+exp(-x(2)))*(nt/2-xval(:)))+x(3)*cosh((+exp(-x(2))+exp(-x(4)))*(nt/2-xval(:)))+x(5)*(1-2*mod(xval(:),2)).*cosh((+exp(-x(6)))*(nt/2-xval(:)))+x(7)*(1-2*mod(xval(:),2)).*cosh((+exp(-x(6))+exp(-x(8)))*(nt/2-xval(:)))+x(9)*(1-2*mod(xval(:),2)).*cosh((+exp(-x(6))+exp(-x(8))+exp(-x(10)))*(nt/2-xval(:)))))));
         case 15 % cosh 3 oscil 3
            guess = [coeff(1) coeff(2) coeff(3) coeff(4) coeff(5) coeff(6) coeff(7) coeff(8) coeff(9) coeff(10) coeff(11) coeff(12) ];
            chisqfunc = @(x)(1.0/(size(yval,2)-numel(guess))*(yval(:)-(+x(1)*cosh((+exp(-x(2)))*(nt/2-xval(:)))+x(3)*cosh((+exp(-x(2))+exp(-x(4)))*(nt/2-xval(:)))+x(5)*cosh((+exp(-x(2))+exp(-x(4))+exp(-x(6)))*(nt/2-xval(:)))+x(7)*(1-2*mod(xval(:),2)).*cosh((+exp(-x(8)))*(nt/2-xval(:)))+x(9)*(1-2*mod(xval(:),2)).*cosh((+exp(-x(8))+exp(-x(10)))*(nt/2-xval(:)))+x(11)*(1-2*mod(xval(:),2)).*cosh((+exp(-x(8))+exp(-x(10))+exp(-x(12)))*(nt/2-xval(:)))))'*(ycorr\(yval(:)-(+x(1)*cosh((+exp(-x(2)))*(nt/2-xval(:)))+x(3)*cosh((+exp(-x(2))+exp(-x(4)))*(nt/2-xval(:)))+x(5)*cosh((+exp(-x(2))+exp(-x(4))+exp(-x(6)))*(nt/2-xval(:)))+x(7)*(1-2*mod(xval(:),2)).*cosh((+exp(-x(8)))*(nt/2-xval(:)))+x(9)*(1-2*mod(xval(:),2)).*cosh((+exp(-x(8))+exp(-x(10)))*(nt/2-xval(:)))+x(11)*(1-2*mod(xval(:),2)).*cosh((+exp(-x(8))+exp(-x(10))+exp(-x(12)))*(nt/2-xval(:)))))));
			otherwise % do case 1.
                wut = 'wut'
        end

		
		% check if we care about constraints.
		if (~exist('constraints', 'var'))
			% If there are no constraints, just fit it.
			if (numel(guess) > 1)
				[temp_guess, fit_output(count,15), succ_code] = fminsearch(chisqfunc, guess, optimset('TolX', fit_xprec,'TolFun', 1e-7,'MaxFunEvals',100000,'MaxIter',100000));
			else
				bound = max([abs(max(corr_fcn)) abs(min(corr_fcn))]);
				[temp_guess, fit_output(count,15)] = fminbnd(chisqfunc, -bound, bound, optimset('TolX', fit_xprec,'MaxFunEvals',100000,'MaxIter',100000));
			end
		else % There are constraints!
			% Any non-zero constraint becomes a gaussian.
			
			% Start with the non-mass constraints.
			for i=1:2:11
				if (constraints(i) ~= 0)
					constr{i} = @(c)((c-coeff(i)).^2/constraints(i));
				else
					constr{i} = @(c)(0);
				end
			end
			
			% Next, the mass constraints.
			
			% First cosh constraint.
			if (constraints(2) ~= 0)
				constr{2} = @(c)((c-coeff(2)).^2/(constraints(2)));
			else
				constr{2} = @(c)(0);
			end
			
			% The second is more tricky because we fit the differences.
			if (constraints(4) ~= 0)
				constr{4} = @(c,cm1)((exp(-c)+exp(-cm1)-exp(-coeff(2))-exp(-coeff(4))).^2/(constraints(4)));
			else
				constr{4} = @(c,cm1)(0);
			end
			
			% And the third!
			if (constraints(6) ~= 0)
				constr{6} = @(c,cm1,cm2)((exp(-c)+exp(-cm1)+exp(-cm2)-exp(-coeff(2))-exp(-coeff(4))-exp(-coeff(6))).^2/(constraints(6)));
			else
				constr{6} = @(c,cm1,cm2)(0);
			end
			
			% Repeat for oscil.
			if (constraints(8) ~= 0)
				constr{8} = @(c)((c-coeff(8)).^2/(constraints(8)));
			else
				constr{8} = @(c)(0);
			end
			
			% The second is more tricky because we fit the differences.
			if (constraints(10) ~= 0)
				constr{10} = @(c,cm1)((exp(-c)+exp(-cm1)-exp(-coeff(8))-exp(-coeff(10))).^2/(constraints(10)));
			else
				constr{10} = @(c,cm1)(0);
			end
			
			% And the third!
			if (constraints(12) ~= 0)
				constr{12} = @(c,cm1,cm2)((exp(-c)+exp(-cm1)+exp(-cm2)-exp(-coeff(8))-exp(-coeff(10))-exp(-coeff(12))).^2/(constraints(12)));
			else
				constr{12} = @(c,cm1,cm2)(0);
			end
			
			% Build up a superfunction!
			switch fitfunc
				case 1 % 1 0 
					chisqfunc_constr = @(x)(chisqfunc(x) + constr{1}(x(1)) + constr{2}(x(2)));
				case 2 % 2 0 
					chisqfunc_constr = @(x)(chisqfunc(x) + constr{1}(x(1)) + constr{2}(x(2)) + constr{3}(x(3)) + constr{4}(x(4),x(2)));
				case 3 % 3 0 
					chisqfunc_constr = @(x)(chisqfunc(x) + constr{1}(x(1)) + constr{2}(x(2)) + constr{3}(x(3)) + constr{4}(x(4),x(2)) + constr{5}(x(5)) + constr{6}(x(6),x(4),x(2)));
				case 4 % 0 1
					chisqfunc_constr = @(x)(chisqfunc(x) + constr{7}(x(1)) + constr{8}(x(2)));
				case 5 % 1 1
					chisqfunc_constr = @(x)(chisqfunc(x) + constr{1}(x(1)) + constr{2}(x(2)) + constr{7}(x(3)) + constr{8}(x(4)));
				case 6 % 2 1 
					chisqfunc_constr = @(x)(chisqfunc(x) + constr{1}(x(1)) + constr{2}(x(2)) + constr{3}(x(3)) + constr{4}(x(4),x(2)) + constr{7}(x(5)) + constr{8}(x(6)));
				case 7 % 3 1 
					chisqfunc_constr = @(x)(chisqfunc(x) + constr{1}(x(1)) + constr{2}(x(2)) + constr{3}(x(3)) + constr{4}(x(4),x(2)) + constr{5}(x(5)) + constr{6}(x(6),x(4),x(2)) + constr{7}(x(7)) + constr{8}(x(8)));	
				case 8 % 0 2
					chisqfunc_constr = @(x)(chisqfunc(x) + constr{7}(x(1)) + constr{8}(x(2)) + constr{9}(x(3)) + constr{10}(x(4),x(2)));
				case 9 % 1 2
					chisqfunc_constr = @(x)(chisqfunc(x) + constr{1}(x(1)) + constr{2}(x(2)) + constr{7}(x(3)) + constr{8}(x(4)) + constr{9}(x(5)) + constr{10}(x(6),x(4)));
				case 10 % 2 2 
					chisqfunc_constr = @(x)(chisqfunc(x) + constr{1}(x(1)) + constr{2}(x(2)) + constr{3}(x(3)) + constr{4}(x(4),x(2)) + constr{7}(x(5)) + constr{8}(x(6)) + constr{9}(x(7)) + constr{10}(x(8),x(6)));
				case 11 % 3 2 
					chisqfunc_constr = @(x)(chisqfunc(x) + constr{1}(x(1)) + constr{2}(x(2)) + constr{3}(x(3)) + constr{4}(x(4),x(2)) + constr{5}(x(5)) + constr{6}(x(6),x(4),x(2)) + constr{7}(x(7)) + constr{8}(x(8)) + constr{9}(x(9)) + constr{10}(x(10),x(8)));	
				case 12 % 0 3
					chisqfunc_constr = @(x)(chisqfunc(x) + constr{7}(x(1)) + constr{8}(x(2)) + constr{9}(x(3)) + constr{10}(x(4),x(2))+constr{11}(x(5))+constr{12}(x(6),x(4),x(2)));
				case 13 % 1 3
					chisqfunc_constr = @(x)(chisqfunc(x) + constr{1}(x(1)) + constr{2}(x(2)) + constr{7}(x(3)) + constr{8}(x(4)) + constr{9}(x(5)) + constr{10}(x(6),x(4))+constr{11}(x(7))+constr{12}(x(8),x(6),x(4)));
				case 14 % 2 3 
					chisqfunc_constr = @(x)(chisqfunc(x) + constr{1}(x(1)) + constr{2}(x(2)) + constr{3}(x(3)) + constr{4}(x(4),x(2)) + constr{7}(x(5)) + constr{8}(x(6)) + constr{9}(x(7)) + constr{10}(x(8),x(6))+constr{11}(x(9))+constr{12}(x(10),x(8),x(6)));
				case 15 % 3 3 
					chisqfunc_constr = @(x)(chisqfunc(x) + constr{1}(x(1)) + constr{2}(x(2)) + constr{3}(x(3)) + constr{4}(x(4),x(2)) + constr{5}(x(5)) + constr{6}(x(6),x(4),x(2)) + constr{7}(x(7)) + constr{8}(x(8)) + constr{9}(x(9)) + constr{10}(x(10),x(8))+constr{11}(x(11))+constr{12}(x(12),x(10),x(8)));	
			end
			
			if (numel(guess) > 1)
				[temp_guess, fit_output(count,15), succ_code] = fminsearch(chisqfunc_constr, guess, optimset('TolX', 1e-10,'TolFun', 1e-7,'MaxFunEvals',100000,'MaxIter',100000));
			else
				bound = max([abs(max(corr_fcn)) abs(min(corr_fcn))]);
				[temp_guess, fit_output(count,15)] = fminbnd(chisqfunc_constr, -bound, bound, optimset('TolX', 1e-12,'MaxFunEvals',100000,'MaxIter',100000));
			end			
			
		end

        if (fit_output(count, 15) > 1e-20 && succ_code == 1)

            fit_output(count,1) = t1;
            fit_output(count,2) = t2;
            
            
			% Convert back!
            if (fitcosh > 0)
                for a=1:fitcosh
                    fit_output(count, 2+2*a) = exp(-temp_guess(2*a));
                    if (fitcosh > 1 && a >1)
                        fit_output(count, 2+2*a) = fit_output(count, 2+2*a) + fit_output(count, 2+2*a-2);
                    end
                end
               fit_output(count, 3:2:(1+(fitcosh*2))) = temp_guess(1:2:(2*fitcosh-1));
            end
            
            if (fitoscil > 0)
                for a=1:fitoscil
                    fit_output(count, 8+2*a) = exp(-temp_guess(2*fitcosh+2*a));
                    if (fitoscil > 1 && a > 1)
                        fit_output(count, 8+2*a) = fit_output(count, 8+2*a) + fit_output(count, 8+2*a-2);
                    end
                end
               fit_output(count, 9:2:(7+(fitoscil*2))) = temp_guess((2*fitcosh+1):2:(2*(fitcosh+fitoscil)-1));
            end
            
            
            fit_output(count,16) = 1.0-chi2cdf(fit_output(count,15)*(size(yval,2)-numel(temp_guess)), size(yval,2)-numel(temp_guess));

            fit_output = fit_output(1,:);
        else
            fit_output = [];
        end
    
    else
        fit_output = [];
    end
    
    clear('choice');
    
end
