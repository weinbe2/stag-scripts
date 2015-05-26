function [fit_output] = get_all_nlfit(corr_fcn, corr_mat, tmin, tmax, nt, fiteven, fitfunc, fitdiag, coeff)
    % First, try the projected.
    fit_output = zeros(2, 13);

    count = 1;

    t1 = tmin;
	t2 = tmax;
	
	% One type of derivative fit takes half-int values. Prepare for that.
	if (fitfunc == 13) % Claudio's derivative
		if (fiteven == 1)
			error('fitfunc 13 (Claudios derivative) is incompatible with even fits for now.');
		end
		xval = (t1-0.5):1:(t2-0.5);
		yval = corr_fcn((ceil(t1)):(ceil(t2)))';
		ycorr = corr_mat((ceil(t1)):(ceil(t2)),(ceil(t1)):(ceil(t2)));
	elseif (fiteven == 0) % fit to all
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
	
	fit_output(count, 13) = cond(ycorr);
    
    
    choice = 'Yes';
    
    if (fit_output(count, 13) > 1e13)
        choice = 'No';
        %choice = questdlg(strcat(['Condition number = ' num2str(fit_output(count, 8)) '. Fit may crash MATLAB! Try anyway?']), 'Warning!');
    end
    
    if (strcmp(choice, 'Yes'))
	
        switch fitfunc
            case 1 % single cosh
                guess = [coeff(1) coeff(2)];
                chisqfunc = @(x)(1.0/(size(yval,2)-2)*(yval(:)-x(1)*cosh(x(2)*(nt/2-xval(:))))'*(ycorr\(yval(:)-x(1)*cosh(x(2)*(nt/2-xval(:))))));

            case 2 % two cosh
                guess = [coeff(1) coeff(2) coeff(3) coeff(4)];
                chisqfunc = @(x)(1.0/(size(yval,2)-4)*(yval(:)-x(1)*cosh(x(2)*(nt/2-xval(:)))-x(3)*cosh(x(4)*(nt/2-xval(:))))'*(ycorr\(yval(:)-x(1)*cosh(x(2)*(24-xval(:)))-x(3)*cosh(x(4)*(24-xval(:))))));

            case 3 % cosh + oscil
                guess = [coeff(1) coeff(2) coeff(5) coeff(6)];
                chisqfunc = @(x)(1.0/(size(yval,2)-4)*(yval(:)-x(1)*cosh(x(2)*(nt/2-xval(:)))-x(3)*(1-2*mod(xval(:),2)).*cosh(x(4)*(nt/2-xval(:))))'*(ycorr\(yval(:)-x(1)*cosh(x(2)*(nt/2-xval(:)))-x(3)*(1-2*mod(xval(:),2)).*cosh(x(4)*(nt/2-xval(:))))));

            case 4 % folded baryon
                guess = [coeff(1) coeff(2)]; 
                chisqfunc = @(x)(1.0/(size(yval,2)-2)*(yval(:)-x(1)*cosh(x(2)*(nt/2-xval(:)))+x(1)*(1-2*mod(xval(:),2)).*cosh(x(2)*(nt/2-xval(:))))'*(ycorr\(yval(:)-x(1)*cosh(x(2)*(nt/2-xval(:)))+x(1)*(1-2*mod(xval(:),2)).*cosh(x(2)*(nt/2-xval(:))))));
            
			case 5 % single cosh + const
                guess = [coeff(1) coeff(2) coeff(3)];
                chisqfunc = @(x)(1.0/(size(yval,2)-3)*(yval(:)-x(1)*cosh(x(2)*(nt/2-xval(:)))-x(3))'*(ycorr\(yval(:)-x(1)*cosh(x(2)*(nt/2-xval(:)))-x(3))));

			case 6 % single const - const, zero fit
				guess = [coeff(1) coeff(2)];
                chisqfunc = @(x)(1.0/(size(yval,2)-2)*(yval(:)-x(1)*cosh(x(2)*(nt/2-xval(:)))+x(1))'*(ycorr\(yval(:)-x(1)*cosh(x(2)*(nt/2-xval(:)))+x(1))));
			
            case 7 % cosh + osc + const!
                guess = [coeff(1) coeff(2) coeff(5) coeff(6) coeff(3)];
                chisqfunc = @(x)(1.0/(size(yval,2)-5)*(yval(:)-x(1)*cosh(x(2)*(nt/2-xval(:)))-x(3)*(1-2*mod(xval(:),2)).*cosh(x(4)*(nt/2-xval(:)))-x(5))'*(ycorr\(yval(:)-x(1)*cosh(x(2)*(nt/2-xval(:)))-x(3)*(1-2*mod(xval(:),2)).*cosh(x(4)*(nt/2-xval(:)))-x(5))));
            case 8 % cosh, fixed input mass.
				guess = [coeff(1)];
				chisqfunc = @(x)(1.0/(size(yval,2)-1)*(yval(:)-x(1)*cosh(coeff(2)*(nt/2-xval(:))))'*(ycorr\(yval(:)-x(1)*cosh(coeff(2)*(nt/2-xval(:))))));
			case 9 % cosh + osc, used for fixed mass fits.
				guess = [coeff(1) coeff(5)];
                chisqfunc = @(x)(1.0/(size(yval,2)-2)*(yval(:)-x(1)*cosh(coeff(2)*(nt/2-xval(:)))-x(2)*(1-2*mod(xval(:),2)).*cosh(coeff(6)*(nt/2-xval(:))))'*(ycorr\(yval(:)-x(1)*cosh(coeff(2)*(nt/2-xval(:)))-x(2)*(1-2*mod(xval(:),2)).*cosh(coeff(6)*(nt/2-xval(:))))));
			case 10 % discrete derivative by 2 of cosh.
				guess = [coeff(1) coeff(2)];
				chisqfunc = @(x)(1.0/(size(yval,2)-2)*(yval(:)-x(1)*(((cosh(2*x(2))-1)/2.0)*cosh(x(2)*(nt/2-xval(:))) - (sinh(2*x(2))/2.0)*sinh(x(2)*(nt/2-xval(:)))))'*(ycorr\(yval(:)-x(1)*(((cosh(2*x(2))-1)/2.0)*cosh(x(2)*(nt/2-xval(:))) - (sinh(2*x(2))/2.0)*sinh(x(2)*(nt/2-xval(:)))))));
			case 11 % discrete derivative by 2 of cosh+oscil.
				guess = [coeff(1) coeff(2) coeff(5) coeff(6)];
				chisqfunc = @(x)(1.0/(size(yval,2)-4)*(yval(:)-x(1)*(((cosh(2*x(2))-1)/2.0)*cosh(x(2)*(nt/2-xval(:))) - (sinh(2*x(2))/2.0)*sinh(x(2)*(nt/2-xval(:))))-x(3)*(1-2*mod(xval(:),2)).*(((cosh(2*x(4))-1)/2.0)*cosh(x(4)*(nt/2-xval(:))) - (sinh(2*x(4))/2.0)*sinh(x(4)*(nt/2-xval(:)))))'*(ycorr\(yval(:)-x(1)*(((cosh(2*x(2))-1)/2.0)*cosh(x(2)*(nt/2-xval(:))) - (sinh(2*x(2))/2.0)*sinh(x(2)*(nt/2-xval(:))))-x(3)*(1-2*mod(xval(:),2)).*(((cosh(2*x(4))-1)/2.0)*cosh(x(4)*(nt/2-xval(:))) - (sinh(2*x(4))/2.0)*sinh(x(4)*(nt/2-xval(:)))))));
			case 12 % unfolded baryon
				guess = [coeff(1) coeff(2) coeff(5) coeff(6)];
				chisqfunc = @(x)(1.0/(size(yval,2)-4)*(yval(:) - x(1)*(exp(-x(2)*xval(:)) - (1-2*mod(xval(:),2)).*exp(-x(2)*(nt-xval(:)))) + x(3).*(1-2*mod(xval(:),2)).*(exp(-x(4)*xval(:)) - (1-2*mod(xval(:),2)).*exp(-x(4)*(nt-xval(:)))))'*(ycorr\(yval(:) - x(1)*(exp(-x(2)*xval(:)) - (1-2*mod(xval(:),2)).*exp(-x(2)*(nt-xval(:)))) + x(3).*(1-2*mod(xval(:),2)).*(exp(-x(4)*xval(:)) - (1-2*mod(xval(:),2)).*exp(-x(4)*(nt-xval(:)))))));
			case 13 % Claudio derivative
                guess = [coeff(1) coeff(2) coeff(5) coeff(6)];
                chisqfunc = @(x)(1.0/(size(yval,2)-4)*(yval(:)+x(1).*sinh(x(2)*(nt/2-xval(:)))+x(3)*(1-2*mod(round(xval(:)+0.5),2)).*cosh(x(4).*(nt/2-xval(:))))'*(ycorr\(yval(:)+x(1).*sinh(x(2)*(nt/2-xval(:)))+x(3).*(1-2*mod(round(xval(:)+0.5),2)).*cosh(x(4)*(nt/2-xval(:))))));
			otherwise % do case 1.
                wut = 'wut'
        end

		if (numel(guess) > 1)
			[temp_guess, fit_output(count,11), succ_code] = fminsearch(chisqfunc, guess, optimset('TolX', 1e-15,'TolFun', 1e-7,'MaxFunEvals',100000,'MaxIter',100000));
		else
			bound = max([abs(max(corr_fcn)) abs(min(corr_fcn))]);
			[temp_guess, fit_output(count,11)] = fminbnd(chisqfunc, -bound, bound, optimset('TolX', 1e-12,'MaxFunEvals',100000,'MaxIter',100000));
		end

        if (fit_output(count, 11) > 1e-20 && succ_code == 1)

            fit_output(count,1) = t1;
            fit_output(count,2) = t2;
            
            switch fitfunc
                case 1 % single cosh
                    fit_output(count, 3) = temp_guess(1);
                    fit_output(count, 4) = abs(temp_guess(2));
                    fit_output(count, 5) = 0.0;
                    fit_output(count, 6) = 0.0;
                    fit_output(count, 7) = 0.0;
                    fit_output(count, 8) = 0.0;
                    fit_output(count, 9) = 0.0;
                    fit_output(count, 10) = 0.0;
                    
                case 2 % two cosh
                    fit_output(count, 3) = temp_guess(1);
                    fit_output(count, 4) = abs(temp_guess(2));
                    fit_output(count, 5) = temp_guess(3);
                    fit_output(count, 6) = abs(temp_guess(4));
                    fit_output(count, 7) = 0.0;
                    fit_output(count, 8) = 0.0;
                    fit_output(count, 9) = 0.0;
                    fit_output(count, 10) = 0.0;
                    
                case 3 % cosh + oscil
                    fit_output(count, 3) = temp_guess(1);
                    fit_output(count, 4) = abs(temp_guess(2));
                    fit_output(count, 5) = 0.0;
                    fit_output(count, 6) = 0.0;
                    fit_output(count, 7) = temp_guess(3);
                    fit_output(count, 8) = abs(temp_guess(4));
                    fit_output(count, 9) = 0.0;
                    fit_output(count, 10) = 0.0;
                case 4 % folded baryon
                    fit_output(count, 3) = temp_guess(1);
                    fit_output(count, 4) = abs(temp_guess(2));
                    fit_output(count, 5) = 0.0;
                    fit_output(count, 6) = 0.0;
                    fit_output(count, 7) = 0.0;
                    fit_output(count, 8) = 0.0;
                    fit_output(count, 9) = 0.0;
                    fit_output(count, 10) = 0.0;
                    
                case 5 % single cosh + const
                    fit_output(count, 3) = temp_guess(1);
                    fit_output(count, 4) = abs(temp_guess(2));
                    fit_output(count, 5) = temp_guess(3);
                    fit_output(count, 6) = 0.0;
                    fit_output(count, 7) = 0.0;
                    fit_output(count, 8) = 0.0;
                    fit_output(count, 9) = 0.0;
                    fit_output(count, 10) = 0.0;
                    
                case 6 % single const - const, zero fit
                    fit_output(count, 3) = temp_guess(1);
                    fit_output(count, 4) = abs(temp_guess(2));
                    fit_output(count, 5) = 0.0;
                    fit_output(count, 6) = 0.0;
                    fit_output(count, 7) = 0.0;
                    fit_output(count, 8) = 0.0;
                    fit_output(count, 9) = 0.0;
                    fit_output(count, 10) = 0.0;
                    
                case 7 % cosh + osc + const!
                    fit_output(count, 3) = temp_guess(1);
                    fit_output(count, 4) = abs(temp_guess(2));
                    fit_output(count, 5) = temp_guess(5);
                    fit_output(count, 6) = 0.0;
                    fit_output(count, 7) = temp_guess(3);
                    fit_output(count, 8) = abs(temp_guess(4));
                    fit_output(count, 9) = 0.0;
                    fit_output(count, 10) = 0.0;
					
				case 8 % cosh with fixed mass
					fit_output(count, 3) = temp_guess(1);
                    fit_output(count, 4) = coeff(2);
                    fit_output(count, 5) = 0.0;
                    fit_output(count, 6) = 0.0;
                    fit_output(count, 7) = 0.0;
                    fit_output(count, 8) = 0.0;
                    fit_output(count, 9) = 0.0;
                    fit_output(count, 10) = 0.0;
					
				case 9 % cosh + osc with fixed masses
					fit_output(count, 3) = temp_guess(1);
                    fit_output(count, 4) = coeff(2);
                    fit_output(count, 5) = 0.0;
                    fit_output(count, 6) = 0.0;
                    fit_output(count, 7) = temp_guess(2);
                    fit_output(count, 8) = coeff(6);
                    fit_output(count, 9) = 0.0;
                    fit_output(count, 10) = 0.0;
				
				case 10 % discrete by 2 derivative of cosh
					fit_output(count, 3) = temp_guess(1);
                    fit_output(count, 4) = abs(temp_guess(2));
                    fit_output(count, 5) = 0.0;
                    fit_output(count, 6) = 0.0;
                    fit_output(count, 7) = 0.0;
                    fit_output(count, 8) = 0.0;
                    fit_output(count, 9) = 0.0;
                    fit_output(count, 10) = 0.0;
					
				case 11 % discrete by 2 derivative of cosh+oscil
                    fit_output(count, 3) = temp_guess(1);
                    fit_output(count, 4) = abs(temp_guess(2));
                    fit_output(count, 5) = 0.0;
                    fit_output(count, 6) = 0.0;
                    fit_output(count, 7) = temp_guess(3);
                    fit_output(count, 8) = abs(temp_guess(4));
                    fit_output(count, 9) = 0.0;
                    fit_output(count, 10) = 0.0;
					
				case 12 % baryon
                    fit_output(count, 3) = temp_guess(1);
                    fit_output(count, 4) = abs(temp_guess(2));
                    fit_output(count, 5) = 0.0;
                    fit_output(count, 6) = 0.0;
                    fit_output(count, 7) = temp_guess(3);
                    fit_output(count, 8) = abs(temp_guess(4));
                    fit_output(count, 9) = 0.0;
                    fit_output(count, 10) = 0.0;
					
				case 13 % derivative fit (Claudio)
                    fit_output(count, 3) = temp_guess(1);
                    fit_output(count, 4) = abs(temp_guess(2));
                    fit_output(count, 5) = 0.0;
                    fit_output(count, 6) = 0.0;
                    fit_output(count, 7) = temp_guess(3);
                    fit_output(count, 8) = abs(temp_guess(4));
                    fit_output(count, 9) = 0.0;
                    fit_output(count, 10) = 0.0;
                    
                otherwise % do case 1.
                    wut = 'wut'
            end

            
            fit_output(count,12) = 1.0-chi2cdf(fit_output(count,11)*(size(yval,2)-numel(temp_guess)), size(yval,2)-numel(temp_guess));

            fit_output = fit_output(1,:);
        else
            fit_output = [];
        end
    
    else
        fit_output = [];
    end
    
    clear('choice');
    
end
