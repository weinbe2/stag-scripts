function [fit_output] = get_all_nlfit_multi(corr_fcn, corr_mat, tmin, tmax, nt, fiteven, fitcut, funcdiff, fitcosh, fitoscil, fitdiag, fit_xprec, coeff, constraints)

	% add the path for jacobianest
	%addpath('.\multifit\support','-end');

	% Get the global dircosh, modcosh functions.
	global dircosh;
	global modcosh;

    % Fit some variable number of cosh or oscillating terms, up to 3 of each!
	fitfunc = fitcosh+4*fitoscil;
	
    fit_output = zeros(2, 17);

    count = 1;

    t1 = tmin;
	t2 = tmax;
	
	if (fiteven == 0) % fit to all
		if (funcdiff == 0) % still all
			xval = t1:t2;
			yval = corr_fcn((t1+1):(t2+1))';
			ycorr = corr_mat((t1+1):(t2+1),(t1+1):(t2+1));
		else % finite diff! one less data.
			xval = t1:(t2-1);
			yval = corr_fcn((t1+1):t2)';
			ycorr = corr_mat((t1+1):t2,(t1+1):t2);
		end
	else % only fit to even...
		if (funcdiff == 0) % no finite diff.
			if (mod(t1,2)==0)
				xval = t1:2:t2;
			else
				xval = (t1+1):2:t2;
			end
			yval = corr_fcn((xval+1))';
			ycorr = corr_mat((xval+1),(xval+1));
		else % finite diff! careful about end.
			if (mod(t1,2)==0)
				xval = t1:2:(t2-1);
			else
				xval = (t1+1):2:(t2-1);
			end
			yval = corr_fcn((xval+1))';
			ycorr = corr_mat((xval+1),(xval+1));
		end
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
		% We get two things here: fitparamfunc, which has the log
		% ordering map for the masses, and fitmassfunc, which doesn't
		% have the log map. The former is used for fits, latter for
		% estimating classical errors. 
        switch fitfunc
         case 1 % cosh 1 oscil 0
            guess = [coeff(1) coeff(2) ];
			fitparamfunc = @(x,xv)(+x(1)*dircosh(+exp(-x(2)),xv,nt));
			fitmassfunc = @(x,xv)(+x(1)*dircosh(x(2),xv,nt));
            %chisqfunc = @(x)(1.0/(size(yval,2)-numel(guess))*(yval(:)-(+x(1)*cosh((+exp(-x(2)))*(nt/2-xval(:)))))'*(ycorr\(yval(:)-(+x(1)*cosh((+exp(-x(2)))*(nt/2-xval(:)))))));
         case 2 % cosh 2 oscil 0
            guess = [coeff(1) coeff(2) coeff(3) coeff(4) ];
			fitparamfunc = @(x,xv)(+x(1)*dircosh(+exp(-x(2)),xv,nt)+x(3)*dircosh(+exp(-x(2))+exp(-x(4)),xv,nt));
			fitmassfunc = @(x,xv)(+x(1)*dircosh(x(2),xv,nt)+x(3)*dircosh(x(4),xv,nt));
            %chisqfunc = @(x)(1.0/(size(yval,2)-numel(guess))*(yval(:)-(+x(1)*cosh((+exp(-x(2)))*(nt/2-xval(:)))+x(3)*cosh((+exp(-x(2))+exp(-x(4)))*(nt/2-xval(:)))))'*(ycorr\(yval(:)-(+x(1)*cosh((+exp(-x(2)))*(nt/2-xval(:)))+x(3)*cosh((+exp(-x(2))+exp(-x(4)))*(nt/2-xval(:)))))));
         case 3 % cosh 3 oscil 0
            guess = [coeff(1) coeff(2) coeff(3) coeff(4) coeff(5) coeff(6) ];
			fitparamfunc = @(x,xv)(+x(1)*dircosh(+exp(-x(2)),xv,nt)+x(3)*dircosh(+exp(-x(2))+exp(-x(4)),xv,nt)+x(5)*dircosh(+exp(-x(2))+exp(-x(4))+exp(-x(6)),xv,nt));
			fitmassfunc = @(x,xv)(+x(1)*dircosh(x(2),xv,nt)+x(3)*dircosh(x(4),xv,nt)+x(5)*dircosh(x(6),xv,nt));
            %chisqfunc = @(x)(1.0/(size(yval,2)-numel(guess))*(yval(:)-(+x(1)*cosh((+exp(-x(2)))*(nt/2-xval(:)))+x(3)*cosh((+exp(-x(2))+exp(-x(4)))*(nt/2-xval(:)))+x(5)*cosh((+exp(-x(2))+exp(-x(4))+exp(-x(6)))*(nt/2-xval(:)))))'*(ycorr\(yval(:)-(+x(1)*cosh((+exp(-x(2)))*(nt/2-xval(:)))+x(3)*cosh((+exp(-x(2))+exp(-x(4)))*(nt/2-xval(:)))+x(5)*cosh((+exp(-x(2))+exp(-x(4))+exp(-x(6)))*(nt/2-xval(:)))))));
         case 4 % cosh 0 oscil 1
            guess = [coeff(7) coeff(8) ];
			fitparamfunc = @(x,xv)(+x(1)*modcosh(+exp(-x(2)), xv, nt));
			fitmassfunc = @(x,xv)(+x(1)*modcosh(x(2), xv, nt));
            %chisqfunc = @(x)(1.0/(size(yval,2)-numel(guess))*(yval(:)-(+x(1)*(1-2*mod(xval(:),2)).*cosh((+exp(-x(2)))*(nt/2-xval(:)))))'*(ycorr\(yval(:)-(+x(1)*(1-2*mod(xval(:),2)).*cosh((+exp(-x(2)))*(nt/2-xval(:)))))));
         case 5 % cosh 1 oscil 1
            guess = [coeff(1) coeff(2) coeff(7) coeff(8) ];
			fitparamfunc = @(x,xv)(+x(1)*dircosh(+exp(-x(2)),xv,nt)+x(3)*modcosh(+exp(-x(4)),xv,nt));
			fitmassfunc = @(x,xv)(+x(1)*dircosh(x(2),xv,nt)+x(3)*modcosh(x(4),xv,nt));
            %chisqfunc = @(x)(1.0/(size(yval,2)-numel(guess))*(yval(:)-(+x(1)*cosh((+exp(-x(2)))*(nt/2-xval(:)))+x(3)*(1-2*mod(xval(:),2)).*cosh((+exp(-x(4)))*(nt/2-xval(:)))))'*(ycorr\(yval(:)-(+x(1)*cosh((+exp(-x(2)))*(nt/2-xval(:)))+x(3)*(1-2*mod(xval(:),2)).*cosh((+exp(-x(4)))*(nt/2-xval(:)))))));
         case 6 % cosh 2 oscil 1
            guess = [coeff(1) coeff(2) coeff(3) coeff(4) coeff(7) coeff(8) ];
			fitparamfunc = @(x,xv)(+x(1)*dircosh(+exp(-x(2)),xv,nt)+x(3)*dircosh(+exp(-x(2))+exp(-x(4)),xv,nt)+x(5)*modcosh(+exp(-x(6)),xv,nt));
			fitmassfunc = @(x,xv)(+x(1)*dircosh(x(2),xv,nt)+x(3)*dircosh(x(4),xv,nt)+x(5)*modcosh(x(6),xv,nt));
            %chisqfunc = @(x)(1.0/(size(yval,2)-numel(guess))*(yval(:)-(+x(1)*cosh((+exp(-x(2)))*(nt/2-xval(:)))+x(3)*cosh((+exp(-x(2))+exp(-x(4)))*(nt/2-xval(:)))+x(5)*(1-2*mod(xval(:),2)).*cosh((+exp(-x(6)))*(nt/2-xval(:)))))'*(ycorr\(yval(:)-(+x(1)*cosh((+exp(-x(2)))*(nt/2-xval(:)))+x(3)*cosh((+exp(-x(2))+exp(-x(4)))*(nt/2-xval(:)))+x(5)*(1-2*mod(xval(:),2)).*cosh((+exp(-x(6)))*(nt/2-xval(:)))))));
         case 7 % cosh 3 oscil 1
            guess = [coeff(1) coeff(2) coeff(3) coeff(4) coeff(5) coeff(6) coeff(7) coeff(8) ];
			fitparamfunc = @(x,xv)(+x(1)*dircosh(+exp(-x(2)),xv,nt)+x(3)*dircosh(+exp(-x(2))+exp(-x(4)),xv,nt)+x(5)*dircosh(+exp(-x(2))+exp(-x(4))+exp(-x(6)),xv,nt)+x(7)*modcosh(+exp(-x(8)),xv,nt));
			fitmassfunc = @(x,xv)(+x(1)*dircosh(x(2),xv,nt)+x(3)*dircosh(x(4),xv,nt)+x(5)*dircosh(x(6),xv,nt)+x(7)*modcosh(x(8),xv,nt));
            %chisqfunc = @(x)(1.0/(size(yval,2)-numel(guess))*(yval(:)-(+x(1)*cosh((+exp(-x(2)))*(nt/2-xval(:)))+x(3)*cosh((+exp(-x(2))+exp(-x(4)))*(nt/2-xval(:)))+x(5)*cosh((+exp(-x(2))+exp(-x(4))+exp(-x(6)))*(nt/2-xval(:)))+x(7)*(1-2*mod(xval(:),2)).*cosh((+exp(-x(8)))*(nt/2-xval(:)))))'*(ycorr\(yval(:)-(+x(1)*cosh((+exp(-x(2)))*(nt/2-xval(:)))+x(3)*cosh((+exp(-x(2))+exp(-x(4)))*(nt/2-xval(:)))+x(5)*cosh((+exp(-x(2))+exp(-x(4))+exp(-x(6)))*(nt/2-xval(:)))+x(7)*(1-2*mod(xval(:),2)).*cosh((+exp(-x(8)))*(nt/2-xval(:)))))));
         case 8 % cosh 0 oscil 2
            guess = [coeff(7) coeff(8) coeff(9) coeff(10) ];
			fitparamfunc = @(x,xv)(+x(1)*modcosh(+exp(-x(2)),xv,nt)+x(3)*modcosh(+exp(-x(2))+exp(-x(4)),xv,nt));
			fitparamfunc = @(x,xv)(+x(1)*modcosh(x(2),xv,nt)+x(3)*modcosh(x(4),xv,nt));
            %chisqfunc = @(x)(1.0/(size(yval,2)-numel(guess))*(yval(:)-(+x(1)*(1-2*mod(xval(:),2)).*cosh((+exp(-x(2)))*(nt/2-xval(:)))+x(3)*(1-2*mod(xval(:),2)).*cosh((+exp(-x(2))+exp(-x(4)))*(nt/2-xval(:)))))'*(ycorr\(yval(:)-(+x(1)*(1-2*mod(xval(:),2)).*cosh((+exp(-x(2)))*(nt/2-xval(:)))+x(3)*(1-2*mod(xval(:),2)).*cosh((+exp(-x(2))+exp(-x(4)))*(nt/2-xval(:)))))));
         case 9 % cosh 1 oscil 2
            guess = [coeff(1) coeff(2) coeff(7) coeff(8) coeff(9) coeff(10) ];
			fitparamfunc = @(x,xv)(+x(1)*dircosh(+exp(-x(2)),xv,nt)+x(3)*modcosh(+exp(-x(4)),xv,nt)+x(5)*modcosh(+exp(-x(4))+exp(-x(6)),xv,nt));
			fitmassfunc = @(x,xv)(+x(1)*dircosh(x(2),xv,nt)+x(3)*modcosh(x(4),xv,nt)+x(5)*modcosh(x(6),xv,nt));
            %chisqfunc = @(x)(1.0/(size(yval,2)-numel(guess))*(yval(:)-(+x(1)*cosh((+exp(-x(2)))*(nt/2-xval(:)))+x(3)*(1-2*mod(xval(:),2)).*cosh((+exp(-x(4)))*(nt/2-xval(:)))+x(5)*(1-2*mod(xval(:),2)).*cosh((+exp(-x(4))+exp(-x(6)))*(nt/2-xval(:)))))'*(ycorr\(yval(:)-(+x(1)*cosh((+exp(-x(2)))*(nt/2-xval(:)))+x(3)*(1-2*mod(xval(:),2)).*cosh((+exp(-x(4)))*(nt/2-xval(:)))+x(5)*(1-2*mod(xval(:),2)).*cosh((+exp(-x(4))+exp(-x(6)))*(nt/2-xval(:)))))));
         case 10 % cosh 2 oscil 2
            guess = [coeff(1) coeff(2) coeff(3) coeff(4) coeff(7) coeff(8) coeff(9) coeff(10) ];
			fitparamfunc = @(x,xv)(+x(1)*dircosh(+exp(-x(2)),xv,nt)+x(3)*dircosh(+exp(-x(2))+exp(-x(4)),xv,nt)+x(5)*modcosh(+exp(-x(6)),xv,nt)+x(7)*modcosh(+exp(-x(6))+exp(-x(8)),xv,nt));
			fitmassfunc = @(x,xv)(+x(1)*dircosh(x(2),xv,nt)+x(3)*dircosh(x(4),xv,nt)+x(5)*modcosh(x(6),xv,nt)+x(7)*modcosh(x(8),xv,nt));
            %chisqfunc = @(x)(1.0/(size(yval,2)-numel(guess))*(yval(:)-(+x(1)*cosh((+exp(-x(2)))*(nt/2-xval(:)))+x(3)*cosh((+exp(-x(2))+exp(-x(4)))*(nt/2-xval(:)))+x(5)*(1-2*mod(xval(:),2)).*cosh((+exp(-x(6)))*(nt/2-xval(:)))+x(7)*(1-2*mod(xval(:),2)).*cosh((+exp(-x(6))+exp(-x(8)))*(nt/2-xval(:)))))'*(ycorr\(yval(:)-(+x(1)*cosh((+exp(-x(2)))*(nt/2-xval(:)))+x(3)*cosh((+exp(-x(2))+exp(-x(4)))*(nt/2-xval(:)))+x(5)*(1-2*mod(xval(:),2)).*cosh((+exp(-x(6)))*(nt/2-xval(:)))+x(7)*(1-2*mod(xval(:),2)).*cosh((+exp(-x(6))+exp(-x(8)))*(nt/2-xval(:)))))));
         case 11 % cosh 3 oscil 2
            guess = [coeff(1) coeff(2) coeff(3) coeff(4) coeff(5) coeff(6) coeff(7) coeff(8) coeff(9) coeff(10) ];
			fitparamfunc = @(x,xv)(+x(1)*dircosh(+exp(-x(2)),xv,nt)+x(3)*dircosh(+exp(-x(2))+exp(-x(4)),xv,nt)+x(5)*dircosh(+exp(-x(2))+exp(-x(4))+exp(-x(6)),xv,nt)+x(7)*modcosh(+exp(-x(8)),xv,nt)+x(9)*modcosh(+exp(-x(8))+exp(-x(10)),xv,nt));
			fitmassfunc = @(x,xv)(+x(1)*dircosh(x(2),xv,nt)+x(3)*dircosh(x(4),xv,nt)+x(5)*dircosh(x(6),xv,nt)+x(7)*modcosh(x(8),xv,nt)+x(9)*modcosh(x(10),xv,nt));
            %chisqfunc = @(x)(1.0/(size(yval,2)-numel(guess))*(yval(:)-(+x(1)*cosh((+exp(-x(2)))*(nt/2-xval(:)))+x(3)*cosh((+exp(-x(2))+exp(-x(4)))*(nt/2-xval(:)))+x(5)*cosh((+exp(-x(2))+exp(-x(4))+exp(-x(6)))*(nt/2-xval(:)))+x(7)*(1-2*mod(xval(:),2)).*cosh((+exp(-x(8)))*(nt/2-xval(:)))+x(9)*(1-2*mod(xval(:),2)).*cosh((+exp(-x(8))+exp(-x(10)))*(nt/2-xval(:)))))'*(ycorr\(yval(:)-(+x(1)*cosh((+exp(-x(2)))*(nt/2-xval(:)))+x(3)*cosh((+exp(-x(2))+exp(-x(4)))*(nt/2-xval(:)))+x(5)*cosh((+exp(-x(2))+exp(-x(4))+exp(-x(6)))*(nt/2-xval(:)))+x(7)*(1-2*mod(xval(:),2)).*cosh((+exp(-x(8)))*(nt/2-xval(:)))+x(9)*(1-2*mod(xval(:),2)).*cosh((+exp(-x(8))+exp(-x(10)))*(nt/2-xval(:)))))));
         case 12 % cosh 0 oscil 3
            guess = [coeff(7) coeff(8) coeff(9) coeff(10) coeff(11) coeff(12) ];
			fitparamfunc = @(x,xv)(+x(1)*modcosh(+exp(-x(2)),xv,nt)+x(3)*modcosh(+exp(-x(2))+exp(-x(4)),xv,nt)+x(5)*modcosh(+exp(-x(2))+exp(-x(4))+exp(-x(6)),xv,nt));
			fitmassfunc = @(x,xv)(+x(1)*modcosh(x(2),xv,nt)+x(3)*modcosh(x(4),xv,nt)+x(5)*modcosh(x(6),xv,nt));
            %chisqfunc = @(x)(1.0/(size(yval,2)-numel(guess))*(yval(:)-(+x(1)*(1-2*mod(xval(:),2)).*cosh((+exp(-x(2)))*(nt/2-xval(:)))+x(3)*(1-2*mod(xval(:),2)).*cosh((+exp(-x(2))+exp(-x(4)))*(nt/2-xval(:)))+x(5)*(1-2*mod(xval(:),2)).*cosh((+exp(-x(2))+exp(-x(4))+exp(-x(6)))*(nt/2-xval(:)))))'*(ycorr\(yval(:)-(+x(1)*(1-2*mod(xval(:),2)).*cosh((+exp(-x(2)))*(nt/2-xval(:)))+x(3)*(1-2*mod(xval(:),2)).*cosh((+exp(-x(2))+exp(-x(4)))*(nt/2-xval(:)))+x(5)*(1-2*mod(xval(:),2)).*cosh((+exp(-x(2))+exp(-x(4))+exp(-x(6)))*(nt/2-xval(:)))))));
         case 13 % cosh 1 oscil 3
            guess = [coeff(1) coeff(2) coeff(7) coeff(8) coeff(9) coeff(10) coeff(11) coeff(12) ];
			fitparamfunc = @(x,xv)(+x(1)*dircosh(+exp(-x(2)),xv,nt)+x(3)*modcosh(+exp(-x(4)),xv,nt)+x(5)*modcosh(+exp(-x(4))+exp(-x(6)),xv,nt)+x(7)*modcosh(+exp(-x(4))+exp(-x(6))+exp(-x(8)),xv,nt));
			fitmassfunc = @(x,xv)(+x(1)*dircosh(x(2),xv,nt)+x(3)*modcosh(x(4),xv,nt)+x(5)*modcosh(x(6),xv,nt)+x(7)*modcosh(x(8),xv,nt));
            %chisqfunc = @(x)(1.0/(size(yval,2)-numel(guess))*(yval(:)-(+x(1)*cosh((+exp(-x(2)))*(nt/2-xval(:)))+x(3)*(1-2*mod(xval(:),2)).*cosh((+exp(-x(4)))*(nt/2-xval(:)))+x(5)*(1-2*mod(xval(:),2)).*cosh((+exp(-x(4))+exp(-x(6)))*(nt/2-xval(:)))+x(7)*(1-2*mod(xval(:),2)).*cosh((+exp(-x(4))+exp(-x(6))+exp(-x(8)))*(nt/2-xval(:)))))'*(ycorr\(yval(:)-(+x(1)*cosh((+exp(-x(2)))*(nt/2-xval(:)))+x(3)*(1-2*mod(xval(:),2)).*cosh((+exp(-x(4)))*(nt/2-xval(:)))+x(5)*(1-2*mod(xval(:),2)).*cosh((+exp(-x(4))+exp(-x(6)))*(nt/2-xval(:)))+x(7)*(1-2*mod(xval(:),2)).*cosh((+exp(-x(4))+exp(-x(6))+exp(-x(8)))*(nt/2-xval(:)))))));
         case 14 % cosh 2 oscil 3
            guess = [coeff(1) coeff(2) coeff(3) coeff(4) coeff(7) coeff(8) coeff(9) coeff(10) coeff(11) coeff(12) ];
			fitparamfunc = @(x,xv)(+x(1)*dircosh(+exp(-x(2)),xv,nt)+x(3)*dircosh(+exp(-x(2))+exp(-x(4)),xv,nt)+x(5)*modcosh(+exp(-x(6)),xv,nt)+x(7)*modcosh(+exp(-x(6))+exp(-x(8)),xv,nt)+x(9)*modcosh(+exp(-x(6))+exp(-x(8))+exp(-x(10)),xv,nt));
			fitmassfunc = @(x,xv)(+x(1)*dircosh(x(2),xv,nt)+x(3)*dircosh(x(4),xv,nt)+x(5)*modcosh(x(6),xv,nt)+x(7)*modcosh(x(8),xv,nt)+x(9)*modcosh(x(10),xv,nt));
            %chisqfunc = @(x)(1.0/(size(yval,2)-numel(guess))*(yval(:)-(+x(1)*cosh((+exp(-x(2)))*(nt/2-xval(:)))+x(3)*cosh((+exp(-x(2))+exp(-x(4)))*(nt/2-xval(:)))+x(5)*(1-2*mod(xval(:),2)).*cosh((+exp(-x(6)))*(nt/2-xval(:)))+x(7)*(1-2*mod(xval(:),2)).*cosh((+exp(-x(6))+exp(-x(8)))*(nt/2-xval(:)))+x(9)*(1-2*mod(xval(:),2)).*cosh((+exp(-x(6))+exp(-x(8))+exp(-x(10)))*(nt/2-xval(:)))))'*(ycorr\(yval(:)-(+x(1)*cosh((+exp(-x(2)))*(nt/2-xval(:)))+x(3)*cosh((+exp(-x(2))+exp(-x(4)))*(nt/2-xval(:)))+x(5)*(1-2*mod(xval(:),2)).*cosh((+exp(-x(6)))*(nt/2-xval(:)))+x(7)*(1-2*mod(xval(:),2)).*cosh((+exp(-x(6))+exp(-x(8)))*(nt/2-xval(:)))+x(9)*(1-2*mod(xval(:),2)).*cosh((+exp(-x(6))+exp(-x(8))+exp(-x(10)))*(nt/2-xval(:)))))));
         case 15 % cosh 3 oscil 3
            guess = [coeff(1) coeff(2) coeff(3) coeff(4) coeff(5) coeff(6) coeff(7) coeff(8) coeff(9) coeff(10) coeff(11) coeff(12) ];
			fitparamfunc = @(x,xv)(+x(1)*dircosh(+exp(-x(2)),xv,nt)+x(3)*dircosh(+exp(-x(2))+exp(-x(4)),xv,nt)+x(5)*dircosh(+exp(-x(2))+exp(-x(4))+exp(-x(6)),xv,nt)+x(7)*modcosh(+exp(-x(8)),xv,nt)+x(9)*modcosh(+exp(-x(8))+exp(-x(10)),xv,nt)+x(11)*modcosh(+exp(-x(8))+exp(-x(10))+exp(-x(12)),xv,nt));
			fitmassfunc = @(x,xv)(+x(1)*dircosh(x(2),xv,nt)+x(3)*dircosh(x(4),xv,nt)+x(5)*dircosh(x(6),xv,nt)+x(7)*modcosh(x(8),xv,nt)+x(9)*modcosh(x(10),xv,nt)+x(11)*modcosh(x(12),xv,nt));
            %chisqfunc = @(x)(1.0/(size(yval,2)-numel(guess))*(yval(:)-(+x(1)*cosh((+exp(-x(2)))*(nt/2-xval(:)))+x(3)*cosh((+exp(-x(2))+exp(-x(4)))*(nt/2-xval(:)))+x(5)*cosh((+exp(-x(2))+exp(-x(4))+exp(-x(6)))*(nt/2-xval(:)))+x(7)*(1-2*mod(xval(:),2)).*cosh((+exp(-x(8)))*(nt/2-xval(:)))+x(9)*(1-2*mod(xval(:),2)).*cosh((+exp(-x(8))+exp(-x(10)))*(nt/2-xval(:)))+x(11)*(1-2*mod(xval(:),2)).*cosh((+exp(-x(8))+exp(-x(10))+exp(-x(12)))*(nt/2-xval(:)))))'*(ycorr\(yval(:)-(+x(1)*cosh((+exp(-x(2)))*(nt/2-xval(:)))+x(3)*cosh((+exp(-x(2))+exp(-x(4)))*(nt/2-xval(:)))+x(5)*cosh((+exp(-x(2))+exp(-x(4))+exp(-x(6)))*(nt/2-xval(:)))+x(7)*(1-2*mod(xval(:),2)).*cosh((+exp(-x(8)))*(nt/2-xval(:)))+x(9)*(1-2*mod(xval(:),2)).*cosh((+exp(-x(8))+exp(-x(10)))*(nt/2-xval(:)))+x(11)*(1-2*mod(xval(:),2)).*cosh((+exp(-x(8))+exp(-x(10))+exp(-x(12)))*(nt/2-xval(:)))))));
			otherwise % do case 1.
                wut = 'wut'
        end
	
		% Get #dof. This is size(yval,2)-numel(guess)+ number of
		% things with a constraint ~1e-20.
		dof = size(yval,2)-numel(guess);
		if (exist('constraints','var'))
			% = 1 if an exact constraint, = 0 otherwise. 
			constr_flags = floor(((constraints < 1e-18) + (constraints > 1e-22))/1.5);
			dof = dof+sum(constr_flags);
		end
	
		% Build the chisq function!
		% This depends on if we need to SVD or not.
		if (fitcut == 0)
			chisqfunc = @(x)(1.0/(dof)*(yval(:)-fitparamfunc(x,xval(:)))'*(ycorr\(yval(:)-fitparamfunc(x,xval(:)))));
		else % cut things!
			[U, S, V] = svd(ycorr);
			
			% S contains the singular values in decreasing order.
			num_vals_save = size(yval, 2) - fitcut; 
            if ((dof-fitcut) <= 0) % fit is meaningless!
                fit_output = [];
                return;
            end
            
			for i=1:num_vals_save % invert first batch.
				S(i,i) = 1/S(i,i);
			end
			for i=(num_vals_save+1):size(yval,2)
				S(i,i) = 0; % remove the rest.
			end
			% rebuild cov inverse.
			ycorr_inv = U*S*V';
			
			% fix degrees of freedom.
			dof = dof-fitcut;
			
			% build chisq.
			chisqfunc = @(x)(1.0/(dof)*(yval(:)-fitparamfunc(x,xval(:)))'*(ycorr_inv*(yval(:)-fitparamfunc(x,xval(:)))));
			
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
				[temp_guess, fit_output(count,15), succ_code] = fminsearch(chisqfunc_constr, guess, optimset('TolX', fit_xprec,'TolFun', 1e-7,'MaxFunEvals',100000,'MaxIter',100000));
			else
				bound = max([abs(max(corr_fcn)) abs(min(corr_fcn))]);
				[temp_guess, fit_output(count,15)] = fminbnd(chisqfunc_constr, -bound, bound, optimset('TolX', fit_xprec,'MaxFunEvals',100000,'MaxIter',100000));
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
            
            
            fit_output(count,16) = 1.0-chi2cdf(fit_output(count,15)*dof, dof);
			
			% Get classical errors!
			%{
			undo_logmap = zeros(numel(temp_guess), 1);
			undo_logmap(1:(2*fitcosh)) = fit_output(3:(2+2*fitcosh));
			undo_logmap((2*fitcosh+1):(2*fitcosh+2*fitoscil)) = fit_output(9:(8+2*fitoscil));
			
			% This uses the method in:
			% http://www.mathworks.com/matlabcentral/newsreader/view_thread/157530
			
			% jacobian matrix
			cerr_function = @(x)(fitmassfunc(x,xval(:)));
			J = jacobianest(cerr_function,undo_logmap);
			
			% Get cov of fit parameters.
			% fit_output(count,15) is the chisq/dof.
			Sigma = fit_output(count,15)*

			% I'll be lazy here, and use inv. Please, no flames,
			% if you want a better approach, look in my tips and
			% tricks doc.
			Sigma = sdr^2*inv(J'*J);

			% Parameter standard errors
			se = sqrt(diag(Sigma))'
			%}

			% End scene!
			fit_output = fit_output(1,:);
        else
            fit_output = [];
        end
    
    else
        fit_output = [];
    end
    
    clear('choice');
    
end