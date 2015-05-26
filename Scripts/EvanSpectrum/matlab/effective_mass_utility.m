function [masses, root, amps] = effective_mass_utility(connected_sum, parse_Nt, K, N, C)
	
    % Okay. Let's try this new effective mass method. This form assumes 2
    % effective masses.

    %K = 2; % We want two states.
    M = K+1; % Don't look too far ahead. M>K
    %N = 20; % N >= M+K-1.
    %C = 0; % number to SVD cut. C<K.


    % Loop over all t.

    masses = zeros(parse_Nt-2*N+1, K);
    root = zeros(parse_Nt-2*N+1, K);
    amps = zeros(parse_Nt-2*N+1, K);

    for t=(N+1):(parse_Nt-N+1)

        % We need N y_n's.

        y = zeros(N, 1);

        for n=1:N
            for j=0:(n-1)
                y(n) = y(n) + 2.^(-n+1).*nchoosek(n-1,j).*connected_sum(t+n-2*j-1);
            end
        end

        % Okay. We maybe have all of the n's.

        % Next, let's set up our matrices.
        hlp = zeros(N-M+1, 1);

        hlp(1:(N-M+1)) = y(M:N);

        Hlp = zeros(N-M+1, K);
        for k=1:K
            Hlp(:,k) = -y((M-k):(N-k));
        end

        % Look at SVD cut?
        [U S V] = svd(Hlp);

        for k=(K-C+1):K
            S(k,k) = 0;
        end

        Hlp = U*S*V';

        % Solve the least squares problem hlp = Hlp_prime p for p. 
        % We do this with MATLAB's \ divide.

        p = Hlp \ hlp;

        % Construct the polynomial.

        p_prime = zeros(K+1, 1);
        p_prime(1) = 1;
        p_prime(2:(K+1)) = p;

        r = sort(roots(p_prime));
        root(t-N,:) = r;

        for k=1:K

            if (r(k) > 1) % cosh
                masses(t-N,k) = acosh(r(k));
            elseif (r(k) < -1) % oscillating cosh
                masses(t-N,k) = acosh(-r(k));
            else % something
                masses(t-N,k) = acosh(r(k)); % it'll come out cmplx.
            end
        end

        % Next, we can get amplitudes via y = X a
        Xmat = zeros(N, K);
        yvec = y;

        for n=1:N
            Xmat(n,:) = r.^n;
        end

        amp = Xmat \ yvec;
        amps(t-N,:) = (amp'./cosh(masses(t-N,:)*(parse_Nt/2+1-t)));

    end

end