% Given a set of operators that have a vev, build the vev-subtracted correlator.
% This helps streamline jackknife analysis. This also forms a full variational
% correlation basis, which can be used for glueballs or to cheat combining
% different <pbp> noise sources.
function disc_corrs = build_vev_correlator(ops, fold)

	parse_Nt = size(ops, 1);
	number_bl = size(ops, 2);
	num_data_in = size(ops, 3);

	% First, get a vev. This is a per-noise-source ensemble average.
	vev_reduce = mean(mean(ops,1),3);
	vev_to_save = mean(mean(ops,1),2);
	vev_reduce_resize = zeros(1, num_data_in);
    vev_reduce_resize(:,:) = vev_to_save(1,:,:);
    
	vev_reduce_copy = repmat(vev_reduce, [parse_Nt 1 num_data_in]);
	%ops = ops - vev_reduce_copy; % Don't vev subtract yet, because of blocking.
	
    ops_P = fft(ops, [], 1);
    ops_Pdag = conj(ops_P);

    % Preallocate
    corr_P = zeros(parse_Nt, number_bl, number_bl, num_data_in);

    % Some loops are inevitable...
    for i=1:number_bl
        for j=1:number_bl
			if (fold == 0)
				corr_P(:,i,j,:) = ops_P(:,i,:).*ops_Pdag(:,j,:);
			else
				corr_P(:,i,j,:) = real(ops_P(:,i,:).*ops_Pdag(:,j,:));
			end
        end
    end

    clear('ops_P'); clear('ops_Pdag');

	% Don't contract on different sources!
	disc_nvevsub = ifft(corr_P, [], 1);
	
	% Don't average over all configs.
	%disc_corrs = mean(disc_nvevsub, 4);
    disc_corrs = disc_nvevsub;
	
end

