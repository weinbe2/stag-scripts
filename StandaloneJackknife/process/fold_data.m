% Fold data, assuming dim1 is time.
function [out_data] = fold_data(in_data, is_stag_baryon, is_sinh)
	parse_Nt = size(in_data, 1);
	
	connected = in_data;
	
	if (~exist('is_stag_baryon', 'var'))
		is_stag_baryon = 0;
	end
	
	if (~exist('is_sinh', 'var'))
		is_sinh = 0;
	end
	
	for t=2:(parse_Nt/2)
		if (is_stag_baryon == 1)
			connected(t,:) = (connected(t,:)+(1-2.*repmat(mod(t,2), [1 size(connected, 2)])).*connected((parse_Nt+2)-t,:))/2;
		elseif (is_sinh == 1)
			connected(t,:) = (connected(t,:)-connected((parse_Nt+2)-t,:))/2;
		else
			connected(t,:) = (connected(t,:)+connected((parse_Nt+2)-t,:))/2;
		end
	end
	for t=(parse_Nt/2+2):parse_Nt
		if (is_stag_baryon == 1)
			connected(t,:) = (1-2.*repmat(mod(t,2), [1 size(connected, 2)])).*connected((parse_Nt+2)-t,:);
		elseif (is_sinh == 1)
			connected(t,:) = -connected((parse_Nt+2)-t,:);
		else
			connected(t,:) = connected((parse_Nt+2)-t,:);
		end
	end
	
	out_data = connected;
	
end