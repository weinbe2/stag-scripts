% Positive parity project data, assuming dim1 is time.
function [out_data] = fold_data(in_data)
	parse_Nt = size(in_data, 1);
	
	connected = in_data;
	
	tmp_meh = connected(1,:); % last term gets special treatment
	for tmp_i = 0:(parse_Nt-3)
		connected(tmp_i+1,:) = connected(tmp_i+1,:)+2*connected(tmp_i+2,:)+connected(tmp_i+3,:);
	end
	connected(parse_Nt-1,:) = connected(parse_Nt-1,:)+2*connected(parse_Nt,:)+tmp_meh(1,:);
	clear('tmp_meh');
	% Shift things forward.
	for tmp_i = (parse_Nt):-1:2
		connected(tmp_i,:) = 0.25*connected(tmp_i-1,:);
	end
	connected(1,:) = zeros(1, size(in_data,2)); % we lose the very edge.
	
	out_data = connected;
	
end