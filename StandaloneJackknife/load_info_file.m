function [fl_l, fl_h, parse_Ns, parse_Nt, beta, m_l, m_h] = load_info_file(fname)
    % This is a modified version which loads info.txt.
    full_fname = strcat(fname, '/info.txt');
    info_data = importdata(full_fname);
    
    % The first item is the 'flavor'.
    flv_array = strsplit(info_data{1}, '->');
    switch flv_array{2}
        case '4plus8'
            fl_l = 4;
            fl_h = 8;
            tmp_array = strsplit(info_data{2}, '->');
            parse_Ns = str2num(tmp_array{2});
            tmp_array = strsplit(info_data{3}, '->');
            parse_Nt = str2num(tmp_array{2});
            tmp_array = strsplit(info_data{4}, '->');
            beta = str2num(tmp_array{2});
            tmp_array = strsplit(info_data{5}, '->');
            m_l = str2num(tmp_array{2});
            tmp_array = strsplit(info_data{6}, '->');
            m_h = str2num(tmp_array{2});
            
            
        case '4'
            fl_l = 4;
            fl_h = 0;
            tmp_array = strsplit(info_data{2}, '->');
            parse_Ns = str2num(tmp_array{2});
            tmp_array = strsplit(info_data{3}, '->');
            parse_Nt = str2num(tmp_array{2});
            tmp_array = strsplit(info_data{4}, '->');
            beta = str2num(tmp_array{2});
            tmp_array = strsplit(info_data{5}, '->');
            m_l = str2num(tmp_array{2});
            tmp_array = strsplit(info_data{6}, '->');
            m_h = 0;
            
        case '8'
            fl_l = 8;
            fl_h = 0;
            tmp_array = strsplit(info_data{2}, '->');
            parse_Ns = str2num(tmp_array{2});
            tmp_array = strsplit(info_data{3}, '->');
            parse_Nt = str2num(tmp_array{2});
            tmp_array = strsplit(info_data{4}, '->');
            beta = str2num(tmp_array{2});
            tmp_array = strsplit(info_data{5}, '->');
            m_l = str2num(tmp_array{2});
            tmp_array = strsplit(info_data{6}, '->');
            m_h = 0;
            
        case '12'
            fl_l = 12;
            fl_h = 0;
            tmp_array = strsplit(info_data{2}, '->');
            parse_Ns = str2num(tmp_array{2});
            tmp_array = strsplit(info_data{3}, '->');
            parse_Nt = str2num(tmp_array{2});
            tmp_array = strsplit(info_data{4}, '->');
            beta = str2num(tmp_array{2});
            tmp_array = strsplit(info_data{5}, '->');
            m_l = str2num(tmp_array{2});
            tmp_array = strsplit(info_data{6}, '->');
            m_h = 0;
    end

end