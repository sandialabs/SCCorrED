function patternCenter = getPatternCenter(fname)
    if ~isempty(regexp(fname, '.ang', 'ONCE'))
        patternCenter = getPatternCenter_ang(fname);
    end
end
%     patternCenter = [0.5, 0.5, 0.7];
%     if ~isempty(regexp(fname, '.ang', 'ONCE'))
%         fid = fopen(fname);
%         pc_found_logic = [0, 0, 0];
%         pc_labels = ["# x-star", "# y-star", "# z-star"];
%         for i = 1:50
%             line = fgetl(fid);
%             match_str = 0;
%             match_str = string(regexp(line, '# [xyz]-star', 'match'));
%             if ~isempty(match_str)
%                 ind = find(contains(pc_labels, match_str));
%                 patternCenter(ind) = str2double(string(line(25:end))); 
%                 pc_found_logic(ind) = 1;
%             end
%             if all(pc_found_logic)
%                 disp(i)
%                 break       
%             end
%         end
%     end
% end


function patternCenter = getPatternCenter_ang(fname) 
    patternCenter = [0.5, 0.5, 0.7];
    fid = fopen(fname);
    pc_found_logic = [0, 0, 0];
    pc_labels = ["# x-star", "# y-star", "# z-star"];
    for i = 1:50
        line = fgetl(fid);
        match_str = 0;
        match_str = string(regexp(line, '# [xyz]-star', 'match'));
        if ~isempty(match_str)
            ind = find(contains(pc_labels, match_str));
            patternCenter(ind) = str2double(string(line(25:end))); 
            pc_found_logic(ind) = 1;
        end
        if all(pc_found_logic)
            break       
        end
    end
end