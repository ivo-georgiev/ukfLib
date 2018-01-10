function str = mtx2carr(mtx)
% convert numeric matrix in string ready for initialization of C arrays
% Example:
% in = [1 2 3; 4 5 6]
% out = {{1,2,3,},{4,5,6}}

[r,c]=size(mtx);
str = arrayfun(@num2str, mtx, 'UniformOutput', false);
if c > 1    
    str = cellfun(@strcat,str, repmat({','},r,c),'UniformOutput',false);
    str = [cellfun(@strcat,repmat({'{'},r,1),str(:,1),'UniformOutput',false) str(:,2:end)];
    str = [str(:,1:end-1) cellfun(@strcat,str(:,end),repmat({'},'},r,1),'UniformOutput',false)];
else
    str = [cellfun(@strcat,repmat({'{'},r,1),str(:,1),'UniformOutput',false) str(:,2:end)];
    str = [str(:,1:end-1) cellfun(@strcat,str(:,end),repmat({'},'},r,1),'UniformOutput',false)];
end;
str = num2cell(str,1);
str = strcat(str{:});
str = strjoin(str');
str = ['{' str(1:end-1) '}'];
str = str(~isspace(str));
str = regexprep(str,',}}','}}');
end