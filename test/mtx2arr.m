function arr = mtx2arr(mtx,name)
% convert numeric matrix in string ready for initialization of C arrays
% Example:
% in = [1 2 3; 4 5 6]
% out = {{1,2,3,},{4,5,6}}

[r,c]=size(mtx);
arr = arrayfun(@num2str, mtx,repmat(20,r,c), 'UniformOutput', false);
if c > 1    
    arr = cellfun(@strcat,arr, repmat({','},r,c),'UniformOutput',false);
    arr(:,1) = cellfun(@strcat,repmat({'{'},r,1),arr(:,1),'UniformOutput',false);
    arr(:,end) = cellfun(@strcat,arr(:,end),repmat({'},'},r,1),'UniformOutput',false);
else
    arr(:,1) = cellfun(@strcat,repmat({'{'},r,1),arr(:,1),'UniformOutput',false);
    arr(:,end) = cellfun(@strcat,arr(:,end),repmat({'},'},r,1),'UniformOutput',false);
end
arr = num2cell(arr,1);
arr = strcat(arr{:});
arr = strjoin(arr');
arr = ['{' arr(1:end-1) '}'];
arr = arr(~isspace(arr));
arr = regexprep(arr,',}}','}}');
arr = ['const double ' name '[' num2str(r) '][' num2str(c) ']' ' = ' arr ';' ];
end
