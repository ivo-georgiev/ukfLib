%square positive defined matrix
max_size = 20;
max_num = 100;
mtx_cell = {};
k=1;c=1;
for i=2:max_size %matrix size i x i
    for j = 1:max_num %matrix number
        rr = rand(i);
        mtx = rr*rr.';
        mtx_name = ['mtx_' num2str(j) '_' num2str(i) 'x' num2str(i)];
        
        mtx_inv = inv(mtx);
        mtx_name_inv = [mtx_name '_inv'];
        
        mtx_chol = chol(mtx);
        mtx_name_chol = [mtx_name '_chol'];
            
        mtx_cell(j+(i-2)*max_num+(k-1)+(c-1),:) = {mtx2arr(mtx,mtx_name)};
        mtx_cell(j+(i-2)*max_num+k+(c-1),:) = {mtx2arr(mtx_inv,mtx_name_inv)};
        mtx_cell(j+(i-2)*max_num+k+c,:) = {mtx2arr(mtx_chol,mtx_name_chol)};
        offset = j+(i-2)*max_num+k+c;
        k=k+1;
        c=c+1;
    end
end
%%%
k=1;c=1;
for i=2:max_size %matrix size i x i
    for j = 1:max_num %matrix number
        mtx_name = ['mtx_' num2str(j) '_' num2str(i) 'x' num2str(i)];
        mtx_name_inv = [mtx_name '_inv'];
        mtx_name_chol = [mtx_name '_chol'];
        
        obj_name = ['obj_' num2str(j) '_' num2str(i) 'x' num2str(i)];
        obj_name_inv = [obj_name '_inv'];
        obj_name_chol = [obj_name '_chol'];
            
        mtx_cell(offset+j+(i-2)*max_num+(k-1)+(c-1),:) = {['Matrix_t ' obj_name '={' num2str(i*i) ',' num2str(i) ',' num2str(i) ',&' mtx_name '[0][0]};'; ]};
        mtx_cell(offset+j+(i-2)*max_num+k+(c-1),:) = {['Matrix_t ' obj_name_inv '={' num2str(i*i) ',' num2str(i) ',' num2str(i) ',&' mtx_name_inv '[0][0]};'; ]};
        mtx_cell(offset+j+(i-2)*max_num+k+c,:) = {['Matrix_t ' obj_name_chol '={' num2str(i*i) ',' num2str(i) ',' num2str(i) ',&' mtx_name_chol '[0][0]};'; ]};
        
        k=k+1;
        c=c+1;
    end
end
fidS = fopen('mtxGen.c','w');
fprintf(fidS,'%s\n',mtx_cell{:});  
fclose(fidS);
