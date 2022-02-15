%square positive defined matrix
max_size = 10;
max_num = 10;
mtx_cell = cell((max_size-1)*max_num,1);
mtx_inv_cell = cell((max_size-1)*max_num,1);
mtx_chol_up_cell = cell((max_size-1)*max_num,1);
mtx_chol_low_cell = cell((max_size-1)*max_num,1);
file_cell = {};
for i=2:max_size %matrix size i x i
    for j = 1:max_num %matrix number
        rr = rand(i);
        mtx_sym = rr*rr.';
        mtx_sym_name = ['mtx_' num2str(j) '_sym_' num2str(i) 'x' num2str(i)];
        
        mtx_inv = inv(mtx_sym);
        mtx_name_inv = [mtx_sym_name '_inv'];
        
        mtx_chol_low = chol(mtx_sym,'lower');
        mtx_name_chol_low = [mtx_sym_name '_chol_low'];
        
        mtx_chol_up = chol(mtx_sym,'lower');
        mtx_name_chol_up = [mtx_sym_name '_chol_up'];
            
        mtx_cell(j+(i-2)*max_num,:) = {mtx2arr(mtx_sym,mtx_sym_name)};
        mtx_inv_cell(j+(i-2)*max_num,:) = {mtx2arr(mtx_inv,mtx_name_inv)};
        mtx_chol_low_cell(j+(i-2)*max_num,:) = {mtx2arr(mtx_chol_low,mtx_name_chol_low)};
        mtx_chol_up_cell(j+(i-2)*max_num,:) = {mtx2arr(mtx_chol_low,mtx_name_chol_up)};
    end
end
file_cell = [mtx_cell; mtx_inv_cell; mtx_chol_low_cell; mtx_chol_up_cell];
%%%
obj_cell = cell((max_size-1)*max_num,1);
obj_inv_cell = cell((max_size-1)*max_num,1);
obj_chol_up_cell = cell((max_size-1)*max_num,1);
obj_chol_low_cell = cell((max_size-1)*max_num,1);

allobj_cell = {'Matrix_t *allSym[] = {'};
invobj_cell = {'Matrix_t *invSym[] = {'};
lcholobj_cell = {'Matrix_t *loCholSym[] = {'};
ucholobj_cell = {'Matrix_t *upCholSym[] = {'};
for i=2:max_size %matrix size i x i
    for j = 1:max_num %matrix number
        mtx_sym_name = ['mtx_' num2str(j) '_sym_' num2str(i) 'x' num2str(i)];
        mtx_name_inv = [mtx_sym_name '_inv'];
        mtx_name_chol_low = [mtx_sym_name '_chol'];
        
        obj_name = ['obj_' num2str(j) '_sym_' num2str(i) 'x' num2str(i)];
        obj_name_inv = [obj_name '_inv'];
        obj_name_chol_low = [obj_name '_chol_low'];
        obj_name_chol_up = [obj_name '_chol_up'];
            
        obj_cell(j+(i-2)*max_num,:) = {['Matrix_t ' obj_name '={' num2str(i*i) ',' num2str(i) ',' num2str(i) ',&' mtx_sym_name '[0][0]};'; ]};
        obj_inv_cell(j+(i-2)*max_num,:) = {['Matrix_t ' obj_name_inv '={' num2str(i*i) ',' num2str(i) ',' num2str(i) ',&' mtx_name_inv '[0][0]};'; ]};
        obj_chol_low_cell(j+(i-2)*max_num,:) = {['Matrix_t ' obj_name_chol_low '={' num2str(i*i) ',' num2str(i) ',' num2str(i) ',&' mtx_name_chol_low '[0][0]};'; ]};
        obj_chol_up_cell(j+(i-2)*max_num,:) = {['Matrix_t ' obj_name_chol_up '={' num2str(i*i) ',' num2str(i) ',' num2str(i) ',&' mtx_name_chol_up '[0][0]};'; ]};
        
        allobj_cell(j+(i-2)*max_num+1,:) = {['&' obj_name ','; ]};
        invobj_cell(j+(i-2)*max_num+1,:) = {['&' obj_name_inv ','; ]};
        lcholobj_cell(j+(i-2)*max_num+1,:) = {['&' obj_name_chol_low ','; ]};
        ucholobj_cell(j+(i-2)*max_num+1,:) = {['&' obj_name_chol_up ','; ]};
    end
end
allobj_cell(j+(i-2)*max_num+2,:) = {'};'};
invobj_cell(j+(i-2)*max_num+2,:) = {'};'};
lcholobj_cell(j+(i-2)*max_num+2,:) = {'};'};
ucholobj_cell(j+(i-2)*max_num+2,:) = {'};'};

file_cell = [file_cell; 
              obj_cell; 
              obj_inv_cell; 
              obj_chol_low_cell;
              obj_chol_up_cell; 
              allobj_cell; 
              invobj_cell; 
              lcholobj_cell; 
              ucholobj_cell];
fidS = fopen('mtxGen.c','w');
fprintf(fidS,'%s\n',file_cell{:});  
fclose(fidS);
