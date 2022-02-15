%square positive defined matrix
max_size = 10;
max_num = 10;
mtx_cell = cell((max_size-1)*max_num,1);
mtx_inv_cell = cell((max_size-1)*max_num,1);
mtx_chol_up_cell = cell((max_size-1)*max_num,1);
mtx_chol_low_cell = cell((max_size-1)*max_num,1);
mtx_magic_cell = cell(max_size-1,1);
mtx_magic_transp_cell = cell(max_size-1,1);
file_cell = {};
for i=2:max_size %matrix size i x i
    mtx_magic = magic(i);
    mtx_magic_name = ['mtx_1_magic_' num2str(i) 'x' num2str(i)];
    
    mtx_magic_transp = transp(mtx_magic);
    mtx_magic_transp_name = [mtx_magic_name '_transp'];
    
    mtx_magic_cell(i-1,:) = {mtx2arr(mtx_magic,mtx_magic_name)};
    mtx_magic_transp_cell(i-1,:) = {mtx2arr(mtx_magic_transp,mtx_magic_transp_name)};
        
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
file_cell = [mtx_cell; mtx_inv_cell; mtx_chol_low_cell; mtx_chol_up_cell;mtx_magic_cell;mtx_magic_transp_cell];
%%%
obj_cell = cell((max_size-1)*max_num,1);
obj_inv_cell = cell((max_size-1)*max_num,1);
obj_chol_up_cell = cell((max_size-1)*max_num,1);
obj_chol_low_cell = cell((max_size-1)*max_num,1);
obj_magic_cell = cell(max_size-1,1);
obj_magic_transp_cell = cell(max_size-1,1);

allobj_cell = {'Matrix_t *allSym[] = {'};
invobj_cell = {'Matrix_t *invSym[] = {'};
lcholobj_cell = {'Matrix_t *loCholSym[] = {'};
ucholobj_cell = {'Matrix_t *upCholSym[] = {'};
magicobj_cell = {'Matrix_t *allMagic[] = {'};
magictrobj_cell = {'Matrix_t *transpMagic[] = {'};
for i=2:max_size %matrix size i x i
    mtx_magic_name = ['mtx_1_magic_' num2str(i) 'x' num2str(i)];
    mtx_magic_transp_name = [mtx_magic_name '_transp'];
    
    obj_magic_name = ['obj_1_magic_' num2str(i) 'x' num2str(i)];
    obj_magic_transp_name = [obj_magic_name '_transp'];
    
    obj_magic_cell(i-1,:) = {['Matrix_t ' obj_magic_name '={' num2str(i*i) ',' num2str(i) ',' num2str(i) ',&' mtx_magic_name '[0][0]};'; ]};
    obj_magic_transp_cell(i-1,:) = {['Matrix_t ' obj_magic_transp_name '={' num2str(i*i) ',' num2str(i) ',' num2str(i) ',&' mtx_magic_transp_name '[0][0]};'; ]};
    
    magicobj_cell(i,:) = {['&' obj_magic_name ','; ]};
    magictrobj_cell(i,:) = {['&' obj_magic_transp_name ','; ]};
    
    for j = 1:max_num %matrix number

        mtx_sym_name = ['mtx_' num2str(j) '_sym_' num2str(i) 'x' num2str(i)];
        mtx_sym_inv_name = [mtx_sym_name '_inv'];
        mtx_sym_chol_low_name = [mtx_sym_name '_chol_low'];
        mtx_sym_chol_up_name = [mtx_sym_name '_chol_up'];
        
        obj_sym_name = ['obj_' num2str(j) '_sym_' num2str(i) 'x' num2str(i)];
        obj_sym_inv_name = [obj_sym_name '_inv'];
        obj_sym_chol_low_name = [obj_sym_name '_chol_low'];
        obj_sym_chol_up_name = [obj_sym_name '_chol_up'];
                    
        obj_cell(j+(i-2)*max_num,:) = {['Matrix_t ' obj_sym_name '={' num2str(i*i) ',' num2str(i) ',' num2str(i) ',&' mtx_sym_name '[0][0]};'; ]};
        obj_inv_cell(j+(i-2)*max_num,:) = {['Matrix_t ' obj_sym_inv_name '={' num2str(i*i) ',' num2str(i) ',' num2str(i) ',&' mtx_sym_inv_name '[0][0]};'; ]};
        obj_chol_low_cell(j+(i-2)*max_num,:) = {['Matrix_t ' obj_sym_chol_low_name '={' num2str(i*i) ',' num2str(i) ',' num2str(i) ',&' mtx_sym_chol_low_name '[0][0]};'; ]};
        obj_chol_up_cell(j+(i-2)*max_num,:) = {['Matrix_t ' obj_sym_chol_up_name '={' num2str(i*i) ',' num2str(i) ',' num2str(i) ',&' mtx_sym_chol_up_name '[0][0]};'; ]};
        
        allobj_cell(j+(i-2)*max_num+1,:) = {['&' obj_sym_name ','; ]};
        invobj_cell(j+(i-2)*max_num+1,:) = {['&' obj_sym_inv_name ','; ]};
        lcholobj_cell(j+(i-2)*max_num+1,:) = {['&' obj_sym_chol_low_name ','; ]};
        ucholobj_cell(j+(i-2)*max_num+1,:) = {['&' obj_sym_chol_up_name ','; ]};
    end
end
allobj_cell(j+(i-2)*max_num+2,:) = {'};'};
invobj_cell(j+(i-2)*max_num+2,:) = {'};'};
lcholobj_cell(j+(i-2)*max_num+2,:) = {'};'};
ucholobj_cell(j+(i-2)*max_num+2,:) = {'};'};
magicobj_cell(i+1,:) = {'};'};
magictrobj_cell(i+1,:) = {'};'};
file_cell = [file_cell; 
              obj_cell; 
              obj_inv_cell; 
              obj_chol_low_cell;
              obj_chol_up_cell;
              obj_magic_cell;
              obj_magic_transp_cell
              allobj_cell; 
              invobj_cell; 
              lcholobj_cell; 
              ucholobj_cell;
              magicobj_cell;
              magictrobj_cell];
fidS = fopen('mtxGen.c','w');
fprintf(fidS,'%s\n',file_cell{:});  
fclose(fidS);
