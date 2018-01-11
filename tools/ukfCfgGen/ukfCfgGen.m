%Initial version of UKF cfg generator
%Status : In progress
clear all
clc
%Initialization section (test with pendulum)
cfgID = 1;
dT = 0.0001;

discreteStateFcn = {'x(1) = x(1) + dT*x(2)';
                    'x(2) = (1 - dT*0.1)*x(2) - dt*16.003263*sin(x(1))'};

measFcn = {'y(1) = x(1)'};

[xL,~] = size(discreteStateFcn);
[yL,~] = size(measFcn);
sL = 2*xL+1;
uL= 2;

alpha = 1;
betha = 2;
kappa = 0;

Ryy0 = zeros(yL,yL);
Pxx0 = diag(ones(1,xL)*0.01);
Qxx = diag(ones(1,xL)*0.01);

%UKF matrix properties
ukfMatrix={
% name         ,size ,value,             ,presence
{'<Sc_vector>' ,[1,3],[alpha,betha,kappa],true}
{'<Wm_weight_vector>' ,[1,sL],zeros(1,sL),true}
{'<Wc_weight_vector>' ,[1,sL],zeros(1,sL),true}
{'<u_system_input>' ,[uL,1],zeros(uL,1),false}
{'<u_prev_system_input>' ,[uL,1],zeros(uL,1),false}
{'<y_meas>' ,[yL,1],zeros(yL,1),true}
{'<y_predicted_mean>' ,[yL,1],zeros(yL,1),true}
{'<x_system_states>',[xL,1],zeros(xL,1),true}
{'<x_system_states_ic>' ,[xL,1],zeros(xL,1),true}
{'<x_system_states_limits>' ,[xL,3],zeros(xL,3),true}
{'<x_system_states_limits_enable>' ,[xL,1],zeros(xL,1),true}
{'<x_system_states_correction>' ,[xL,1],zeros(xL,1),true}
{'<X_sigma_points>' ,[xL,sL],zeros(xL,sL),true}
{'<Y_sigma_points>' ,[yL,sL],zeros(yL,1),true}
{'<Pxx_error_covariance>' ,[xL,xL],zeros(xL,xL),true}
{'<Pxx0_init_error_covariance>' ,[xL,xL],Pxx0,true}
{'<Qxx_process_noise_cov>' ,[xL,xL],Qxx,true}
{'<Ryy0_init_out_covariance>' ,[yL,yL],Ryy0,true}
{'<Pyy_out_covariance>' ,[yL,yL],zeros(yL,yL),true}
{'<Pyy_out_covariance_copy>' ,[yL,yL],zeros(yL,yL),true}
{'<Pxy_cross_covariance>' ,[xL,yL],zeros(xL,yL),true}
{'<K_kalman_gain>' ,[xL,yL],zeros(xL,yL),true}
{'<Pxx_covariance_correction>' ,[xL,xL],zeros(xL,xL),true}
{'<I_identity_matrix>',[yL,yL],zeros(yL,yL),true}
}



sourceStateFcn = cell(xL,1);
for mIdx=1:xL
    sourceStateFcn(mIdx) = discreteStateFcn(mIdx);
    
    eqIdx = strfind(sourceStateFcn{mIdx},'=');
    sourceStateFcn(mIdx) = {[strrep(sourceStateFcn{mIdx}(1:eqIdx),['x(' num2str(mIdx) ')' ], [ 'pX_m->val[nCol*' num2str(mIdx-1) '+sigmaIdx]' ]) sourceStateFcn{mIdx}(eqIdx+1:end)]};
    
    for nIdx=1:xL
        sourceStateFcn(mIdx) = {['  ' strrep(sourceStateFcn{mIdx},['x(' num2str(nIdx) ')' ], [ 'pX_p->val[nCol*' num2str(nIdx-1) '+sigmaIdx]' ])]};
    end;
end

sourceMeasFcn = cell(yL,1);
for pIdx=1:yL
    sourceMeasFcn(pIdx) = measFcn(pIdx);
    
    eqIdx = strfind(sourceMeasFcn{pIdx},'=');
    sourceMeasFcn(pIdx) = {[strrep(sourceMeasFcn{pIdx}(1:eqIdx),['y(' num2str(pIdx) ')' ], [ 'pY_m->val[nCol*' num2str(pIdx-1) '+sigmaIdx]' ]) sourceMeasFcn{pIdx}(eqIdx+1:end)]};
    
    for qIdx=1:yL
        sourceMeasFcn(pIdx) = {['    ' strrep(sourceMeasFcn{pIdx},['x(' num2str(qIdx) ')' ], [ 'pX_m->val[nCol*' num2str(qIdx-1) '+sigmaIdx]' ])]};
    end;
end

newCfgSource = ['ukfCfg' num2str(cfgID) '.c'];
newCfgHeader = ['ukfCfg' num2str(cfgID) '.h'];

[fidS,msg] = fopen('ukfCfgTemplate.c');

str = textscan(fidS,'%s', 'delimiter', '\n');

c = str{1};
%put header ID
tmp1 = cell(1,length(c))';
tmp1(:) = {'<cfgId>'};

tmp2 = cell(1,length(c))';
tmp2(:) = {num2str(cfgID)};
c = cellfun(@strrep,c,tmp1,tmp2,'UniformOutput',false);

ptrString = 'static tPredictFcn PredictFcn[xL] = {';
for i = 1:xL 
    %State transition prototype
    endIdx = find(~cellfun(@isempty,strfind(c, '<STATE TRANSITION PROTOTYPE:END>')));
    c(endIdx+1:end+1,:) = c(endIdx:end,:);
    c(endIdx,:) = {['static void Fx' num2str(i) '(tMatrix * pu_p, tMatrix * pX_p, tMatrix * pX_m,uint8 sigmaIdx, float64 dT);']};
    
    % state transition ptr array
    ptrString = [ptrString '&Fx' num2str(i) ','];
    
    %State transition body
    endIdx = find(~cellfun(@isempty,strfind(c, '<STATE TRANSITION:END>')));
    c(endIdx+1:end+1,:) = c(endIdx:end,:);
    c(endIdx,:) = {['void Fx' num2str(i) '(tMatrix * pu_p, tMatrix * pX_p, tMatrix * pX_m,uint8 sigmaIdx, float64 dT)']};
    
    endIdx = endIdx+1;
    c(endIdx+1:end+1,:) = c(endIdx:end,:);
    c(endIdx,:) = {'{'};
    
    endIdx = endIdx+1;
    c(endIdx+1:end+1,:) = c(endIdx:end,:);
    c(endIdx,:) = {['    ' 'const uint8 nCol = pX_m->ncol;']};
    
    endIdx = endIdx+1;
    c(endIdx+1:end+1,:) = c(endIdx:end,:);
    c(endIdx,:) = sourceStateFcn(xL);
    
    endIdx = endIdx+1;
    c(endIdx+1:end+1,:) = c(endIdx:end,:);
    c(endIdx,:) = {'}'}
end

ptrString = ptrString(1:end-1);
ptrString = [ptrString '};'];

endIdx = find(~cellfun(@isempty,strfind(c, '<STATE TRANSITION PTR ARRAY:END>')));
c(endIdx+1:end+1,:) = c(endIdx:end,:);
c(endIdx,:) = {ptrString};

ptrString = 'static tObservFcn  ObservFcn[yL] = {';
for j = 1:yL
    %measurement prototype
    endIdx = find(~cellfun(@isempty,strfind(c, '<MEASUREMENT PROTOTYPE:END>')));
    c(endIdx+1:end+1,:) = c(endIdx:end,:);
    c(endIdx,:) = {['static void Hy' num2str(j) '(tMatrix * pu, tMatrix * pX_m, tMatrix * pY_m,uint8 sigmaIdx);']};
    
    %measurement ptr array
     ptrString = [ptrString '&Hy' num2str(j) ','];
    
    %measurement body
    endIdx = find(~cellfun(@isempty,strfind(c, '<MEASUREMENT FUNCTION:END>')));
    c(endIdx+1:end+1,:) = c(endIdx:end,:);
    c(endIdx,:) = {['void Hy' num2str(j) '(tMatrix * pu, tMatrix * pX_m, tMatrix * pY_m,uint8 sigmaIdx)']};
    
    endIdx = endIdx+1;
    c(endIdx+1:end+1,:) = c(endIdx:end,:);
    c(endIdx,:) = {'{'};
    
    endIdx = endIdx+1;
    c(endIdx+1:end+1,:) = c(endIdx:end,:);
    c(endIdx,:) = {['    ' 'const uint8 nCol = pY_m->ncol;']}
    
    endIdx = endIdx+1;
    c(endIdx+1:end+1,:) = c(endIdx:end,:);
    c(endIdx,:) = sourceMeasFcn(j);
    
    endIdx = endIdx+1;
    c(endIdx+1:end+1,:) = c(endIdx:end,:);
    c(endIdx,:) = {'}'}   
end

ptrString = ptrString(1:end-1);
ptrString = [ptrString '};'];

endIdx = find(~cellfun(@isempty,strfind(c, '<MEASUREMENT PTR ARRAY:END>')))
c(endIdx+1:end+1,:) = c(endIdx:end,:);
c(endIdx,:) = {ptrString};

l=cell(length(c),1);
replace = cell(length(c),1);
for k = 1:length(ukfMatrix)
    if(true == cell2mat(ukfMatrix{k}(4)))
        l(:)=ukfMatrix{k}(1);
        v = ukfMatrix{k}(2);
        row = v{1}(1);
        col = v{1}(2);
        replace(:) = {mtx2carr(cell2mat(ukfMatrix{k}(3)))};      
        c = cellfun(@strrep, c, l,replace,'UniformOutput',false);
    else
        rowIdx =  cellfun(@strfind,c,repmat(ukfMatrix{k}(1),length(c),1),'UniformOutput',false);
        rowIdx = find(~cellfun(@isempty,rowIdx))
        c{rowIdx}(1:2)= '//'
    end
end

delete(newCfgSource);
fid = fopen(newCfgSource,'w');

for rowIdx=1:length(c)
fprintf(fid,'%s\n',c{rowIdx});    
end

