%Initial version of UKF cfg generator
%Status : In progress
clear all
clc
%Initialization section (test with pendulum)
cfgID = 1;
dT = 0.0001;

discreteStateFcn = {'x(1) = x(1) + dT*x(2)';
                    'x(2) = (1 - dT*0.1)*x(2) - dt*16.003263*sin(x(1))'}

[xL,~] = size(discreteStateFcn);

sourceStateFcn = cell(xL,1);
for mIdx=1:xL
    sourceStateFcn(mIdx) = discreteStateFcn(mIdx)
    
    eqIdx = strfind(sourceStateFcn{mIdx},'=');
    sourceStateFcn(mIdx) = {[strrep(sourceStateFcn{mIdx}(1:eqIdx),['x(' num2str(mIdx) ')' ], [ 'pX_m->val[nCol*' num2str(mIdx-1) '+ sigmaIdx]' ]) sourceStateFcn{mIdx}(eqIdx+1:end)]}
    
    for nIdx=1:xL
        sourceStateFcn(mIdx) = {strrep(sourceStateFcn{mIdx},['x(' num2str(nIdx) ')' ], [ 'pX_p->val[nCol*' num2str(nIdx-1) '+sigmaIdx]' ])}
    end;
end

newCfgSource = ['ukfCfg' num2str(cfgID) '.c']
newCfgHeader = ['ukfCfg' num2str(cfgID) '.h']

copyfile('ukfCfgTemplate.c', newCfgSource)
[fidS,msg] = fopen(newCfgSource)

str = textscan(fidS,'%s', 'delimiter', '\n')

%search for state transition section
%beginStTrIdx = find(strcmp(str{1}, '<STATE TRANSITION:BEGIN>'), 1, 'first');



