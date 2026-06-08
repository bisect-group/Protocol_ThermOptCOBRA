function [a,modModel] = ThermOptCC(model,tol,getDirection)
% Identifies thermodynamically feasible flux directions for all the
% reactions in the input model
%
% USAGE: 
%   [a,modModel] = ThermOptCC(model,tol,getDirection)
%
% INPUTS:
%     model:        COBRA model structure for which thermodynamic feasibility
%                   of the reactions has to be identified
%     tol:          Tolerance value (User defined non-zero value).
%     getDirection: (Optional) Boolean. If true, returns detailed direction.
%                   If false, returns 'Blocked' or 'Unblocked' and optimizes execution.
%
% OUTPUTS:
%     a:          A cell describing the thermodynamically feasible direction 
%                 of the reactions in the given input model                 
%     modModel:   Modified model that has no irreversible reactions that
%                 carry flux in reverse direction.
%
% .. Author:
%       - Pavan Kumar S, BioSystems Engineering and control (BiSECt) lab, IIT Madras

if nargin < 3 || isempty(getDirection)
    getDirection = false;
end

% checking for reactions irreversible in reverse direction
IrR = model.ub<=0; % reactions irreversible in reverse direction
temp = model.ub(IrR);
model.S(:,IrR) = -model.S(:,IrR);
model.ub(IrR) = -1*model.lb(IrR);
model.lb(IrR) = -1*temp;
modModel = model;
% normalising the bounds to max of 1000
temp1 = max(abs(model.lb));
temp2 = max(abs(model.ub));
model.lb(abs(model.lb)==temp1) = sign(model.lb(abs(model.lb)==temp1))*1000;
model.ub(abs(model.ub)==temp2) = sign(model.ub(abs(model.ub)==temp2))*1000;


if getDirection
    % the model should not have any reaction with flux only in reverse
    % direction
    [workModel, ~, rev2irrev] = convertToIrreversible(model);
else
    workModel = model;
    rev2irrev = num2cell(1:numel(model.rxns))';
end
[~,n] = size(workModel.S);
% isolate the model's TICs from the working model
model_tic = workModel;
% converting all the positive lower bounds to zero lower bounds
model_tic.lb(model_tic.lb>0)=0;



ind=findExcRxns(model_tic);
% blocking all the exchange reactions
model_tic.lb(ind) = 0; model_tic.ub(ind) = 0;

if exist('spectraCC','file')
    tic_idx = spectraCC(model_tic,tol); 
else
    tic_idx = fastcc(model_tic,tol,0);
end

if isempty(tic_idx)
    Nsnp = [];
else
    model_tic = removeRxns(model_tic,model_tic.rxns(setdiff(1:numel(model_tic.rxns),tic_idx)));
    Nsnp = fast_snp(model_tic.S, model_tic.lb, model_tic.ub, tol);
end

irrev2rev = zeros(n, 1);
for i = 1:length(rev2irrev)
    irrev2rev(rev2irrev{i}) = i;
end

core = 1:n; 
runtime=120;
while 1
    [IDS,stat] = findConsistentIDS(workModel,core,Nsnp,tic_idx,tol,runtime);

    if isempty(IDS)
        if stat==3
            runtime = runtime+60;
            continue
        else
            break
        end
    else
        core = setdiff(core,IDS);
    end
end

consIrr = double(~ismember([1:n],core)); % consistent irreversible reactions
if getDirection
    a = cellfun(@(x)getConsDir(consIrr,x),rev2irrev,'UniformOutput',0);
else
    a = cellfun(@(x)getConsDirBlocked(consIrr,x),rev2irrev,'UniformOutput',0);
end
end

function consDir = getConsDir(consIrr,x)
temp = consIrr(x);
if sum(temp)==2
    consDir = 'Reversible'; % consistent in both the direction
elseif sum(temp)==0
    consDir = 'Blocked'; % blocked reaction
elseif temp(1)==1
    consDir = 'Forward'; % consistent only in forward direction
elseif temp(2)==1
    consDir = 'Reverse'; % consistent only in reverse direction
end
end

function consDir = getConsDirBlocked(consIrr,x)
temp = consIrr(x);
if sum(temp)>0
    consDir = 'Unblocked';
else
    consDir = 'Blocked';
end
end


