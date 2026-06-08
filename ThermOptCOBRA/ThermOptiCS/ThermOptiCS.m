function [Model,bCoreRxns] = ThermOptiCS(model,core,tol,doTOCC,timeLimit)
% Thermodynamically feasible context-specific model building.
%
% USAGE: 
%   [Model,bCoreRxns] = ThermOptiCS(model,core,tol,doTOCC,timeLimit)
%
% INPUTS:
%     model:     COBRA model structure for which act as the GSMM
%     core:      Reaction IDs which are defined to be core (These
%                reactions will be present in the final model)
%     tol:       Tolerance value (User defined non-zero value).
%
% OPTIONAL INPUTS:
%     doTOCC:    Boolean value.
%                  1: ThermOptCC will be run on the model and blocked
%                     reactions will be removed (Default)
%                  0: ThermOptCC will not be run on the model (This
%                     requires the input model to be thermodynamically-flux
%                     consistent)
%     timeLimit: Upper time limit in seconds to solve the minReac
%                optimization problem (Default: 120s)
%
% OUTPUTS:
%     Model:     Context-specific model
%     bCoreRxns: Core reactions that are thermodynamically blocked (or infeasible)
%
% .. Author:
%       - Pavan Kumar S, BioSystems Engineering and control (BiSECt) lab, IIT Madras


% checking for reactions irreversible in reverse direction
IrR = model.ub<=0; % reactions irreversible in reverse direction
temp = model.ub(IrR);
model.S(:,IrR) = -model.S(:,IrR);
model.ub(IrR) = -1*model.lb(IrR);
model.lb(IrR) = -1*temp;

% normalising the bounds to max of 1000
temp1 = max(abs(model.lb));
temp2 = max(abs(model.ub));
model.lb(abs(model.lb)==temp1) = sign(model.lb(abs(model.lb)==temp1))*1000;
model.ub(abs(model.ub)==temp2) = sign(model.ub(abs(model.ub)==temp2))*1000;

if ~exist('doTOCC','var')||isempty(doTOCC)
    doTOCC=1;
end

if doTOCC
    [a, ~] = ThermOptCC(model,tol,false);
else
    a = repmat({'Unblocked'},numel(model.rxns),1);
end

if ~exist('timeLimit','var')||isempty(timeLimit)
    timeLimit=120;
end

bRxns = find(ismember(a,'Blocked')); % thermodynamically blocked reactions
bCoreRxns = intersect(bRxns,core); % getting the blocked core reactions
if ~isempty(bCoreRxns)
    warning('Few reactions in the core reaction list are thermodynamically blocked')
    fprintf('\n Thermodynamically blocked reactions will be removed from core reaction list\n')
    core =  setdiff(core,bRxns);
end

% Compute Nsnp and tic_idx directly on the reversible model
model_tic = model;
model_tic.lb(model_tic.lb>0)=0;
ind=findExcRxns(model_tic);
model_tic.lb(ind) = 0; model_tic.ub(ind) = 0;

if exist('spectraCC','file')
    tic_idx = spectraCC(model_tic,tol); 
else
    tic_idx = fastcc(model_tic,tol,0);
end

if isempty(tic_idx)
    Nsnp = [];
else
    model_tic_reduced = removeRxns(model_tic,model_tic.rxns(setdiff(1:numel(model_tic.rxns),tic_idx)));
    Nsnp = fast_snp(model_tic_reduced.S, model_tic_reduced.lb, model_tic_reduced.ub, tol);
end

minNetwork = [];
P = setdiff(1:numel(model.rxns), core); % penalty set
while ~isempty(core)
    [IDScc,x0] = findConsistentIDS_TOCS(model,core,Nsnp,tic_idx,tol,timeLimit);
    
    if isempty(IDScc)
        error('One or more core reactions cannot carry non-zero flux (MILP infeasible).');
    end
    
    IDS = findMinNetwork(model,intersect(core,IDScc),P,Nsnp,tic_idx,tol,x0,timeLimit);
    minNetwork = [minNetwork(:); IDS(:)];
    
    core = setdiff(core,IDS);
    P = setdiff(P,IDS);
end
minNetwork = unique(minNetwork);
Model = removeRxns(model,model.rxns(setdiff(1:numel(model.rxns),minNetwork)));

end

