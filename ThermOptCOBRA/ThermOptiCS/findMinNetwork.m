function IDS = findMinNetwork(model,core,P,Nsnp,tic_idx,tol,x0,timeLimit)
% Identifies the thermodynamically consistent minimal network
% 
% USAGE: 
%   IDS = findMinNetwork(model,core,P,Nsnp,tic_idx,tol,x0,timeLimit)
%
% INPUTS:
%     model:     COBRA model structure
%     core:      Reaction IDs for which thermodynamic consistency has been verified
%     P:         Penalty set (reactions to minimize)
%     Nsnp:      Sparse null space matrix from Fast-SNP
%     tic_idx:   Indices of TIC reactions in the model
%     tol:       Minimum positive number defined as non-zero
%     x0:        Initial feasible point for the MILP problem
%     timeLimit: timeLimit for solving the minNetwork MILP problem
%
% OUTPUTS:
%     IDS:       Indices of reactions in the minimal network

if ~exist('timeLimit','var') || isempty(timeLimit)
    timeLimit = 120;
end

[m,n] = size(model.S);

if isempty(Nsnp)
    K = 0;
    n_tic = 0;
else
    K = size(Nsnp, 2);
    n_tic = length(tic_idx);
end

% Identify reversible reactions (lb < 0)
rev_rxns = find(model.lb < 0);
n_rev = length(rev_rxns);

% Total variables: v (n), z+ (n), z- (n_rev), G (n_tic)
n_vars = 2*n + n_rev + n_tic;

% Map z- variables to full reaction indices
M_rev = sparse(rev_rxns, 1:n_rev, ones(n_rev, 1), n, n_rev);

% Reverse lookup for z- variables
rev_idx = zeros(n, 1);
rev_idx(rev_rxns) = 1:n_rev;

% Minimum flux magnitudes for all reactions
min_fwd = max(tol, model.lb);
min_rev = tol * ones(n, 1);

% Objective: Minimize active reactions in penalty set P
c_z_plus = zeros(n,1); c_z_plus(P) = 1;
c_z_minus = zeros(n_rev,1); 
p_rev = ismember(rev_rxns, P);
c_z_minus(p_rev) = 1;
f = [zeros(n,1); c_z_plus; c_z_minus; zeros(n_tic,1)];

% Equalities
% 1. Steady state: S*v = 0
Aeq1 = [model.S, sparse(m, n + n_rev + n_tic)];
beq1 = zeros(m,1);
csenseeq1 = repmat('E', m, 1);

% 2. Loopless: Nsnp^T * G = 0
if K > 0
    Aeq2 = [sparse(K, 2*n + n_rev), Nsnp'];
    beq2 = zeros(K, 1);
    csenseeq2 = repmat('E', K, 1);
else
    Aeq2 = [];
    beq2 = [];
    csenseeq2 = [];
end

Aeq = [Aeq1; Aeq2];
beq = [beq1; beq2];
csenseeq = [csenseeq1; csenseeq2];

% Inequalities
temp_v = speye(n);

% 1. v_i - ub_i*z_i^+ + min_rev_i*z_i^- <= 0
Aineq1 = [temp_v, -1*sparse(diag(model.ub)), sparse(diag(min_rev))*M_rev, sparse(n, n_tic)];
bineq1 = zeros(n,1);
csenseineq1 = repmat('L', n, 1);

% 2. v_i - min_fwd_i*z_i^+ - lb_i*z_i^- >= 0
Aineq2 = [temp_v, -1*sparse(diag(min_fwd)), -1*sparse(diag(model.lb))*M_rev, sparse(n, n_tic)];
bineq2 = zeros(n,1);
csenseineq2 = repmat('G', n, 1);

% 3. z_i^+ + z_i^- <= 1 (only for reversible reactions)
if n_rev > 0
    Aineq3 = [sparse(n_rev, n), M_rev', speye(n_rev), sparse(n_rev, n_tic)];
    bineq3 = ones(n_rev, 1);
    csenseineq3 = repmat('L', n_rev, 1);
else
    Aineq3 = [];
    bineq3 = [];
    csenseineq3 = [];
end

% 4. Fast-SNP Loopless Eq 5
if n_tic > 0
    % G_j + 1001 z_k^+ <= 1000
    Aineq4 = sparse(n_tic, n_vars);
    for j = 1:n_tic
        k = tic_idx(j);
        Aineq4(j, n + k) = 1001; % z_k^+
        Aineq4(j, 2*n + n_rev + j) = 1; % G_j
    end
    bineq4 = 1000 * ones(n_tic, 1);
    csenseineq4 = repmat('L', n_tic, 1);
    
    % G_j - 1001 z_k^- >= -1000  (only if k is reversible)
    tic_rev_flags = ismember(tic_idx, rev_rxns);
    num_tic_rev = sum(tic_rev_flags);
    if num_tic_rev > 0
        Aineq5 = sparse(num_tic_rev, n_vars);
        row = 1;
        for j = 1:n_tic
            k = tic_idx(j);
            if rev_idx(k) > 0
                Aineq5(row, 2*n + rev_idx(k)) = -1001; % z_k^-
                Aineq5(row, 2*n + n_rev + j) = 1; % G_j
                row = row + 1;
            end
        end
        bineq5 = -1000 * ones(num_tic_rev, 1);
        csenseineq5 = repmat('G', num_tic_rev, 1);
    else
        Aineq5 = [];
        bineq5 = [];
        csenseineq5 = [];
    end
else
    Aineq4 = [];
    bineq4 = [];
    csenseineq4 = [];
    Aineq5 = [];
    bineq5 = [];
    csenseineq5 = [];
end

% 5. Force core reactions to be active: z_i^+ + z_i^- >= 1
n_core = length(core);
if n_core > 0
    Aineq6 = sparse(n_core, n_vars);
    for r = 1:n_core
        i = core(r);
        Aineq6(r, n + i) = -1; % -z_i^+
        if rev_idx(i) > 0
            Aineq6(r, 2*n + rev_idx(i)) = -1; % -z_i^-
        end
    end
    bineq6 = -1 * ones(n_core, 1);
    csenseineq6 = repmat('L', n_core, 1); % -z_i^+ - z_i^- <= -1
else
    Aineq6 = [];
    bineq6 = [];
    csenseineq6 = [];
end

% Bounds
lb = [model.lb; zeros(n, 1); zeros(n_rev, 1); -1000 * ones(n_tic, 1)];
ub = [model.ub; ones(n, 1); ones(n_rev, 1); 1000 * ones(n_tic, 1)];

% Set up MILP problem
MILPproblem.A=[Aeq; Aineq1; Aineq2; Aineq3; Aineq4; Aineq5; Aineq6];
MILPproblem.b=[beq; bineq1; bineq2; bineq3; bineq4; bineq5; bineq6];
MILPproblem.lb=lb;
MILPproblem.ub=ub;
MILPproblem.c=f;
MILPproblem.osense=1; % minimise
MILPproblem.vartype = [repmat('C',n,1); repmat('B',n,1); repmat('B',n_rev,1); repmat('C',n_tic,1)];
MILPproblem.csense = [csenseeq; csenseineq1; csenseineq2; csenseineq3; csenseineq4; csenseineq5; csenseineq6];

if ~exist('x0', 'var') || isempty(x0)
    MILPproblem.x0 = [];
else
    MILPproblem.x0 = x0;
end

changeCobraSolverParams('MILP', 'feasTol', 1e-9);
changeCobraSolverParams('MILP', 'intTol', 1e-9);

solution = solveCobraMILP(MILPproblem, 'timeLimit', timeLimit);

if solution.stat == 1 || solution.stat == 3
    x = solution.full;
    z_plus = x(n+1 : 2*n);
    z_minus_full = zeros(n, 1);
    v = x(1:n);
    if n_rev > 0
        z_minus_full(rev_rxns) = x(2*n+1 : 2*n+n_rev);
    end
    active_mask = (z_plus >= 0.5) | (z_minus_full >= 0.5) | (abs(v)>=tol*1e-3);
    IDS = find(active_mask);
else
    IDS = [];
end
end
