
clear

% Read parameters file
p = Parameters();

% Create grids for household assets and aggregate capital
grids = struct();
grids.k = create_grid(p.kmin,p.kmax,p.nk,p.kcurve);
grids.K = linspace(p.Kmin,p.Kmax,p.nK)';

% Guess policy function k'(k,K,s,z)
kpol0 = 0.9 * grids.k .* ones([1,p.nK,p.nl,p.nz]);

% Guess capital law of motion
lom0_K = struct();
lom0_K.alpha = reshape([0,0],[1,1,1,2]);
lom0_K.beta = reshape([1,1],[1,1,1,2]);

% Solve
[kpol, lom] = iterate_lom(p,grids,kpol0,lom0_K);
[~,K_t] = simulate(p,grids,kpol);


% Iterate over law of motion
function [kpol, lom_K] = iterate_lom(p,grids,kpol0,lom0_K)
    lom_K = lom0_K;
    
    for i=1:p.maxiters
        % Solve for policy function
        kpol = solve_policy(p,grids,kpol0,lom_K);

        % Simulate
        b1 = simulate(p,grids,kpol);

        % Check difference
        b0 = [lom_K.alpha(:)';lom_K.beta(:)'];
        bnorm = norm(b1(:)-b0(:),Inf);
        delta = 0.5;
        lom_K.alpha(1) = delta*b0(1,1) + (1-delta)*b1(1,1);
        lom_K.alpha(2) = delta*b0(1,2) + (1-delta)*b1(1,2);
        lom_K.beta(1) = delta*b0(2,1) + (1-delta)*b1(2,1);
        lom_K.beta(2) = delta*b0(2,2) + (1-delta)*b1(2,2);
        
        if bnorm < 1e-5
            disp('Converged')
            return
        elseif mod(i,1)==0
        	fprintf('LOM iteration = %d, norm=%f\n',i,bnorm)
        end
    end
    error('No convergence')
end

function kpol = solve_policy(p,grids,kpol0,lom_K)
    kpol = kpol0;
    for i=1:p.maxiters
        
        kpol_update = update_policy(p,grids,kpol,lom_K);
        knorm = norm(kpol(:)-kpol_update(:),Inf);
        
        % Update decision rule
        kpol = p.del_update * kpol_update...
            + (1-p.del_update) * kpol;
        
        if knorm < p.tol
            fprintf('\tConverged\n')
            return
        elseif mod(i,100)==0
        	fprintf('\tPolicy iteration = %d, norm=%f\n',i,knorm)
        end
    end
    error('No convergence')
end

function kpol_euler = update_policy(p,grids,kpol,lom_K)
    % Updates the policy function
    
    % Conjectured values for K'(K,z), repeated onto k and l dimensions.
    % dim(Kpvec) = (nk*nK*nl,nz)
    Kp = capital_conjecture(lom_K,grids.K');
    Kpvec = repmat(Kp,[p.nk,1,p.nl,1]);
    Kpvec = reshape(Kpvec,[],p.nz);
    
    % Prices r(K,z) and w(K,z)
    % dim(r) = dim(w) = (1,nk,1,nz)
    [r, w] = compute_prices(p, grids.K');
    
    % Prices r(K',z') and w(K',z')
    [rp, wp] = compute_prices(p, Kp);
    
    % Repeat r' and w' onto (k,l,z,l') dimensions and reshape.
    % Rows correspond with current state and columns correspond with
    % (l',z').
    % dim(rp) = dim(wp) = (nk*nK*nl*nz,nl*nz)
    rp = reshape(rp,[1,p.nK,1,1,1,p.nz]);
    rp = repmat(rp,[p.nk,1,p.nl,p.nz,p.nl,1]);
    rp = reshape(rp,[],p.nl*p.nz);
    wp = reshape(wp,[1,p.nK,1,1,1,p.nz]);
    wp = repmat(wp,[p.nk,1,p.nl,p.nz,p.nl,1]);
    wp = reshape(wp,[],p.nl*p.nz);
    
    % Next periods labor supply, dimension (nl,nz)
    lp = repmat(p.l(:),[1,p.nz]);

    % Interpolate k'(k'(k,K,l,z),K',l',z')
    kpp = zeros([p.nk*p.nK*p.nl*p.nz,p.nl,p.nz]);
    for il=1:p.nl % l'
        for iz=1:p.nz % z'
            kp_lz = kpol(:,:,il,iz);
            kinterp = griddedInterpolant(...
                {grids.k,grids.K},kp_lz,'spline');
            kpp(:,il,iz) = kinterp(kpol(:),Kpvec(:));
        end
    end
    
    % Construct c' vector
    kpp = reshape(kpp,[],p.nl*p.nz);
    cp = (1+rp-p.delta).*kpol(:) + wp.*lp(:)' - kpp;
    
    % Expected marginal utility
    A = kron(p.pimat,ones(p.nk*p.nK,1));
    eMU = sum(((1+rp-p.delta)./cp) .* A, 2);
    
    % k' decision rule on LHS of Euler
    eMU = reshape(eMU,p.dims);
    kpol_euler = (1+r-p.delta).*grids.k + w.*shiftdim(lp,-2) - 1./(p.beta*eMU);
    kpol_euler = max(kpol_euler,p.kmin);
    kpol_euler = min(kpol_euler,p.kmax);
end

function Kp = capital_conjecture(lom_K, K)
    logKp = lom_K.alpha + lom_K.beta .* log(K);
    Kp = exp(logKp);
end

function [r, w] = compute_prices(p, K)
    klratio = K ./ reshape(p.L,[1,1,1,p.nz]);

    % r(K,z), dimension (1,nK,1,nz)
    r = p.z * p.alpha .* klratio .^ (p.alpha-1);
    r = reshape(r, [1,p.nK,1,p.nz]);
    
    % w(K,z), dimension (1,nK,1,nz)
    w = p.z * (1-p.alpha) .* klratio .^(p.alpha);
    w = reshape(w, [1,p.nK,1,p.nz]);
end

function [b,K_t] = simulate(p,grids,kpol)

    % Create policy interpolants
    kinterp = cell(p.nl,p.nz);
    for il=1:p.nl % l'
        for iz=1:p.nz % z'
            kp_lz = kpol(:,:,il,iz);
            kinterp{il,iz} = griddedInterpolant(...
                {grids.k,grids.K},kp_lz,'spline');
        end
    end
    
    % Time periods to simulate
    T = p.sim_tburn + p.sim_T;
    
    rng(9520);
    
    % Aggregate productivity: iz=1 bad, iz=2 good
    zrand = rand(T,1);
    iz0 = 1;
    
    % Initial capital distribution
    k = linspace(p.kmin,p.kmax,p.sim_nHH)';
    
    % Initial employment status
    employed = rand(p.sim_nHH,1) < p.L(iz0);
    l = (~employed) * 1 + (employed) * 2;
    
    % Employment transitions conditional on agg transition
	pi_l_trans = cell(2,2);
%     pi_l_trans{1,1} = p.pimat([1,3],[1,3]);
%     pi_l_trans{1,2} = p.pimat([1,3],[2,4]);
%     pi_l_trans{2,1} = p.pimat([2,4],[1,3]);
%     pi_l_trans{2,2} = p.pimat([2,4],[2,4]);

    pi_l_trans{1,1} = p.pimat([1,2],[1,2]);
    pi_l_trans{1,2} = p.pimat([1,2],[3,4]);
    pi_l_trans{2,1} = p.pimat([3,4],[1,2]);
    pi_l_trans{2,2} = p.pimat([3,4],[3,4]);
    
    for i=1:2
        for j=1:2
            pi_l_trans{i,j} = pi_l_trans{i,j} ./ sum(pi_l_trans{i,j},2);
            pi_l_trans{i,j} = cumsum(pi_l_trans{i,j},2);
        end
    end
    
    K_t = zeros(T,1);
    z_t = zeros(T,1);
    for t=1:T
        K_t(t) = mean(k);
        z_t(t) = iz0;

        % Update next aggregate productivity
        iz1 = 1 + (zrand(t) < p.pi_z(1,iz0));
                
        % Draw idiosyncratic shocks -- do I need to switch the order?
        lrand = rand(p.sim_nHH,1);
        [~,l] = max(lrand<=pi_l_trans{iz0,iz1}(l,:),[],2);
        
        % Interpolate
        Kvec = K_t(t) * ones(p.sim_nHH,1);
        
        employed = (l==1);
        k(~employed) = kinterp{1,iz0}(k(~employed),Kvec(~employed));
        k(employed) = kinterp{2,iz0}(k(employed),Kvec(employed));
        
        % Update current productivity
        iz0 = iz1;
    end
    
    K_t = K_t(p.sim_tburn+1:T);
    z_t = z_t(p.sim_tburn+1:T);

    logKt = log(K_t);
    
    b = zeros(2,2);
    for iz=[1,2]
        condition = (z_t == iz);
        icond = find(condition);
        icond = icond(1:end-1);

        nc = numel(icond);
        X = [ones(nc,1), logKt(icond)];
        Y = logKt(icond+1);
        b(:,iz) = (X'*X) \ (X'*Y);
    end
end
