
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
lom0_K = [0,0;1,1];

% Solve
[kpol,lom,results] = iterate_lom(p,grids,kpol0,lom0_K);
[~,K_t,iz_t] = simulate(p,grids,kpol);

% Compute conjectured path for capital stock
T = numel(K_t);
Kapprox = zeros(numel(K_t),1);
Kapprox(1) = K_t(1);
for it = 1:T-1
     Kp = reshape(capital_conjecture(lom, Kapprox(it)),[2,1]);
     Kapprox(it+1) = Kp(iz_t(it));
end

% Capital time series
plot(1:numel(K_t),K_t)
hold on
plot(1:numel(K_t),Kapprox)
hold off
ylabel("K")
xlabel("t")
legend("K_t, simulated","K_t, approximated from LoM")
set(gcf,'color','w');
savefig(gcf,'k_t.png')

% Print results
fprintf('Bad state: alpha=%f, beta=%f\n', lom(1,1), lom(2,1));
fprintf('\tR-squared=%f\n', results{1}.r2);
fprintf('Good state: alpha=%f, beta=%f\n', lom(1,2), lom(2,2));
fprintf('\tR-squared=%f\n', results{2}.r2);

% Iterate over law of motion
function [kpol, lom_K, sim_results] = iterate_lom(p,grids,kpol0,lom0_K)
    lom_K = lom0_K;
    
    for i=1:p.maxiters
        % Solve for policy function
        kpol = solve_policy(p,grids,kpol0,lom_K);

        % Simulate
        sim_results = simulate(p,grids,kpol);
        lom1_K = [sim_results{1}.b,sim_results{2}.b];

        % Check difference
        bnorm = norm(lom1_K(:)-lom_K(:),Inf);

        if bnorm < 1e-5
            disp('Converged')
            return
        elseif mod(i,1)==0
        	fprintf('LOM iteration = %d, norm=%f\n',i,bnorm)
        end
        
        % Update
        lom_delta = 0.3;
        lom_K = lom_delta*lom1_K + (1-lom_delta)*lom_K;
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
        elseif mod(i,250)==0
        	fprintf('\tPolicy iteration = %d, norm=%f\n',i,knorm)
        end
    end
    error('No convergence')
end

function kpol_euler = update_policy(p,grids,kpol,lom_K)
    % Updates the policy function
    
    % Conjectured values for K'(K,z), repeated onto k and l dimensions.
    % dim(Kpvec) = (nk,nK,nl,nz)
    Kp = capital_conjecture(lom_K,grids.K');
    Kpvec = repmat(Kp,[p.nk,1,p.nl,1]);
    
    % Prices r(K,z) and w(K,z)
    % dim(r) = dim(w) = (1,nK,1,nz)
    [r, w] = compute_prices(p, grids.K');
    
    % Prices r(K',z') and w(K',z')
    % Dimensions (1,nK,1,nz)
    [rp, wp] = compute_prices(p, Kp);
    
    % Repeat along dimensions of length 1
    rp = reshape(rp,[1,p.nK,1,1,1,p.nz]);
    rp = repmat(rp,[p.nk,1,p.nl,p.nz,p.nl,1]);
    wp = reshape(wp,[1,p.nK,1,1,1,p.nz]);
    wp = repmat(wp,[p.nk,1,p.nl,p.nz,p.nl,1]);
  
    % Reshape so columns are next period's shocks
    rp = reshape(rp,[],p.nl*p.nz);
    wp = reshape(wp,[],p.nl*p.nz);
    
    % Next periods labor supply, dimension (nl,nz)
    lp = repmat(p.l(:),[1,p.nz]);

    % Interpolate k'(k'(k,K,l,z),K',l',z')
    kpp = zeros([p.nk*p.nK*p.nl*p.nz,p.nl,p.nz]);
    for il=1:p.nl % l'
        for iz=1:p.nz % z'
            kinterp = griddedInterpolant(...
                {grids.k,grids.K},kpol(:,:,il,iz),'spline');
            kpp(:,il,iz) = kinterp(kpol(:),Kpvec(:));
        end
    end
    
    % Construct c' vector
    kpp = reshape(kpp,[],p.nl*p.nz);
    cp = (1+rp-p.delta).*kpol(:) + wp.*lp(:)' - kpp;
    cp = max(cp,1e-8);
    
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
    logKp =  [ones(numel(K),1), log(K(:))] * lom_K;
    Kp = reshape(exp(logKp),[1,numel(K),1,2]);
end

function [r, w] = compute_prices(p, K)
    klratio = K ./ reshape(p.L,[1,1,1,p.nz]);

    % r(K,z), dimension (1,nK,1,nz)
    r = reshape(p.z,[1,1,1,p.nz]) .* p.alpha .* klratio .^ (p.alpha-1);
    r = reshape(r, [1,p.nK,1,p.nz]);
    
    % w(K,z), dimension (1,nK,1,nz)
    w = reshape(p.z,[1,1,1,p.nz]) .* (1-p.alpha) .* klratio .^(-p.alpha);
    w = reshape(w, [1,p.nK,1,p.nz]);
end

function [regs,K_t,iz_t] = simulate(p,grids,kpol)
    rng(1324);

    % Create policy interpolants
    kinterp = cell(p.nl,p.nz);
    for il=1:p.nl % l'
        for iz=1:p.nz % z'
            kinterp{il,iz} = griddedInterpolant(...
                {grids.k,grids.K},kpol(:,:,il,iz),'spline');
        end
    end
    
    % Time periods to simulate
    T = p.sim_tburn + p.sim_T;
    
    % Aggregate productivity shocks: iz=1 bad, iz=2 good
    zrand = rand(T,1);
    iz0 = 1;
    
    % Idiosyncratic shocks
    lrand = rand(p.sim_nHH,T);
    
    % Initial capital distribution
    k = ones(p.sim_nHH,1) * 35;
    
    % Initial employment status
    employed = rand(p.sim_nHH,1) < p.L(iz0);
    i_employed = (~employed) * 1 + (employed) * 2;
    
    % Employment transitions conditional on agg transition
    indices = @(k) int8((k==1)*[1,2]+(k==2)*[3,4]);
	pi_l_trans = cell(2,2);
	for i=1:2
        for j=1:2
            pi_l_trans{i,j} = p.pimat(indices(i),indices(j));
            pi_l_trans{i,j} = pi_l_trans{i,j} ./ sum(pi_l_trans{i,j},2);
            pi_l_trans{i,j} = cumsum(pi_l_trans{i,j},2);
        end
    end
    
    K_t = zeros(T,1);
    iz_t = zeros(T,1);
    for t=1:T
        K_t(t) = mean(k);
        iz_t(t) = iz0;

        % Interpolate
        Kvec = K_t(t) * ones(p.sim_nHH,1);
        k(~employed) = kinterp{1,iz0}(k(~employed),Kvec(~employed));
        k(employed) = kinterp{2,iz0}(k(employed),Kvec(employed));
        k = max(k,p.kmin);
        k = min(k,p.kmax);
        
        % Update next aggregate productivity
        iz1 = 1 + (zrand(t) > p.pi_z(iz0,1));
        
        % Draw idiosyncratic shocks conditional on aggregate transition
        [~,i_employed] = max(lrand(:,t)<=pi_l_trans{iz0,iz1}(i_employed,:),[],2);
        employed = (i_employed==1);
        
        % Update current productivity
        iz0 = iz1;
    end
    
    K_t = K_t(p.sim_tburn+1:T);
    iz_t = iz_t(p.sim_tburn+1:T);

    logKt = log(K_t);
    regs = cell(p.nz,1);
    for iz=1:p.nz
        % Regression conditional on aggregate state at time t
        icond = find(iz_t == iz);
        icond = icond(1:end-1);
        regs{iz} = regress(logKt(icond),logKt(icond+1));
    end
end

function results = regress(covariates,Y)
        X = [ones(size(covariates,1),1), covariates];
        results = struct();
        
        % Estimate by least squares
        results.b = (X'*X) \ (X'*Y);
        
        % Compute r2 and sigma
        u = Y - X * results.b;
        rss = u' * u;
        results.r2 = 1 - rss/((Y-mean(Y))'*(Y-mean(Y)));
        results.sigma = sqrt(results.r2 / numel(Y)-size(X,2));
end
