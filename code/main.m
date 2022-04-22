
clear

p = Parameters();
kgrid = create_grid(p.kmin,p.kmax,p.nk,p.kcurve);
aggKgrid = create_grid(p.Kmin,p.Kmax,p.nK,p.Kcurve);

% Guess policy function k'(k,K,s,z)
kpol0 = repmat(0.9*kgrid, [1,p.nK,p.nl,p.nz]);

% Guess capital law of motion
lom_K = struct();
lom_K.alpha = shiftdim([0;0],-3);
lom_K.beta = shiftdim([1;1],-3);

% Solve for policy function
kpol = iterate(p,kgrid,aggKgrid,kpol0,lom_K);

disp('done')

function kpol = iterate(p,kgrid,aggKgrid,kpol0,lom_K)
    kpol = kpol0;
    for i=1:p.maxiters
        
        kpol_update = solve_policy(p,kgrid,aggKgrid,kpol,lom_K);
        knorm = norm(kpol(:)-kpol_update(:),Inf);
        kpol = kpol_update;
        
        if knorm < p.tol
            disp('Converged')
            return
        else
            fprintf('Iteration = %d, norm=%f\n',i,knorm)
        end
    end
    error('No convergence')
end

function kpol_update = solve_policy(p,kgrid,aggKgrid,kpol,lom_K)
    % Solves the policy function given the aggregate states and the
    % conjecture for the law of motion.
    
    % Conjectured K', agg.Kp
    Kp = capital_conjecture(lom_K,aggKgrid');
    Kpvec = repmat(Kp,[p.nk,1,p.nl,1]);
    Kpvec = reshape(Kpvec,[],p.nz);
    
    % Prices
    [r, w] = compute_prices(p, aggKgrid');
    [rp, wp] = compute_prices(p, Kp);
    rp = reshape(rp,[1,p.nK,1,1,1,p.nz]);
    rp = repmat(rp,[p.nk,1,p.nl,p.nz,p.nl,1]);
    rp = reshape(rp,[],p.nl*p.nz);
    wp = reshape(wp,[1,p.nK,1,1,1,p.nz]);
    wp = repmat(wp,[p.nk,1,p.nl,p.nz,p.nl,1]);
    wp = reshape(wp,[],p.nl*p.nz);
    
    % Labor supply, dimension (nl,nz)
    l = reshape(p.l, [p.nl,1]);
    l = repmat(l,[1,p.nz]);
    lp = l(:)';

    % Interpolate k'(k'(k,K,l,z),K',l',z')
    kpp = zeros([p.nk*p.nK*p.nl*p.nz,p.nl,p.nz]);
    for il=1:p.nl % l'
        for iz=1:p.nz % z'
            kp_lz = kpol(:,:,il,iz);
            kinterp = griddedInterpolant(...
                {kgrid,aggKgrid},kp_lz,'spline');
            kpp(:,il,iz) = kinterp(kpol(:),Kpvec(:));
        end
    end
    
    % Construct c' vector
    kpp = reshape(kpp,[],p.nl*p.nz);
    cp = (1+rp+p.delta).*kpol(:) + wp.*lp - kpp;
    
    % Expected marginal utility
    A = kron(p.pimat,ones(p.nk*p.nK,1));
    eMU = sum((1./cp) .* A, 2);
    
    % k' decision rule on LHS of Euler
    eMU = reshape(eMU,p.dims);
    kpol_euler = (1+r+p.delta).*kgrid + w.*shiftdim(l,-2) - 1./(p.beta*eMU);
    kpol_euler = max(kpol_euler,p.kmin);
    kpol_euler = min(kpol_euler,p.kmax);
    
    
    % Update decision rule
    kpol_update = p.del_update * kpol_euler...
        + (1-p.del_update) * kpol;
end

function Kp = capital_conjecture(lom_K, K)
    logKp = lom_K.alpha + lom_K.beta .* log(K);
    Kp = exp(logKp);
end

function [r, w] = compute_prices(p, K)
    klratio = K ./ p.L;

    % r(K,z), dimension (1,nK,1,nz)
    r = p.z * p.alpha .* klratio .^ (p.alpha-1);
    r = reshape(r, [1,p.nK,1,p.nz]);
    
    % w(K,z), dimension (1,nK,1,nz)
    w = p.z * (1-p.alpha) .* klratio .^(p.alpha);
    w = reshape(w, [1,p.nK,1,p.nz]);
end
