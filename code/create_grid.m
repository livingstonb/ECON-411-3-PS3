function kgrid = create_grid(lb,ub,n,pow)
    kgrid = linspace(0,1,n)';
    kgrid = lb + (ub-lb) * kgrid.^pow;
end