
classdef Parameters
    properties
        nk = 101
        kcurve = 7
        kmin = 1e-8
        kmax = 100
        nK = 8
        Kcurve = 7
        Kmin = 10
        Kmax = 80
        
        delta = 0.025
        beta = 0.99
        alpha = 0.36
        
        l = [0,1]
        nl
        z = [0.99,1.01]
        nz
        
        pimat
        
        dims
        
        maxiters = 1e4
        tol = 1e-8
        del_update = 0.4
        
        L = 0.96
    end  
    
    methods
        function obj = Parameters()
            obj.nl = numel(obj.l);
            obj.nz = numel(obj.z);
            obj.dims = [obj.nk,obj.nK,obj.nl,obj.nz];
            
            obj.z = shiftdim(obj.z,-2);
            
            obj.pimat = [...
                0.8507, 0.1229, 7/12, 3/32
                0.1159, 0.8361, 1/32, 7/20
                0.0243, 0.0021, 7/24, 1/32
                0.0091, 0.0389, 3/32, 21/40
                ];
        end
    end
end
