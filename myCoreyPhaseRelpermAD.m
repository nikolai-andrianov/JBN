function fn = myCoreyPhaseRelpermAD(n, sr, kwm, sr_tot)
    fn = @(s, varargin) myCoreyRelperm(s, n, sr, kwm, sr_tot);
end

function kr = myCoreyRelperm(s, n, sr, kwm, sr_tot)
    den = 1 - sr_tot;
    sat = ((s - sr)./den);
    sat = max(min(sat, 1), 0);
    
    kr = kwm*sat.^n;
    
    if isa(sat, 'ADI')
       if any(~isfinite(kr.val)) 
          a = 1; 
       end
    end
end
