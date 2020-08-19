% Last updated: 1 Oct 2019

function Dfull = PA_dis(submat, P, pct)

Ns = length(submat);
D = zeros(Ns);

for ll = 1:(Ns-1)
    for rr = (ll+1):Ns
        
        F = norml2(submat{ll}, 1);
        G = norml2(submat{rr}, 1);
        
        [theta diff dmax] = PA_calc(F,G,pct);
        D(ll,rr) = (sum(1.-cos(theta))+diff)/dmax;
    end
end

Dfull = D + D';

end