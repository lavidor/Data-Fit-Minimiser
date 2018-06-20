function [a , chi2Min] = Chi2FindMin(x,y, e, f, a, aLim,aLimMin)
    k = length(a);
    if(length(aLimMin)== 1)
        aLimMin = aLimMin*ones(1,k);
    end
    if(length(aLim)< k)
        aLim =a/5;
    end
    aN = 3*ones(1,k);
    cont = 1;
    while(cont > 0)
        cont = 0;
        aLim = aLim/2;
        [chiGrid , aVectors , chi2Min , minIndex] = Chi2Grid(x,y, e, f, a, aLim, aN);
        i = 0;
        while (i<k)
            i = i+1;
            a(i) = aVectors(i , minIndex{i});
            if(aLim(i)>aLimMin(i))
                cont = 1;
            else
                aN(i) = 1;
            end
        end
        if(all(chiGrid == chi2Min))
            cont = 0;
        end
    end
end