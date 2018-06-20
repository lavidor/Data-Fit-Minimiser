function p = Chi2CalcDirect(y,eErr,yTheory)
    if (eErr == 0)
        eErr = ones(1,length(y));
    end
    n = min([length(yTheory(:)),length(y(:)),length(eErr(:))]);
    if(n > 0)
        p = sum(((y(1:n)-yTheory(1:n))./eErr(1:n)).^2);
    end
end
