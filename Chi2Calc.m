function p = Chi2Calc(x,y,eErr,f)
    p = Chi2CalcDirect(y,eErr,f(x));
end
