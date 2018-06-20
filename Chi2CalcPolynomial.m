function p = Chi2CalcPolynomial(x,y,eErr,a)
    yTheory = zeros(1,length(x));
    parameterNum = length(a);
    i = 0;
    while (i<parameterNum)
        i = i+1;
        yTheory = yTheory+a(i)*x.^(i+1);
    end
    p = Chi2CalcDirect(y,eErr,yTheory);
end
