function p = Chi2CalcLinear(x,y,eErr,fArr,a)
    yTheory = zeros(1,length(x));
    parameterNum = min(length(fArr),length(a));
    i = 0;
    while (i<parameterNum)
        i = i+1;
        yTheory = yTheory+a(i)*fArr{i}(x);
    end
    p = Chi2CalcDirect(y,eErr,yTheory);
end
