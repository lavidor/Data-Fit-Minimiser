function [p , covM] = ParameterFinderLinearWrapper(name,x,y,eErr,f)
    p = 0;
    covM = 0;
    f2 = f;
    if(strcmp(name,'Polynomial'))
        polinomeOrder = f;
        if(length(f)>1)
            polinomeOrder = length(f);
        end
        f2 = cell(1,polinomeOrder+1);
        i = 0;
        while (i < polinomeOrder+1)
            f2{i+1} = @(x) x.^i;
            i = i+1;
        end
        [p , covM] = ParameterFinderLinearMain(x,y,eErr,f2);
    elseif(strcmp(name,'Quadratic'))
        [p , covM] = ParameterFinderLinearWrapper('Polynomial',x,y,eErr2,2);
    elseif(STRCMP(name,'Linear'))
        [p , covM] = ParameterFinderLinearWrapper('Polynomial',x,y,eErr2,1);
    elseif(STRCMP(name,'Constant'))
        [p , covM] = ParameterFinderLinearWrapper('Polynomial',x,y,eErr2,0);
    else
        [p , covM] = ParameterFinderLinearMain(x,y,eErr,f2);
    end
end
