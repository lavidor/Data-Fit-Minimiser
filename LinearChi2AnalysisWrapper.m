function [chi2min , pErrors , covM]  = LinearChi2AnalysisWrapper(functionsName,orderNumOrFunctions,x,y,eErr,accuracyLevel)
    f2 = orderNumOrFunctions;
    if(strcmp(functionsName,'Polynomial'))
        f2 = cell(1,orderNumOrFunctions+1);
        i = 0;
        while (i<orderNumOrFunctions+1)
            f2{i+1} = @(x) x.^i;
            i = i+1;
        end
    elseif(strcmp(functionsName,'Fourier'))
        f2 = cell(1,2*orderNumOrFunctions(1)+1);
        f2{1} = @(x) 1;
        i = 0;
        while (i<orderNumOrFunctions(1)+1)
            i = i+1;
            f2{2*i} = @(x) sin((2*pi()*i/orderNumOrFunctions(2))*x);
            f2{2*i+1} = @(x) cos((2*pi()*i/orderNumOrFunctions(2))*x);
        end
    end
    [chi2min , pErrors , covM]  = LinearChi2Analysis(x,y,eErr,f2,accuracyLevel);
end
