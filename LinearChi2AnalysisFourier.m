function [chi2min , a0Errors, apErrors, bpErrors , covM]  = LinearChi2AnalysisFourier(orderOfFourier,l,x,y,eErr,accuracyLevel,errorAccuracy)
    if(l == 0)
        l = 2*pi();
    end
        f = cell(1,2*orderOfFourier+1);
        f{1} = @(x) 1;
        i = 1;
        while (i < orderOfFourier+1)
            f{2*i} = @(x) sin((2*pi()*i/l)*x);
            f{2*i+1} = @(x) cos((2*pi()*i/l)*x);
            i = i+1;
        end
    
    [chi2min , pErrors , covM] = LinearChi2Analysis(x,y,eErr,f,accuracyLevel,errorAccuracy);
    a0Errors = pErrors(:,1);
    apErrors = pErrors(:,3:2:end);
    bpErrors = pErrors(:,2:2:end);
    
end
