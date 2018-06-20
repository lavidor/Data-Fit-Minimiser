function [chi2min , pErrors , covM]  = LinearChi2AnalysisUserDesignFunc(x,y,eErr,accuracyLevel,errorAccuracy)
    %define here
    f = cell(1,5);
    f{1} = @(x) 1;
    f{2} = @(x) exp(x);
    f{3} = @(x) exp(-x);
    f{4} = @(x) 1/x;
    f{5} = @(x) sqrt(x);
    %stop defining
    
    accLevel = 1;
    if(exists(accuracyLevel) == 1)
        accLevel = accuracyLevel;
    end
    accError = 2^-15;
    if(exists(errorAccuracy) == 1)
        accError = errorAccuracy;
    end
    eErr2 = ones(length(x));
    if(exists(eErr) == 1)
        if(eErr ~= 0)
            eErr2 = eErr;
        end
    end
    [chi2min , pErrors , covM] = LinearChi2Analysis(x,y,eErr2,f,accLevel,accError);
    
end
