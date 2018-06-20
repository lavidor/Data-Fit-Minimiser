function [chi2min , pErrors , covM]  = LinearChi2AnalysisPolynomial(orderOfPolynomial,x,y,eErr,accuracyLevel,errorAccuracy)
        f = cell(1,orderOfPolynomial+1);
        i = 0;
        while (i < orderOfPolynomial+1)
            f{i+1} = @(x) x.^i;
            i = i+1;
        end
    
    [chi2min , pErrors , covM] = LinearChi2Analysis(x,y,eErr,f,accuracyLevel,errorAccuracy);
end
