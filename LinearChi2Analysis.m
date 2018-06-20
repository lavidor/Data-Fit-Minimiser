function [chi2min , pErrors , covM]  = LinearChi2Analysis(x,y,eErr,f,accuracyLevel)
    accLevel = 1;
    if(accuracyLevel > 0)
        accLevel = accuracyLevel;
    end
    eErr2 = eErr;
    if (eErr == 0)
        eErr2 = ones(length(x));
    end
    [pErrors , covM] = ParameterFinderLinearMain(x,y,eErr2,f);
    chi2min = Chi2CalcLinear(x,y,eErr2,f,pErrors(1,:));
    deltaChi2 = AccuracyTableChi2(length(pErrors(1,:)), accLevel);
    pErrors(2,:) = pErrors(2,:).*(deltaChi2.^(1/2));
end
