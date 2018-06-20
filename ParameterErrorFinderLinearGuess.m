function pErrors = ParameterErrorFinderLinearGuess(x,y,eErr,f,bestP,pErrorGuess,deltaChi2,errorAccuracyLimit)
    chi2Limit = Chi2CalcLinear(x,y,eErr,f,bestP) + deltaChi2;
    parameterNum = min([length(f),length(bestP),length(errorAccuracyLimit)]);
    pErrors = zeros(3,parameterNum);
    pErrors(1,:) = bestP;
    if (eErr == 0)
        eErr = ones(1,length(x));
    end
    if(parameterNum>1)
        f2 = cell(1,parameterNum-1);
        f2(1,:) = f(1,2:parameterNum);
        i = 0;
        while(i<parameterNum)
            i = i+1;
            [matrixF , fX] = ParameterFinderLinearCalculateMatrix(x,eErr,f2);
            k = 0;
            while(k<2)
                k = k+1;
                dP = errorAccuracyLimit(i)*2;
                pCurrent = bestP(i)+ (3-2*k)*pErrorGuess(i);
                [~,chi2Current] = ParameterFinderLinearKnownMatrix(matrixF,fX,y-pCurrent*f{i}(x),eErr);
                smallerChi = 1;
                if(chi2Current > chi2Limit)
                    dP = -dP;
                    smallerChi = 0;
                end
                pPrev = pCurrent;
                pCurrent = pCurrent+dP;
                dP = dP*2;
                cont = 1;
                while(cont > 0)
                    [~,chi2Current] = ParameterFinderLinearKnownMatrix(matrixF,fX,y-pCurrent*f{i}(x),eErr);
                    if(or(and(chi2Current < chi2Limit,smallerChi == 1),and(chi2Current > chi2Limit,smallerChi == 0)))
                        pPrev = pCurrent;
                        pCurrent = pCurrent+dP;
                        dP=dP*2;
                    else
                        cont = 0;
                    end
                end
                pNext = pCurrent;
                if(smallerChi == 0)
                    pNext = pPrev;
                    pPrev = pCurrent;
                end
                pCurrent = (pNext+pPrev)/2;
                while(and(and(pCurrent ~= pNext,pCurrent ~= pPrev),abs(pNext-pPrev) > errorAccuracyLimit(i)))
                    [~,chi2Current] = ParameterFinderLinearKnownMatrix(matrixF,fX,y-pCurrent*f{i}(x),eErr);
                    if(chi2Current < chi2Limit)
                        pPrev = pCurrent;
                    else
                        pNext = pCurrent;
                    end
                    pCurrent = (pNext+pPrev)/2;
                end
                pErrors(k,i) = pNext - bestP(i);
            end
            if(i < parameterNum)
                f2(1,i) = f(1,i);
            end
        end
    elseif(parameterNum==1)
        qP = zeros(1,3);
        qP(1) = sum(((f{1}(x))./eErr).^2);
        qP(2) = -2*sum(((f{1}(x)).*y)./(eErr.^2));
        qP(3) = sum((y./eErr).^2);
        pErrors(1,1) = -qP(2)/2 + (1/2)*(qP(2)*qP(2)-4*qP(1)*qP(3)+4*qP(1)*chi2Limit).^(1/2) - bestP(1);
        pErrors(2,1) = -qP(2)/2 - (1/2)*(qP(2)*qP(2)-4*qP(1)*qP(3)+4*qP(1)*chi2Limit).^(1/2) - bestP(1);
        pErrors = real(pErrors);
    end
end
