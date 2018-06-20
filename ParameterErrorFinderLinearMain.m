function pErrors = ParameterErrorFinderLinearMain(x,y,eErr,f,bestP,deltaChi2,errorAccuracyLimit)
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
            k = 3;
            while(k>-1)
            k = k-2;
            dP = k*errorAccuracyLimit(i)*2^50;
            cont = 1;
            pCurrent = bestP(i);
            pPrev = pCurrent;
            while(cont > 0)
                [~,chi2Current] = ParameterFinderLinearKnownMatrix(matrixF,fX,y-pCurrent*f{i}(x),eErr);
                %chi2Current = Chi2CalcLinear(x,y-pCurrent*f{i}(x),eErr,f2,currentBestP);
                if(chi2Current < chi2Limit)
                    pPrev = pCurrent;
                	pCurrent = pCurrent+dP;
                	dP=dP*2;
                else
                	cont = 0;
                end
            end
            pNext = pCurrent;
            pCurrent = (pNext+pPrev)/2;
            while(and(and(pCurrent ~= pNext,pCurrent ~= pPrev),abs(pNext-pPrev) > errorAccuracyLimit(i)))
                [~,chi2Current] = ParameterFinderLinearKnownMatrix(matrixF,fX,y-pCurrent*f{i}(x),eErr);
                %chi2Current = Chi2CalcLinear(x,y-pCurrent*f{i}(x),eErr,f2,currentBestP);
                if(chi2Current < chi2Limit)
                    pPrev = pCurrent;
                else
                	pNext = pCurrent;
                end
                pCurrent = (pNext+pPrev)/2;
            end
            pErrors((5-k)/2,i) = pNext - bestP(i);
            end
            if(i < parameterNum)
                f2(1,i) = f(1,i);
            end
        end
    elseif(parameterNum==1)
        yTheory = ones(1,length(y));
        i = 1;
        k = 3;
        while(k>-1)
            k = k-2;
            dP = k*errorAccuracyLimit(i)*2^50;
            cont = 1;
            pCurrent = bestP(i);
            pPrev = pCurrent;
            while(cont > 0)
                chi2Current = Chi2CalcDirect(y,eErr,pCurrent*yTheory);
                if(chi2Current < chi2Limit)
                    pPrev = pCurrent;
                	pCurrent = pCurrent+dP;
                	dP=dP*2;
                else
                	cont = 0;
                end
            end
            pNext = pCurrent;
            pCurrent = (pNext+pPrev)/2;
            while(and(and(pCurrent ~= pNext,pCurrent ~= pPrev),abs(pNext-pPrev) > errorAccuracyLimit(i)))
                pCurrent = (pNext+pPrev)/2;
                chi2Current = Chi2CalcDirect(y,eErr,pCurrent*yTheory);
                if(chi2Current < chi2Limit)
                    pPrev = pCurrent;
                else
                	pNext = pCurrent;
                end
            end
            pErrors((5-k)/2,i) = pNext - bestP(i);
        end
    end
end
