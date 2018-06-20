function [p , covM] = ParameterFinderLinearMain(x,y,eErr,f)
    parameterNum = length(f);
    p = zeros(2,parameterNum);
    covM = zeros(parameterNum);
    if (eErr == 0)
        eErr = ones(1,length(x));
    end
    n = min([length(x(:)),length(y(:)),length(eErr(:))]);
    if(n > 0)
        xR = x(1:n);
        yR = y(1:n);
        eR = eErr(1:n);
        fX = zeros(parameterNum,n);
        i = 0;
        while(i < parameterNum)
            i = i+1;
            fX(i,:) = (f{i}(xR))./eR;
        end
        yTemp = yR./eR;
        matrixF = fX*fX';
        vectorY = fX*yTemp';
        p(1,:) = (matrixF\vectorY)';
        covM = (matrixF)^-1;
        i = 0;
        while(i < parameterNum)
            i = i+1;
            p(2,i) = (covM(i,i)).^(1/2);
        end
    end
end
