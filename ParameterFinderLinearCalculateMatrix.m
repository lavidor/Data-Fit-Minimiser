function [matrixF , fX] = ParameterFinderLinearCalculateMatrix(x,eErr,f)
    parameterNum = length(f);
    if (eErr == 0)
        eErr = ones(1,length(x));
    end
    matrixF = 0;
    fX = 0;
    n = min([length(x(:)),length(eErr(:))]);
    if(n > 0)
        xR = x(1:n);
        eR = eErr(1:n);
        fX = zeros(parameterNum,n);
        i = 0;
        while(i < parameterNum)
            i = i+1;
            fX(i,:) = (f{i}(xR))./eR;
        end
        matrixF = fX*fX';
    end
end
