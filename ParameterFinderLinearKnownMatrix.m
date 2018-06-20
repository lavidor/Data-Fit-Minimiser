function [p,chi2] = ParameterFinderLinearKnownMatrix(matrixF,fX,y,eErr)
    parameterNum = length(fX);
    yR = y(1:length(fX));
    eR = eErr(1:length(fX));
    yTemp = yR./eR;
    vectorY = fX*yTemp';
    p = zeros(2,parameterNum);
    p(1,:) = (matrixF\vectorY)';
    chi2 = sum((p*fX - yTemp).^2);
    covM = (matrixF)^-1;
    i = 0;
    while(i < parameterNum)
        i = i+1;
        p(2,i) = (covM(i,i)).^(1/2);
    end
end
