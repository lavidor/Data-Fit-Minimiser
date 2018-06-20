function [p , covM] = Linear3ParameterFinder(x,y,eErr,f)
    parameterNum = length(f);
    p = zeros(2,parameterNum);
    covM = zeros(parameterNum);
    eErr2 = eErr;
    if (eErr == 0)
        eErr2 = ones(1,length(x));
    end
    fR = f;
    n = min([length(x(:)),length(y(:)),length(eErr2(:))]);
    if(n > 0)
        xR = x(1:n);
        yR = y(1:n);
        eR = eErr2(1:n);
        fX = zeros(parameterNum,n);
        i = 0;
        while(i < parameterNum)
            i = i+1;
            fX(i,:) = (fR{i}(xR))./eR;
        end
        yTemp = yR./eR;
        matrixF = fX*fX';
        vectorY = fX*yTemp';
        p = (matrixF\vectorY)';
        covM = (matrixF)^-1;
        
        ff = matrixF(1,1);
        fy = vectorY(1);
        
        if(parameterNum == 1)
            denom = ff;
            p(1,1) = fy;
            if(denom ~= 0)
                p(1,1) = p(1,1)/denom;
            end
        else
            gg = matrixF(2,2);
            fg = matrixF(1,2);
            gy = vectorY(2);

            if(parameterNum == 2)
                denom = fg*fg-ff*gg;
                p(1,1) = fg*gy-gg*fy;
                p(1,2) = fg*fy-ff*gy;
                if(denom ~= 0)
                    p(1,1) = p(1,1)/denom;
                    p(1,2) = p(1,2)/denom;
                end
            else
                hh = matrixF(3,3);
                fh = matrixF(1,3);
                gh = matrixF(2,3);
                hy = vectorY(3);
                
                if(parameterNum == 3)
                    denom = ff*gg*hh-fh*fh*gg+2*fg*gh*fh-fg*fg*hh-gh*gh*ff;
                    p(1,1) = fh*gh*gy-fg*hh*gy+gg*hh*fy-gh*gh*fy+fg*gh*hy-fh*gg*hy;
                    if(denom ~= 0)
                        p(1,1) = p(1,1)/denom;
                    end
                    denom = hh*fg-fh*gh;
                    p(1,2) = p(1,1)*fh*fh - p(1,1)*ff*hh + hh*fy - fh*hy;
                    if(denom ~= 0)
                        p(1,2) = p(1,2)/denom;
                    end
                    denom = fh;
                    p(1,3) = fy - p(1,1)*ff - p(1,2)*fg;
                    if(denom ~= 0)
                        p(1,3) = p(1,3)/denom;
                    end
                end
            end
        end
    end
    k = 0;
    while(k<parameterNum)
        k = k+1;
        p(2,k) = sqrt(covM(k,k));
    end
end
