function [chi2min , dof, p , chiGrid , contourGridError , minIndex , aVectors , aN] = ParameterFinderNonLinearGridMain(x,y,eErr,f,aEst,aLim,repNum,accLevel,gridNumMatrix)
    contourGridError = 0;
    contourSizeCounter = 1;
    k = length(aEst);
    if(accLevel<=0)
        accLevel = 1;
    end
    if (repNum <= 0)
        repNum = 2;
    end
    dChi2 = AccuracyTableChi2(k, accLevel);
    gridN = gridNumMatrix;
    if(length(gridNumMatrix(1,:)) == 1)
        gridN = 10*ones(repNum,k);
        if(length(gridNumMatrix(:,1)) == 1)
            if(gridNumMatrix > 0)
                gridN = gridNumMatrix*ones(repNum,k);
            end
        else
            i = 0;
            while (i < length(gridNumMatrix(:,1)))
                i = i+1;
                gridN(i,:) = gridNumMatrix(i,1)*ones(1,k);
            end
        end
    elseif(length(gridNumMatrix(:,1)) == 1)
        gridN = 10*ones(repNum,k);
        i = 0;
        while (i < length(gridNumMatrix(1,:)))
            i = i+1;
            gridN(:,i) = gridNumMatrix(1,i)*ones(repNum,1);
        end
    end
    chiGrid = 0;
    chiGrid2 = 0; 
    chi2min = -1;
    if (eErr == 0)
        eErr = ones(1,length(x));
    end
    n = min([length(x(:)),length(y(:)),length(eErr(:))]);
    dof = n - k;
    if(n > 0)
        xR = x(1:n);
        yR = y(1:n);
        eR = eErr(1:n);
        repCounter = 0;
        while(repCounter < repNum)
            repCounter = repCounter+1;
            contourSizeCounter = 1;
            p = zeros(4,k);
            aN = gridN(repCounter,:);
            [chiGrid , aVectors , chi2min , minIndex] = Chi2Grid(xR,yR, eR, f, aEst, aLim, aN);
            m = 0;
            aLim2 = aLim;
            while(m < k)
                m = m+1;
                aEst(m) = aVectors(m,minIndex{m});
                if(aN(1,m) == 1)
                    aLim2(m) = 0;
                elseif(aN(1,m) == 2)
                    aLim2(m) = (aVectors(m,2)-aVectors(m,1))/2;
                    aEst(m) = (aVectors(m,1)+aVectors(m,2))/2;
                else
                    aLim2(m) = (aVectors(m,2)-aVectors(m,1));
                    if(minIndex{m} >= aN(1,m))
                        minIndex{m} = minIndex{m}-1;
                    end
                    if(minIndex{m} <= 1)
                        minIndex{m} = 2;
                    end
                    aEst(m) = aVectors(m,minIndex{m});
                end
            end
            [p(1,:) , chi2min] = Chi2FindMin(xR,yR, eR, f, aEst, aLim2,10^-15);
            aEst(1,:) = p(1,:);
            chi2Lim = chi2min + dChi2;
            coord = cell(1,k);
            aLim = zeros(1,k);
            m = 0;
            while(m < k)
                m = m+1;
                coord{m} = 0;
            end
            chiLimGrid = chiGrid.*0+chi2Lim;
            chiGrid2Bool = chiGrid<chiLimGrid;
            chiGrid2Bool(minIndex{:}) = 1;
            chiGrid2 = 0*chiGrid2Bool;
            i = 1;
            while(i > 0)
                if(coord{i} < aN(1,i))
                    coord{i} = coord{i}+1;
                    if(i < k)
                        i = i+1;
                    else
                        if(chiGrid2Bool(coord{:}) == 0)
                            if(sumNeighbours(chiGrid2Bool , coord , aN(1,:))>0)
                                contourSizeCounter = contourSizeCounter+1;
                                chiGrid2(coord{:}) = chiGrid(coord{:});
                            end
                        end
                        if(or(chiGrid2Bool(coord{:}) == 1 , chiGrid2(coord{:}) > 0))
                            m = 0;
                            while(m < k)
                                m = m+1;
                                eTemp = aVectors(m,coord{m}) - aVectors(m,minIndex{m});
                                if(abs(eTemp) > aLim(m))
                                    aLim(m) = abs(eTemp);
                                    p(2,m) = abs(eTemp);
                                end
                                if(eTemp > p(3,m))
                                    p(3,m) = eTemp;
                                end
                                if(eTemp < p(4,m))
                                    p(4,m) = eTemp;
                                end
                            end
                        end
                    end
                else
                    coord{i} = 0;
                    i = i-1;
                end
            end
            chiGrid2(minIndex{:}) = -chi2min;
            m = 0;
            while(m < k)
                m = m+1;
                if(aEst(1,m) + aLim(m) > aVectors(m,aN(1,m)))
                    aEst(1,m) = (aEst(1,m) - aLim(m) + aVectors(m,aN(1,m)))/2;
                    aLim(m) = aVectors(m,aN(1,m)) - aEst(1,m);
                end
                if(aEst(1,m) - aLim(m) < aVectors(m,1))
                    aEst(1,m) = (aEst(1,m) + aLim(m)+ aVectors(m,1))/2;
                    aLim(m) = aEst(1,m)-aVectors(m,1);
                end
            end
        end
        contourGridError = GridToContour(chiGrid2,aVectors,aN,contourSizeCounter);
    end
end

function s = sumNeighbours(gridC , coordArr , nArr)
    s = 0;
    i = 0;
    while i < length(coordArr)
        i = i+1;
        if(coordArr{i} > 1)
            coordArr{i} = coordArr{i}-1;
            s = s + gridC(coordArr{:});
            coordArr{i} = coordArr{i}+1;
        end
        if(coordArr{i} < nArr(i))
            coordArr{i} = coordArr{i}+1;
            s = s + gridC(coordArr{:});
            coordArr{i} = coordArr{i}-1;
        end
    end
end
