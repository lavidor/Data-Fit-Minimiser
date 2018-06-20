function [chiGrid ,contourError, aVectors , chi2Min , minIndex] = Chi2GridAndErrorGrid(x,y, e, f, a, aLim, aN, chi2Lim)
    chi2Min = -1;
    contourSize = 1;
    k = length(a);
    aVectors = zeros(k,max(aN));
    if(length(aN) == 1)
        if(aN <= 0)
            aN = 10*ones(1,k);
        else
            aN = aN*ones(1,k);
        end
    end
    i = 0;
    while (i<k)
        i = i+1;
        aVectors(i,1:aN(i)) = linspace(a(i)-aLim(i) , a(i)+aLim(i) , aN(i));
    end
    chiGrid = zeros(aN.*ones(1,k));
    chi2Grid = zeros(aN.*ones(1,k));
    chi2GridBool = zeros(aN.*ones(1,k));
    if(k < 2)
        chiGrid = zeros(1,aN);
    end
    coord = cell(1,k);
    i = 0;
    while(i < k)
        i = i+1;
        coord{i} = 0;
    end
    aArr = zeros(1,k);
    i = 1;
    while(i > 0)
        if(coord{i} < aN(i))
            coord{i} = coord{i}+1;
            aArr(i) = aVectors(i,coord{i});
            if(i < k)
                i = i+1;
            else
                chiGrid(coord{:}) = Chi2CalcDirect(y,e,f(x,aArr));
                if(chi2Min < 0)
                    chi2Min = chiGrid(coord{:})+1;
                end
                if(chi2Min > chiGrid(coord{:}))
                    chi2Min = chiGrid(coord{:});
                    minIndex = coord;
                end
                if(chiGrid(coord{:}) < chi2Lim)
                    chi2GridBool(coord{:}) = 1;
                end
            end
        else
            coord{i} = 0;
            i = i-1;
        end
    end
    i = 0;
    while(i < k)
        i = i+1;
        coord{i} = 0;
    end
    i = 1;
    while(i > 0)
        if(coord{i} < aN(i))
            coord{i} = coord{i}+1;
            if(i < k)
                i = i+1;
            else
                if(chi2GridBool(coord{:}) == 0)
                    if(sumNeighbours(chi2GridBool,coord,aN)>0)
                        contourSize = contourSize+1;
                        chi2Grid(coord{:}) = chiGrid(coord{:});
                    end
                end
            end
        else
            coord{i} = 0;
            i = i-1;
        end
    end
    chi2Grid(minIndex{:}) = -chi2Min;
    contourError = GridToContour(chi2Grid,aVectors,aN,contourSize);
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