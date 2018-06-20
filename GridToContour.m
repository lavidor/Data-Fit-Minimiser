function contour = GridToContour(grid,aVectors,aN,contourSizeCounter)
    k = length(aN);
    contour =zeros(contourSizeCounter,k+1);
    coord = cell(1,k);
    i = 0;
    while(i < k)
        i = i+1;
        coord{i} = 0;
    end
    i2 = 1;
    i = 1;
    while(i > 0)
        if(coord{i} < aN(i))
            coord{i} = coord{i}+1;
            if(i < k)
                i = i+1;
            else
                if(grid(coord{:}) > 0)
                    i2 = i2+1;
                    i3 = 0;
                    contour(i2,1) = grid(coord{:});
                    while(i3<k)
                        i3 = i3+1;
                        contour(i2,i3+1) = aVectors(i3,coord{i3});
                    end
                end
                if(grid(coord{:}) < 0)
                    i3 = 0;
                    contour(1,1) = -grid(coord{:});
                    while(i3<k)
                        i3 = i3+1;
                        contour(1,i3+1) = aVectors(i3,coord{i3});
                    end
                end
            end
        else
            coord{i} = 0;
            i = i-1;
        end
    end
end