function p = AccuracyTableChi2(numberOfParameters, accuracyLevel)
    SigmaLevels = [1,2,3,4,5,6,7,8,9 ;
                    1,1.28,1.64,1.96,2,2.58,3,3.29,4];
    PValueLevels = [1,2,3,4,5,6,7,8,9 ;
                    0.317, 0.2, 0.1, 0.05, 0.0455, 0.01, 0.0027, 0.001, 0.00006];
    AValueLevels = [1,2,3,4,5,6,7,8,9 ;
                    0.683, 0.8, 0.9, 0.95, 0.9545, 0.99, 0.9973, 0.999, 0.9999];
    T =         [1, 1.64, 2.71, 3.84, 4, 6.63, 9, 10.83, 16;
                 2.3, 3.22, 4.61, 5.99, 6.18, 9.21, 11.83, 13.82, 19.33;
                 3.53, 4.64, 6.25, 7.81, 8.02, 11.34, 14.16, 16.27, 22.06;
                 4.72, 5.99, 7.78, 9.49, 9.72, 13.28, 16.25, 18.47, 24.5;
                 5.89, 7.29, 9.24, 11.07, 11.31, 15.09, 18.21, 20.52, 26.77;
                 7.04, 8.56, 10.64, 12.59, 12.85, 16.81, 20.06, 22.46, 28.91;
                 8.18, 9.8, 12.02, 14.07, 14.34, 18.48, 21.85, 24.32, 30.96;
                 9.3, 11.03, 13.36, 15.51, 15.79, 20.09, 23.57, 26.12, 32.93;
                 10.42, 12.24, 14.68, 16.92, 17.21, 21.67, 25.26, 27.88, 34.85;
                 11.54, 13.44, 15.99, 18.31, 18.61, 23.21, 26.9, 29.59, 36.72];
    
    safeNum = max(min(round(numberOfParameters,0),length(T(:,1))),1);
    safeOrd = max(min(round(accuracyLevel,0),length(T(1,:))),1);
    p = T(safeNum,safeOrd);
end
