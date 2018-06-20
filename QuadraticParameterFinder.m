function [p , covM] = QuadraticParameterFinder(x,y,eErr)
    [p , covM] = ParameterFinderLinearWrapper('Polinomial',x,y,eErr,2);
end
