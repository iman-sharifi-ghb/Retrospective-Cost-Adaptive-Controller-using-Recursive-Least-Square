function Out = NonLinDynamic(x,u)
    
    global T
    Da = 0.72;
    B = 8;
    gama = 20;
    beta = 0.3;
    
    x1 = x(1);
    x2 = x(2);

    x1_prime = x1+T*(-x1 + Da*(1-x1)*exp(x2/(1+x2/gama)));
    
    x2_prime = x2+T*(-x2 + B*Da*(1-x1)*exp(x2/(1+x2/gama))) + beta*(u-x2);

    Out      = [x1_prime;x2_prime];
    
end