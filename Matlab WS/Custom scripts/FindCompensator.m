function coefs = FindCompensator(func, coefs, lb, ub, A, b, Aeq, beq, nonlcon)
%FINDCOMPENSATOR Summary of this function goes here

    arguments
        func 
        coefs 
        lb = []
        ub = []
        A = []
        b = [] 
        Aeq = []
        beq = []
        nonlcon= [] 
    end

    %options = optimset('PlotFcns',@optimplotfval, 'MaxFunEvals', 10000, 'MaxIter', 10000);
    %options = optimset('Display', 'iter', 'MaxFunEvals', 10000, 'MaxIter', 10000);
    options = optimset('MaxFunEvals', 100000, 'MaxIter', 100000);
    
    [coefs, ~, ~, exitflag, ~] = fminimax(func, coefs, A, b, Aeq, beq, lb, ub, nonlcon, options);
    %[coefs, ~, exitflag, ~] = fmincon(func, coefs, A, b, Aeq, beq, lb, ub, nonlcon, options);
    
    if (exitflag ~= 1)
        fprintf('ERROR - Could not find solution\n');
    else
        fprintf('OK - Solution found within Tol\n');
    end
end

