function values = op_power(values, operation, power)
    % Operation power.
    % Calculates (operation^power)(values).
    %
    % Examples:
    % op_power(2,          @(x,y)(x*y), 10) = 2^10 = 1024
    % op_power(3,          @(x,y)(x*y), 4)  = 3^4  = 81
    % op_power(5,          @(x,y)(x+y), 7)  = 5*7  = 35
    % op_power([1 2; 2 1], @(A,B)(A*B), 2)  = [5 4; 4 5] - matrix power
    if (power == 0), error('not supported'); end;
    if (power ~= floor(power)), error('power must be integer'); end;
    bi = de2bi(power - 1);
    powN = values;
    for i=1:length(bi)
        if bi(i), values = operation(values, powN); end;
        powN = operation(powN, powN);
    end
end