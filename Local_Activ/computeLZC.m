function LZComplexity = computeLZC(series)
%COMPUTELZC Computes the Lempel-Ziv complexity of the input sequence.
% This function does not apply any filtering; if required, filter the
% signal beforehand.

if any(isnan(series))
    LZComplexity = NaN;
else
    n = length(series);
    b = n / log2(n);

    % Binarize the series using the median as threshold
    threshold = median(series);
    for j = 1:n
        if series(j) >= threshold
            series(j) = 1;
        else
            series(j) = 0;
        end
    end

    % Initialize values to compute the L-Z complexity
    c = 1;                  % initial complexity
    S = series(1);
    Q = series(2);

    for i = 2:n
        SQ = [S, Q];
        SQ_pi = SQ(1:(length(SQ) - 1));

        % Indices at which Q starts within SQ_pi
        k = findstr(Q, SQ_pi);

        if length(k) == 0
            % Q is not present in the series: new sequence
            c = c + 1;
            if (i + 1) > n
                break;
            else
                S = [S, Q];             % build S
                Q = series(i + 1);      % update Q
            end
        else
            % Q is present in the series: not a new sequence
            if (i + 1) > n
                break;
            else
                Q = [Q, series(i + 1)]; % update Q
            end
        end
    end

    LZComplexity = c / b;
end
end
