m = 2048;
A = rand(m, m);
errors = zeros(20, 1);

[L, U, P] = lu(A);

for k = 1:20
    Bk = [rand(m, 1) rand(m, 1)] * [[1; zeros(m-1, 1)] [zeros(m-1, 1); 1]].';
    fk = rand(m, 1);

    xk = stepsolve(P, L, U, Bk, fk);
    xktrue = (A + Bk) \ fk;

    errors(k, :) = norm(xktrue - xk) / norm(xk);
end

error = max(errors)

function xk = stepsolve(P, L, U, Bk, fk)
    temp1 = L \ (P * fk);
    Ainvf = U \ temp1;
    
    X = [Bk(:, 1) Bk(:, end)];
    Y = [[1; zeros(length(X)-1, 1)] [zeros(length(X)-1, 1); 1]].';

    temp2 = L \ (P * X);
    AinvX = U \ temp2;

    xk = Ainvf - AinvX * ((eye(2, 2) + Y * AinvX) \ (Y * Ainvf));
end