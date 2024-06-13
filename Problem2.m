m = 4000;
c = rand(m,1);
C = toeplitz(c, [c(1); flip(c(2:end))]);
a = rand(m, 1);
b = rand(m, 1);
v = rand(m, 1);

t1 = tic;
xtrue = (C + a * b.') \ v;
trueTime = toc(t1)

t2 = tic;
x = solveSM(C, a, b, v);
SMtime = toc(t2)

accuracy = norm(x - xtrue)


function x = solveSM(C, a, b, v)

    Linv = 1 ./ fft(C(:, 1));
    Cinvv = ifft(Linv .* fft(v));
    Cinva = ifft(Linv .* fft(a));

    x = Cinvv - ((1 / (1 + b.' * Cinva)) * (Cinva * b.')) * Cinvv;

end