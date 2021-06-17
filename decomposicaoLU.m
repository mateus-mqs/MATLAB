clear all; clc;
A = [1 -3 2; -2 8 -1; 4 -6 5];
m = size(A,1);
n = size(A,2);
L = zeros(m,n);
U = zeros(m,n);
pivo = [0 0 0];
for i = 1:m
    pivo(i) = i;
end
PdU = 1; Info = 0;
for j = 1:min(m,n)
    p = j;
    Amax = abs(A(j,j));
    for k = j+1:m
       if abs(A(k,j)) > Amax
           Amax = abs(A(k,j));
           p = k;
       end
    end
    if p ~= j
        for k = 1:n
            t = A(j,k);
            A(j,k) = A(p,k);
            A(p,k) = t;
        end
        i = pivo(j);
        pivo(j) = pivo(p);
        pivo(p) = i;
        PdU = -PdU;
    end
    PdU = PdU*A(j,j);
    if abs(A(j,j)) ~= 0
        r = 1/A(j,j);
        for i = j+1:m
            Mult = A(i,j)*r;
            A(i,j)= Mult;
            for k = j+1:n
                A(i,k) = A(i,k) - Mult*A(j,k);
            end
        end
    else
        if Info == 0
            Info = j;
        end
    end
end

for i = 1:m
    for j = 1:n
        if i == j
            L(i,j) = 1;
            U(i,j) = A(i,j);
        elseif i < j
            L(i,j) = 0;
            U(i,j) = A(i,j);
        elseif i > j
            L(i,j) = A(i,j);
            U(i,j) = 0;
        end
    end
end

disp(pivo);
disp(U);
disp(L);
disp(PdU);
disp(Info);