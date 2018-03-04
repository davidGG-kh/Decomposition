A = [  1  1  0  0  0  0  0 ;
    -1  0  1  1  0  0  0 ;
    0 -1 -1  0  1  1  0 ;
    0  0  0 -1 -1  0  1 ;
    0  0  0  0  0 -1 -1 ];
c = [  1; 1; 1; 1; 1; 1; 1];
s = [ 0.2; 0.6; 0; 0; -0.8 ];
[m,n] = size(A);
x0 = linprog(zeros(n,1),...
[eye(n) 
    -eye(n)], [c;zeros(n,1)],...
    A, s);

%%


alpha = 1;

mu = rand(m,1);
x = x0;
for k = 1:100
    for j=1:n
        if mu'*A(:,j)>-1/c(j)
            x(j) = 0;
        else
            x(j) = c(j) - sqrt(-c(j)/(mu'*A(:,j)));
        end
    end
    g = A*x - s;
    mu = mu + alpha * g;
    fprintf('iteration %d : obj = %g, primal = %g\n', k, sum(x./(c-x)), norm(A*x-s));
end