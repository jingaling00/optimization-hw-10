% Problem 2
Q = zeros(20);
%c = [-1.8; -1.4; -1.2; -1; 0; 0; 0]; % Rn vector
%A = [0 0 6 10 1 0 0;
%    0 12 8 6 0 1 0;
%    16 4 2 0 0 0 1]; % R(m x n) matrix, full row rank
%b = [17600; 21120; 1320]; % Rm vector
A = A;
b = b;
c = c;
[m,n] = size(A);

% hyperparamaters
delta = 0.1;
eps_max = 0.99;
eps_min = 0.01;

% initial vectors
x = ones(n,1);
y = ones(m,1); 
z = ones(n,1);

k = 0;
eps = 0;
beta = 0;
h = 0.00001;
F = inf;
tol = 1e-6;

while true
    F_un = [ -Q*x - c + A.'*y + z;
              A*x - b;
              x .* z ];
  
    if mod(k,2) == 0
        eps = eps_max;
    else
        eps = eps_min;
    end

    F = [-Q*x - c + A.'*y + z;
        A*x - b;
        z .* x - eps * beta * ones(n,1)];

    if norm(F) < tol
        fprintf('Done, converged with %d iterations\n', k);
        fprintf('Optimal X* = \n');
        disp(x);
        fprintf('Associated KKT multipliers:\n');
        fprintf('y* = \n'); disp(y); 
        fprintf('z* = \n'); disp(z);
        break;
    end

    % Jacobian of F
    nab_F = [-Q, A.', eye(n);
             A, zeros(m,m), zeros(m,n);
             diag(z), zeros(n,m), diag(x)];

    %  nabla F [delta x; delta y; delta z] = -F 
    delta_vector = -nab_F \ F;

    dx = delta_vector(1:n);
    dy = delta_vector(n+1:n+m);
    dz = delta_vector(n+m+1:end);

    % finding Newton step size
    alpha = 0;
    % alpha = 0, h, 2h, ... up to 1
    for t = 0:floor(1/h)
        alpha_trial = t * h;
        if alpha_trial > 1 + 1e-12 
            break;
        end

        x_trial = x + alpha_trial*dx;
        z_trial = z + alpha_trial*dz;

        % must keep strict positivity
        if any(x_trial <= 0) || any(z_trial <= 0)
            break;
        end

        % check the conditions for N_inf(delta)
        beta_trial = (x_trial' * z_trial) / n;
        if all(x_trial .' * z_trial >= delta * beta_trial)
            alpha = alpha_trial;
        else
            break;
        end
    end

    % update iterate
    x = x + alpha*dx;
    y = y + alpha*dy;
    z = z + alpha*dz;

    k = k + 1;
end

