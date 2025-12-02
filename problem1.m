% Problem 1 - long step interior point algorithm
Q = []; % symmetric PD matrix 
c = []; % Rn vector
A = []; % R(m x n) matrix, full row rank
b = []; % Rm vector
n = length(c);
[m,] = size(A);

% hyperparamaters
delta = ;
eps_max = ;
eps_min = ;

% initial vectors
x = ones(n);
y = ones(m); 
z = ones(n);

ok = TRUE;
k = 0;
eps = 0;
beta = 0;
h = 0.00001;

while norm(F) < 0.000001
    if mod(k,2) == 0
        eps = eps_max;
    else
        eps = eps_min;
    end

    F = [-Q*x - c + A.'*y + z;
        A*x - b;
        z .* x - tau * ones(n)];

    % Jacobian of F
    nab_F = [Q, A.', eye(n);
             A, zeros(m,n), zeros(m,n);
             diag(z), zeros(n), diag(x)];

    %  ∇F [Δx; Δy; Δz] = -F 
    delta_vector = nab_F \ F;

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
        condition = x_trial .* z_trial - beta_trial*delta*ones(n,1);

        if condition >= 0
            alpha = alpha_trial;  % still inside neighborhood
        else
            break;     % outside neighborhood
        end
    end

    % update iterate
    x = x + alpha*dx;
    y = y + alpha*dy;
    z = z + alpha*dz;

    k = k + 1;
        
end

