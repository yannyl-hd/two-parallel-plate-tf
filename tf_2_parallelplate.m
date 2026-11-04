
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Midterm - Fluid Flow between 2 Parallel Plates
% By: Hannah Douthit
% AE 601
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. Scaling and Equations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% du/dt = -1/rho * dP/dx + nu * d^2u/dy^2
% dU/dtau = G + d^2U/dY^2

clear; clc;
h =1;

%change spacing
dy = 0.1;       %grid spacing; change 0.2 0.1 0.05
% r = 0.25;    % change r value with constant grid spacing
% dt = r*dy^2;
dt = 0.000625;    %change for each run rk4
t = 1;

%nondim
Y = -h:dy:h;
tau = 0:dt:t;

%bc
ubot = 0;
utop = 0;

%ic
U_i = zeros(size(Y));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2. Exact Solution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tau_final = tau(end);
Yex = linspace(-1,1,1000);
Uss = -Yex.^2/2 + 1/2;

N = 500;
tsum = zeros(size(Yex));

for n = 0:N

    m = 2*n + 1;

    term = ((-1)^n/m^3) .* ...
        cos(m*pi*Yex/2) .* ...
        exp(-(m^2*pi^2/4)*tau_final);

    tsum = tsum + term;

end

Utrans = -16/pi^3 * tsum;

U_ex = Uss + Utrans;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3. Numerical Solutions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

r = dt / (dy)^2;
Q = [];

% i. Explicit

N = length(Y);
Uta = zeros(N, 1);
Utn = U_i';
forcing = zeros(N,1);
forcing(2:N-1) = 1;

%matrix
LHS = r;
M = 1-2*r;
RHS = r;
A = zeros(N,N);

%boundary
A(1,1) = 1;
A(N,N) = 1;

for i = 2:N-1
    A(i,i) = M;
    if i > 1
        A(i, i-1) = LHS;
    end

    if i < N
        A(i, i+1) = RHS;
    end
end

%initialize
tstart = 0;
tend = 1;
count = tstart;
tau_hist = [];
ws = [];

%loop
while count < tend
    Uta = A*Utn + dt*forcing;
    count = count + dt;
    Utn = Uta;
    %flow rate
    tau_hist(end+1) = count;
    Q(end+1) = trapz(Y,Uta);
    %wall shear
    ws(end+1) = (Uta(2)-Uta(1))/dy;
end


% ii. Implicit
% matrix
LHSi = -r;
Mi = 1 + 2*r;
RHSi = -r;

A_i = zeros(N,N);

A_i(1,1) = 1;
A_i(N,N) = 1;

for i = 2:N-1
    A_i(i,i) = Mi;
    A_i(i,i-1) = LHSi;
    A_i(i,i+1) = RHSi;
end

%initialize while loop
forcing_i = zeros(N,1);
forcing_i(2:N-1) = 1;

Uta_i = zeros(N,1);
Utn_i = U_i';
Q_i = [];
ws_i = [];
tau_hist_i = [];

count = tstart;

%loop
while count < tend

    rhs_i = Utn_i + dt*forcing_i;
    Uta_i = A_i \ rhs_i;

    count = count + dt;
    Utn_i = Uta_i;

    %flow rate
    tau_hist_i(end+1) = count;
    Q_i(end+1) = trapz(Y,Uta_i);

    %shear stress
    ws_i(end+1) = (Uta_i(2)-Uta_i(1))/dy;
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4. 4th- Order Runge Kutta
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%matrix
N = length(Y);
L = zeros(N,N);

for i = 2:N-1
    L(i,i) = -2;
    L(i,i-1) = 1;
    L(i,i+1) = 1;
end

L = L/dy^2;

%initialize
source = zeros(N,1);
source(2:N-1) = 1;

% loop solving for each of the four point slops
F = @(U) source + L*U;

Utnrk = U_i';
Utark = zeros(N,1);

Qr = [];
wsr = [];
tau_hist_r = [];

count = tstart;

while count < tend

    k1 = F(Utnrk);

    Uhalf = Utnrk + (dt/2)*k1;
    k2 = F(Uhalf);

    Uqua = Utnrk + (dt/2)*k2;
    k3 = F(Uqua);

    Ufull = Utnrk + dt*k3;
    k4 = F(Ufull);

    Utark = Utnrk + (dt/6)*(k1 + 2*k2 + 2*k3 + k4);

    Utnrk = Utark;
    count = count + dt;

    tau_hist_r(end+1) = count;
    Qr(end+1) = trapz(Y,Utark);
    wsr(end+1) = (Utark(2)-Utark(1))/dy;

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 5. Results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%using tau_his to get Qex and Tauex for plotting purposes
Q_exact = zeros(size(tau_hist));
ws_exact = zeros(size(tau_hist));

Nterms = 500;

for k = 1:length(tau_hist)

    tau_k = tau_hist(k);
    sum_Q = 0;
    sum_ws = 0;

    for n = 0:Nterms

        m = 2*n + 1;

        sum_Q = sum_Q + (1/m^4) * ...
            exp(-(m^2*pi^2/4)*tau_k);

        sum_ws = sum_ws + (1/m^2) * ...
            exp(-(m^2*pi^2/4)*tau_k);

    end

    Q_exact(k) = 2/3 - (64/pi^4)*sum_Q;
    ws_exact(k) = 1 - (8/pi^2)*sum_ws;

end


% Velocity Profile

plot(Uta,Y,'r-','LineWidth',1.5)
hold on
plot(Uta_i,Y,'b--','LineWidth',1.5)
plot(Utark,Y,'g-.','LineWidth',1.5)
plot(U_ex,Yex,'k:','LineWidth',1.5)
grid on
xlabel('Dimensionless velocity, U','FontSize',13,'FontWeight','bold')
ylabel('Dimensionless height, Y','FontSize',13,'FontWeight','bold')
legend('Explicit Euler','Implicit Euler','RK4','Exact', ...
    'Location','best')
xlim([0 0.5])


% Flow Rate over Time
figure
plot(tau_hist,Q,'r-','LineWidth',1.5)
hold on
plot(tau_hist_i,Q_i,'b--','LineWidth',1.5)
plot(tau_hist_r,Qr,'g-.','LineWidth',1.5)
plot(tau_hist,Q_exact,'k:','LineWidth',1.5)

xlabel('Dimensionless time, \tau','FontSize',13,'FontWeight','bold')
ylabel('Dimensionless flow rate, Q','FontSize',13,'FontWeight','bold')

legend('Explicit Euler','Implicit Euler','RK4','Exact', ...
    'Location','best')
grid on
pbaspect([1 1 1])


% Wall Shear over Time
figure
plot(tau_hist,ws,'r-','LineWidth',1.5)
hold on
plot(tau_hist_i,ws_i,'b--','LineWidth',1.5)
plot(tau_hist_r,wsr,'g-.','LineWidth',1.5)
plot(tau_hist,ws_exact,'k:','LineWidth',1.5)

xlabel('Dimensionless time, \tau','FontSize',13,'FontWeight','bold')
ylabel('Dimensionless wall shear, \tau_w^*','FontSize',13,'FontWeight','bold')

legend('Explicit Euler','Implicit Euler','RK4','Exact', ...
    'Location','best')

grid on
pbaspect([1 1 1])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 5. Error Comparison
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
U_ex_num = interp1(Yex(:),U_ex(:),Y(:));
% velocity
Ee = max(abs(Uta(:) - U_ex_num(:)));
Ei = max(abs(Uta_i(:) - U_ex_num(:)));
Erk = max(abs(Utark(:) - U_ex_num(:)));
% flow rate
Ee_Q = max(abs(Q - Q_exact));
Ei_Q = max(abs(Q_i - Q_exact));
Erk_Q = max(abs(Qr - Q_exact));

% wall shear
Ee_ws = max(abs(ws - ws_exact));
Ei_ws = max(abs(ws_i - ws_exact));
Erk_ws = max(abs(wsr - ws_exact));

% Error plot
Uta = Uta(:);
Uta_i = Uta_i(:);
Utark = Utark(:);
U_ex_num = U_ex_num(:);

error_explicit = abs(Uta-U_ex_num);
error_implicit = abs(Uta_i-U_ex_num);
error_RK4 = abs(Utark-U_ex_num);

figure
plot(Y(:),error_explicit,'r-','LineWidth',1.5)
hold on
plot(Y(:),error_implicit,'b--','LineWidth',1.5)
plot(Y(:),error_RK4,'g-.','LineWidth',1.5)
hold off

xlabel('Y')
ylabel('Absolute error')
legend('Explicit Euler','Implicit Euler','RK4')
grid on

% grid study values
% spatial loglog plot
dy_values = [0.05, 0.1, 0.2];

error_explicit_grid = [9.5170e-05, 3.7987e-04, 4.4648e-04];
error_implicit_grid = [7.1231e-05, 2.8439e-04, 2.2000e-03];
error_RK4_grid = [1.1963e-05, 4.7637e-05, 8.9235e-04];


figure
loglog(dy_values,error_explicit_grid,'r-o','LineWidth',1.5)
hold on
loglog(dy_values,error_implicit_grid,'b-s','LineWidth',1.5)
loglog(dy_values,error_RK4_grid,'g-^','LineWidth',1.5)

xlabel('\Delta Y','FontSize',13,'FontWeight','bold')
ylabel('Maximum velocity-profile error','FontSize',13,'FontWeight','bold')

legend('Explicit Euler','Implicit Euler','RK4', ...
    'Location','best')
grid on

% Temporal grid-spacing study by changing r

r_values = [0.10, 0.25, 0.40, 0.50];

error_explicit_r = [8.8933e-05, 3.8031e-04, ...
                    3.1133e-04, 4.4648e-04];

error_implicit_r = [3.5561e-04, 2.8394e-04, ...
                    7.5537e-04, 8.8855e-04];

error_RK4_r = [2.2229e-04, 4.8084e-05, ...
               2.2229e-04, 2.2229e-04];

figure
semilogy(r_values,error_explicit_r,'r-o','LineWidth',1.5)
hold on
semilogy(r_values,error_implicit_r,'b-s','LineWidth',1.5)
semilogy(r_values,error_RK4_r,'g-^','LineWidth',1.5)

xlabel('Diffusion number, r','FontSize',13,'FontWeight','bold')
ylabel('Maximum velocity-profile error','FontSize',13,'FontWeight','bold')

legend('Explicit Euler','Implicit Euler','RK4','Location','best')
grid on

%RK4 study
dt_values = [0.005, 0.0025, 0.00125, 0.000625];

RK4_dt_errors = [2.2184e-04, 4.8084e-05, ...
    8.6671e-05, 1.5420e-04];

figure
loglog(dt_values,RK4_dt_errors,'-o','LineWidth',1.5)

xlabel('\Delta\tau','FontSize',13,'FontWeight','bold')
ylabel('Maximum velocity-profile error','FontSize',13,'FontWeight','bold')

grid on
