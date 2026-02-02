%% MAE 283A Final Project
clear; close all; clc;

S = load('rand_Pinv_Pbess.mat');
t = S.t(:);
u = S.Pin_INV(:);
y = S.Pout_BESS(:);

N  = length(t);
Ts = median(diff(t));
fs = 1/Ts;

u0 = detrend(u,0);
y0 = detrend(y,0);

z99 = 2.576;                         % ~99% two-sided normal quantile
outPDF = 'MAE283A_FinalProject_Figures.pdf';
if exist(outPDF,'file'), delete(outPDF); end
export_fig = @(fig) exportgraphics(fig, outPDF, 'Append', true);

nfft  = 2^nextpow2(min(N,4096));
nseg  = min(N,2048); if mod(nseg,2), nseg=nseg-1; end
win   = hann(nseg);
nover = floor(0.5*nseg);

%% Fig 1: input autocorrelation
maxLag1 = min(N-1, 5*nseg);
[Ru,lags] = xcorr(u0, maxLag1, 'biased');
tau = lags * Ts;

fig = figure('Color','w');
plot(tau, Ru, 'LineWidth', 1.2); grid on;
xlabel('tau (s)'); ylabel('Ruu');
title('Fig 1: Input autocorrelation');
export_fig(fig);

%% Fig 2: input spectrum
[Pu,f] = pwelch(u0, win, nover, nfft, fs);
w = 2*pi*f;

fig = figure('Color','w');
plot(w, 10*log10(Pu+eps), 'LineWidth', 1.2); grid on;
xlabel('omega (rad/s)'); ylabel('PSD (dB)');
title('Fig 2: Input spectrum');
xlim([0 pi/Ts]);
export_fig(fig);

%% Fig 3: FRF estimate (H1)
[Puy,~] = cpsd(u0, y0, win, nover, nfft, fs);
[Puu,~] = cpsd(u0, u0, win, nover, nfft, fs);
Ghat = Puy ./ Puu;

fig = figure('Color','w');
subplot(2,1,1);
semilogx(w, 20*log10(abs(Ghat)+eps), 'LineWidth', 1.2); grid on;
ylabel('mag (dB)');
title('Fig 3: FRF estimate');

subplot(2,1,2);
semilogx(w, unwrap(angle(Ghat))*180/pi, 'LineWidth', 1.2); grid on;
xlabel('omega (rad/s)'); ylabel('phase (deg)');
export_fig(fig);

%% Fig 4: FIR impulse + 99% CI
nFIR = 150;
idx = (nFIR+1):N;

Phi = zeros(length(idx), nFIR+1);
for k = 0:nFIR
    Phi(:,k+1) = u0(idx-k);
end
yreg = y0(idx);

gFIR = Phi \ yreg;
epsFIR = yreg - Phi*gFIR;

sigma2 = (epsFIR'*epsFIR) / max(length(idx)-(nFIR+1), 1);
Cov = sigma2 * pinv(Phi'*Phi);
ci = z99 * sqrt(max(diag(Cov),0));

k = (0:nFIR).';
fig = figure('Color','w');
plot(k, gFIR, 'LineWidth', 1.2); hold on; grid on;
plot(k, gFIR+ci, '--', 'LineWidth', 1.0);
plot(k, gFIR-ci, '--', 'LineWidth', 1.0);
xlabel('k (samples)'); ylabel('g(k)');
title('Fig 4: FIR impulse');
legend('g','g+CI','g-CI','Location','best');
export_fig(fig);

%% Fig 5: FIR residual vs input (normalized, with 99% bounds)
maxLag5 = 2*nFIR;

u_align   = u0(idx);
eps_align = epsFIR;

% biased cross-correlation
Reu = xcorr(eps_align, u_align, maxLag5, 'biased');

% normalize using zero-lag variances (also biased)
Ree0 = xcorr(eps_align, 0, 'biased');   % = mean(eps^2)
Ruu0 = xcorr(u_align,   0, 'biased');   % = mean(u^2)
ReuN = Reu / sqrt(Ree0 * Ruu0);

lags5 = (-maxLag5:maxLag5).';

sel  = (lags5 >= 0);
tau5 = lags5(sel) * Ts;

M = length(u_align);
conf5 = z99 / sqrt(M);

fig = figure('Color','w');
stem(tau5, ReuN(sel), 'filled'); grid on; hold on;
yline(conf5,'--'); yline(-conf5,'--');
xlabel('tau (s)'); ylabel('corr');
title('Fig 5: FIR residual vs input');
export_fig(fig);


%% Fig 6: Hankel singular values (ERA from FIR/Markov parameters)
% gFIR is your FIR estimate: gFIR(1)=g(0)=D, gFIR(2)=g(1), ...

Kmark = min(120, length(gFIR));   % use first Kmark Markov params (adjust if needed)
g = gFIR(1:Kmark);

% Choose Hankel size so indices are valid
L = floor((Kmark-1)/2);           % ensures 2*L+1 <= Kmark

% Hankel matrices built from g(1), g(2), ... (i.e., excluding D = g(0))
% This matches the demo idea: H uses theta(2:...) and Hbar uses theta(3:...)
H0 = hankel(g(2:L+1),   g(L+1:2*L));      % size L x L
H1 = hankel(g(3:L+2),   g(L+2:2*L+1));    % shifted Hankel (same size)

[U,S,V] = svd(H0,'econ');
sv = diag(S);

fig = figure('Color','w');
semilogy(sv,'o-','LineWidth',1.2); grid on;
xlabel('index'); ylabel('singular value');
title('Fig 6: Hankel singular values');
export_fig(fig);

%% Fig 7: SS impulse vs FIR (ERA realization)
% Pick order (simple threshold on singular values)
nSS = 6;

Ur = U(:,1:nSS);
Vr = V(:,1:nSS);
Sr = S(1:nSS,1:nSS);

% Realization (matches demo)
sqrtSr = sqrt(Sr);
invSqrtSr = diag(1 ./ sqrt(diag(Sr)));

A = invSqrtSr * (Ur' * H1 * Vr) * invSqrtSr;

H1fac = Ur * sqrtSr;
H2fac = sqrtSr * Vr';

D = gFIR(1);          % D = g(0)
C = H1fac(1,:);       % first row
B = H2fac(:,1);       % first column

% Build SS impulse response gSS(k)
Kimp = nFIR;
k = (0:nFIR).';

gSS = zeros(Kimp+1,1);
gSS(1) = D;

Ak = eye(nSS);
for kk = 2:(Kimp+1)
    gSS(kk) = C * Ak * B;   % k=1 => C*I*B = C B
    Ak = Ak * A;
end

fig = figure('Color','w');
plot(k, gFIR(1:Kimp+1), 'LineWidth', 1.2); hold on; grid on;
plot(k, gSS, '--', 'LineWidth', 1.2);
xlabel('k (samples)'); ylabel('g(k)');
title(sprintf('Fig 7: SS vs FIR (nSS=%d)', nSS));
legend('FIR','SS','Location','best');
export_fig(fig);


%% Fig 8: SS residual vs input (R_{e u} with 99% bounds)  [PROJECT FORMAT]
rho = max(abs(eig(A)));     % discrete-time stability check
fprintf('nSS=%d, max|eig(A)|=%.4f\n', nSS, rho);

if rho >= 0.995
    warning('Realized A is (near) unstable -> reduce nSS (try nSS-1, nSS-2, ...)');
end

sysSS = ss(A,B,C,D,Ts);

% Discrete-time simulation (no time vector)
ySS   = lsim(sysSS, u0);
epsSS = y0 - ySS;

maxLag8 = 2*nSS;

% Biased cross-correlation estimate (matches demo style)
[Reu8, lags8] = xcorr(epsSS, u0, maxLag8, 'biased');

sel  = (lags8 >= 0);
tau8 = lags8(sel) * Ts;
Reu8 = Reu8(sel);

Ncorr = length(u0);
conf8 = z99 / sqrt(Ncorr);

fig = figure('Color','w');
stem(tau8, Reu8, 'filled'); grid on; hold on;
yline(conf8,'--'); yline(-conf8,'--');
xlabel('tau (s)');
ylabel('$\hat R_{\epsilon u}(\tau)$','Interpreter','latex');
title(sprintf('Fig 8: SS residual vs input (nSS=%d)', nSS));
export_fig(fig);


%% Fig 9: ARX residual vs input (PROJECT FORMAT)
% Uses existing: u0, y0, Ts, N, z99, export_fig
% Chooses low-order ARX model that minimizes violations of 99% bounds
% on normalized residual-input correlation for 0 <= tau <= 2n.

data = iddata(y0, u0, Ts);

naSet = 1:6;
nbSet = 1:6;
nkSet = 0:3;

bestScore = Inf;
bestM     = [];
bestOrd   = [1 1 0];
bestn     = 1;

% For confidence bounds use the effective sample size used in residual test
M = length(u0);
conf9 = z99 / sqrt(M);

for na = naSet
    for nb = nbSet
        for nk = nkSet

            m = arx(data, [na nb nk]);

            % prediction-error residual epsilon(t)
            E = resid(data, m);
            e = E.OutputData(:);

            % Define "model order n" for the validation window 0 <= tau <= 2n
            % (includes delay effect)
            n = max(na, nb + nk);
            maxLag = 2*n;

            % Biased cross-correlation estimate
            Reu = xcorr(e, u0, maxLag, 'biased');

            % Normalize to correlation using zero-lag variances (biased)
            Ree0 = xcorr(e,  0, 'biased');   % mean(e^2)
            Ruu0 = xcorr(u0, 0, 'biased');   % mean(u^2)
            ReuN = Reu / sqrt(max(Ree0,eps) * max(Ruu0,eps));

            lags = (-maxLag:maxLag).';

            % Keep only 0 <= tau <= 2n
            sel = (lags >= 0);
            Reu_pos = ReuN(sel);

            % Score = total amount of 99% bound violation + tiny complexity penalty
            viol = sum(max(0, abs(Reu_pos) - conf9));
            score = viol + 0.01*(na + nb + nk);

            if score < bestScore
                bestScore = score;
                bestM   = m;
                bestOrd = [na nb nk];
                bestn   = n;
            end
        end
    end
end

% Recompute for best model and plot
E9 = resid(data, bestM);
e9 = E9.OutputData(:);

maxLag9 = 2*bestn;

Reu9 = xcorr(e9, u0, maxLag9, 'biased');
Ree0 = xcorr(e9, 0, 'biased');
Ruu0 = xcorr(u0, 0, 'biased');
Reu9N = Reu9 / sqrt(max(Ree0,eps) * max(Ruu0,eps));

lags9 = (-maxLag9:maxLag9).';
sel9  = (lags9 >= 0);
tau9  = lags9(sel9) * Ts;

fig = figure('Color','w');
stem(tau9, Reu9N(sel9), 'filled'); grid on; hold on;
yline(conf9,'--'); yline(-conf9,'--');
xlabel('tau (s)');
ylabel('corr');
title(sprintf('Fig 9: ARX residual vs input (na=%d nb=%d nk=%d), n=%d', ...
      bestOrd(1), bestOrd(2), bestOrd(3), bestn));
export_fig(fig);

fprintf('Fig 9 best ARX order: na=%d nb=%d nk=%d (n=%d), score=%.4g, bound=±%.4g\n', ...
    bestOrd(1), bestOrd(2), bestOrd(3), bestn, bestScore, conf9);

%% Fig 10: y vs ysim vs ypred
ysim = sim(bestM, data).OutputData(:);
ypred = predict(bestM, data, 1).OutputData(:);

fig = figure('Color','w');
plot(t, y0, 'LineWidth', 1.1); hold on; grid on;
plot(t, ysim, '--', 'LineWidth', 1.1);
plot(t, ypred, ':', 'LineWidth', 1.3);
xlabel('t (s)'); ylabel('y');
title('Fig 10: y vs ysim vs ypred');
legend('y','ysim','ypred','Location','best');
export_fig(fig);

disp('DONE: exported MAE283A_FinalProject_Figures.pdf');
