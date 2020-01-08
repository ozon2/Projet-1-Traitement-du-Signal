load('donnees1.mat');
load('donnees2.mat');

% Initialisation des données
f1 = 0;
f2 = 46000;
Tslot = 0.04;
Fe = 128000;
Te = 1/Fe;
Ns = 10;
Ts = Ns*Te;

%% 3.2.1 Modulation bande base

% Création des signaux m1(t) et m2(t)
m1 = kron(2*bits_utilisateur1 -1, ones(1,Ns));
m2 = kron(2*bits_utilisateur2 -1, ones(1,Ns));

% Tracé des signaux m1(t) et m2(t)
t = linspace(0, Tslot, length(m1));
figure; subplot(2,1,1); plot(t, m1);
title("Signal m1");
xlabel("Temps (s)")
ylabel("Amplitude m1(t)");
ylim([-1.5,1.5]);

subplot(2,1,2); plot(t, m2);
title("Signal m2");
xlabel("Temps (s)")
ylabel("Amplitude m2(t)");
ylim([-1.5,1.5]);

% Calcul des transformées de fourier de m1 et m2
M1 =fft(m1);
M2 =fft(m2);

% Calcul des densités spectrales de puissance
Sm1 = abs(M1).^2;
Sm2 = abs(M2).^2;

% Tracé des densités spectrales de puissance
f = linspace(-Fe/2, Fe/2, length(M1));
figure; subplot(2,1,1); semilogy(f, fftshift(Sm1));
title("DSP de m1");
xlabel("Fréquence (f)")
ylabel("Amplitude Sm1(t)");

subplot(2,1,2); semilogy(f, fftshift(Sm2));
title("DSP de m2");
xlabel("Frequence (f)")
ylabel("Amplitude Sm2(t)");

%% 3.2.2 Construction du signal MF-TDMA

N = length(m1);
% Construction des slots
Signal1 = zeros (1,5*N);
Signal1(N+1:2*N) = m1;

Signal2 = zeros (1,5*N);
Signal2(4*N+1:5*N) = m2;

% Tracé des signaux
t = linspace(0, 5*Tslot, length(Signal1));
figure; subplot(2,1,1); plot(t, Signal1);
title("Signal Utilisateur 1");
xlabel("Temps (s)")
ylabel("Amplitude");
ylim([-1.5,1.5]);

subplot(2,1,2); plot(t, Signal2);
title("Signal Utilisateur 2");
xlabel("Temps (s)")
ylabel("Amplitude");
ylim([-1.5,1.5]);

% Modulation d'amplitude
x1 = Signal1.*cos(2*pi*f1*t);
x2 = Signal2.*cos(2*pi*f2*t);

% signal MF-TDMA
SNRdb = 100;
Ps = mean (abs(x1 + x2).^2);
Pb = Ps/(10^(SNRdb/10));
n = randn (1, length(Signal1))*sqrt(Pb);
x = x1 + x2 + n;

% Tracé de x
figure; plot(t, x);
title("Signal Complet");
xlabel("Temps (s)");
ylabel("Amplitude");

% Calcul de la transformée de fourier de x
X =fft(x);

% Calcul de la DSP de x
Sx = abs(X).^2;

% Tracé de la DSP de x
f = linspace(-Fe/2, Fe/2, length(X));
figure; semilogy(f, fftshift(X));
title("DSP de x");
xlabel("Fréquence (f)")
ylabel("Amplitude Sx(f)");



