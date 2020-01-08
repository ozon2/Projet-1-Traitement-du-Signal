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
xlabel("Fréquence (Hz)")
ylabel("Amplitude Sm1(t)");

subplot(2,1,2); semilogy(f, fftshift(Sm2));
title("DSP de m2");
xlabel("Frequence (Hz)")
ylabel("Amplitude Sm2(t)");

%% 3.2.2 Construction du signal MF-TDMA

N = length(m1);
% Construction des slots
Signal1 = zeros (1,5*N);
Signal1(N+1:2*N) = m1;

Signal2 = zeros (1,5*N);
Signal2(4*N+1:5*N) = m2;

NP = 2^nextpow2(length(Signal1));

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
X =fft(x, NP);

% Calcul de la DSP de x
Sx = abs(X).^2;

% Tracé de la DSP de x
f = linspace(-Fe/2, Fe/2, length(Sx));
figure; semilogy(f, fftshift(abs(Sx)));
title("DSP de x");
xlabel("Fréquence (Hz)")
ylabel("Amplitude Sx(f)");

%% 4.1.1 Synthèse du filtre passe-bas

% Création filtre passe bas et TF
fc = f2/2;
N = 101;
k = (-(N-1)/2 : (N-1)/2);
filtre_bas = 2*(fc/Fe)*sinc(2*k*(fc/Fe));
TF_filtre_bas = fft(filtre_bas, NP);

% Tracé de la réponse impulsionnelle
figure;
plot (k,filtre_bas);
title("Tracé de la réponse impulsionnelle du filtre bas");
xlabel("Temps (s)");
ylabel("Réponse impulsionnelle");

% Tracé de la réponse en fréquence
figure; 
semilogy (linspace(-Fe/2, Fe/2, length(TF_filtre_bas)), fftshift(abs(TF_filtre_bas)));
title("Tracé de la réponse en fréquence du filtre bas");
xlabel("Frequence (Hz)");
ylabel("Réponse en fréquence");

% Verification du filtre
Sxnormalise = ((1/max(abs(Sx))) * abs(Sx));
figure;
semilogy (linspace(-Fe/2, Fe/2, length(TF_filtre_bas)), fftshift(abs(TF_filtre_bas))); hold;
semilogy (f, fftshift(abs(Sxnormalise)));
title("Verification du filtre passe-bas");
xlabel("Frequence (Hz)");
ylabel("Amplitude")
legend("Filtre passe-bas","DSP Normalisée")

%% 4.1.2 Synthèse du ltre passe-haut

% Création filtre passe haut et TF
filtre_haut = -filtre_bas;
filtre_haut((N-1)/2+1) = 1 - filtre_bas((N-1)/2+1);
TF_filtre_haut = fft(filtre_haut, NP);

% Tracé de la réponse impulsionnelle
figure;
plot (k,filtre_haut);
title("Tracé de la réponse impulsionnelle du filtre haut");
xlabel("Temps (s)");
ylabel("Réponse impulsionnelle");

% Tracé de la réponse en fréquence
figure; 
semilogy (linspace(-Fe/2, Fe/2, length(TF_filtre_haut)), fftshift(abs(TF_filtre_haut)));
title("Tracé de la réponse en fréquence du filtre haut");
xlabel("Frequence (Hz)");
ylabel("Réponse en fréquence");

% Verification du filtre
figure;
semilogy (linspace(-Fe/2, Fe/2, length(TF_filtre_haut)), fftshift(abs(TF_filtre_haut))); hold;
semilogy (f, fftshift(abs(Sxnormalise)));
title("Verification du filtre passe-haut");
xlabel("Frequence (Hz)");
ylabel("Amplitude")
legend("Filtre passe-bas","DSP Normalisée")

%% 4.2 Retour en bande de base













