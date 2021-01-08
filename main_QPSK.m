%% GUENEGO CRESCENCE

clc;  % Efface la console
close all;  % Ferme les figures ouvertes
clear;  % Efface les variables de l'environnement de travail

fe = 10000; % Frequence d'echantillonnage
Te = 1/fe; % periode d'echantillonnage
Ds = 1000; % Debit symboles
Ts = 1/Ds; % periode symboles
M = 4; % Nombre de symboles dans la modulation
n_b = log2 (M) ; % Nombre de bits par symboles
Ns= 5000; % nombre de symboles a emettre par paquet
Fse = Ts/Te; % Facteur de sur echantillonnage
Nfft = 512; % nombre de points des transformees de Fourier


%% émetteur
sb = randi([0, M-1], Ns,1); % generation de nombre aleatoire entre 0 et M-1 (periode Te)
sb_bit = de2bi(sb); % on genere la tramme binaire associe a la tramme de symboles

%% association bit/symbole
ss = pskmod(sb, M, pi/M, 'gray'); % modulation en QPSK en mapping de gray
ss_surech = upsample(ss, Fse); % surechantillonage a Te

%% filtre de mise en forme
g = ones(1, Fse); % reponse impultionnelle du filtre de mise en forme (1 quand 0<t<Ts) FSE = Ts/Te => nbr de coef a 1
T_filtre = size(g); T_filtre = T_filtre(2); % on recupere la taille du filtre de canal
sl = conv(ss_surech, g); % application du filtre au signal

%% filtre de canal
hl(1) = 1; % hl est un dirac 
y = conv(sl, hl); % application du filtre du canal

%% bruit
y_size = size(y);
nl = zeros(1, y_size(2)); % génération du bruit
yl = y + nl; % ajout du bruit

%% Recepteur
% filtre de reception
ga = conj(flip(g)); %ga est l'inverse du filtre g;
rl = conv(yl, ga); % application du filtre de reception

%% Echantillonage à Ts
rl_n = rl(T_filtre:Fse:(Ns-1)*Fse+T_filtre); % sous echantillonage, en prenant en compte le retard du aux filtres du canal et de reception

%% association symbole/bit
s_out = pskdemod(rl_n, M, pi/M, 'gray'); % on demodule le signal echantillone rl_n en qpsk (M = 4), on specifie la phase a l'origine pi/M (dans la constelation) et le codage gray
s_out_bit = de2bi(s_out); % On genere la trame binaire associee a s_out


%% Affichage des resultats

TEB = biterr(s_out_bit(:), sb_bit(:)) / (Ns * n_b)  % TEB
scatterplot(yl); % affiche la constellation du signal module

%% 3.2.2 - 1
figure('Name', 'sl(t) et rl(t)');
n = linspace(0, 10*Fse-1, 10*Fse); % On genere les echantillons temporels
plot(n*Te, real(sl(1:10*Fse)), n*Te, real(rl(1:10*Fse))); % on trace la partie reel de sl(t) et rl(t) entre 0 et 10Ts-Te
title('partie reel de sl(t) et rl(t) pour t de 0 à 10Ts-Te'); 
xlabel('temps en secondes');
ylabel('amplitude');
legend({'sl(t)','rl(t)'},'Location','southwest');

%% 3.2.2 - 2

[pw, fpw] = pwelch(sl,Nfft,0,Nfft,fe,'centered'); % diagramme de welch avec Nfft points et pas de recouvrement entre deux fenetres

x1 = 1/Te*(-0.49:1/(Nfft-1):0.49); % on genere la frequence de la fft
dsp = Ts/Te*(Fse*Ts*sinc(x1*Ts)).^2; % On affiche la DSP theorique (TF carree d'une porte de largeur Ts) on ne rajoute pas le pi car dans matlab le sinus cardinal est defini : sin(pi t) / (pi t)


figure('Name', 'Périodogramme de Welch et DSP théorique de sl(t)');
plot(fpw,10*log10(pw),x1, 10*log10(dsp),'g');
xlabel('fréquence en Hz');
ylabel('dB / (rad/sample)');
legend({'Périodogramme de Welch','DSP théorique'},'Location','southwest');
title('Périodogramme de Welch et DSP théorique de sl(t)');



%% ajout manuel d'erreur
% if s_out_bit (1,2) == 0
%     s_out_bit (1,2) = 1; 
% else
%     s_out_bit (1,2) = 0; 
% end
% if s_out_bit (1,1) == 0
%     s_out_bit (1,1) = 1; 
% else
%     s_out_bit (1,1) = 0; 
% end




