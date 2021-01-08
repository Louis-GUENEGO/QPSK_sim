%% GUENEGO CRESCENCE

clear ; % Efface les variables de l'environnement de travail
close all ; % Ferme les figures ouvertes
clc ; % Efface la console

%% Initialisation des parametres
fe = 10000; % Frequence d'echantillonnage
Te = 1/fe; % periode d'echantillonnage
Ds = 1000; % Debit symboles
Ts = 1/Ds; % periode symboles
M = 4; % Nombre de symboles dans la modulation
n_b = log2 (M) ; % Nombre de bits par symboles
Ns= 5000; % nombre de symboles a emettre par paquet
Fse = Ts/Te; % Facteur de sur echantillonnage
Nfft = 512; % nombre de points des transformees de Fourier
f0 = 2500; % fréquence de la porteuse

eb_n0_dB = 0:0.5:10; % rapport signal sur bruit (Eb = energie par bits, N0 amplitude du bruit)
eb_n0 = 10.^( eb_n0_dB /10) ; % Liste des Eb/N0
sigma2 = 0 ; % Variance du bruit complexe en bande de base

TEB = zeros ( size ( eb_n0 ) ); % Tableau des TEB ( resultats )
Pb = qfunc ( sqrt (2*eb_n0 ) ) ; % Tableau des probabilites d'erreurs theoriques     %qfunc ( sqrt (1/2* eb_n0 ) ) ;

%% émetteur
sb = randi([0, M-1], Ns,1); % generation de nombre aleatoire entre 0 et M-1 (periode Te)
sb_bit = de2bi(sb); % on genere la tramme binaire associe a la tramme de symboles

%% association bit/symbole
ss = pskmod(sb, M, pi/M, 'gray'); % modulation en QPSK en mapping de gray
sigA2 = var(ss); % Variance des symboles
ss_surech = upsample(ss, Fse); % surechantillonage a Te

%% filtre de mise en forme
g = rcosdesign(0.5, 8, Fse); % cosinus sur eleve roll-off 0.5 et temps de propagation 4Ts
T_filtre = size(g); T_filtre = T_filtre(2); % on recupere la taille du filtre de canal
Eg = sum(g.^2); % Energie du filtre de mise en forme
sl = conv(ss_surech, g); % application du filtre au signal

%% Modulation
t = (0:1:length(sl)-1)';
TF = exp(j*2*pi*f0/fe*t);
s = real(sl.*TF);

%% filtre de canal
hl(1) = 1; % hl est un dirac 
y = conv(s, hl); % application du filtre du canal
yl = y;

%% Recepteur
% démodulation
fcos = 2*cos(2*pi*f0/fe*t);
fsin = 2*cos(2*pi*f0/fe*t+pi/2);
yI = yl .* fcos;
yQ = yl .* fsin;
yl_demod = yI + j*yQ;

% filtre de reception
ga = conj(flip(g)); %ga est l'inverse du filtre g; 
rl = conv(yl_demod, ga); % application du filtre de reception

%% Echantillonage à Ts
rl_n = rl(T_filtre:Fse:(Ns-1)*Fse+T_filtre); % sous echantillonage, en prenant en compte le retard du aux filtres du canal et de reception

%% association symbole/bit
s_out = pskdemod(rl_n, M, pi/M, 'gray'); % on demodule le signal echantillone rl_n en qpsk (M = 4), on specifie la phase a l'origine pi/M (dans la constelation) et le codage gray
s_out_bit = de2bi(s_out); % On genere la trame binaire associee a s_out

%% constelations
scatterplot(ss);
scatterplot(rl_n);

%% densite spectrale de puissance
[pw, fpw] = pwelch(s,Nfft,0,Nfft,fe,'centered'); % diagramme de welch avec Nfft points et pas de recouvrement entre deux fenetres
v_dsp = (-(length(g)-1)/2:1:(length(g)-1)/2);
dsp = ((fftshift(abs(fft(real(g.*exp(j*2*pi*f0/fe*v_dsp)),Nfft)))/300).^2)'; % On affiche la DSP theorique (TF carree du filtre)
figure('Name', 'Périodogramme de Welch et DSP théorique de s(t)');
fpw2 = (-5000:(10000/512):5000-1)'; 
plot(fpw,10*log10(pw),fpw2, 10*log10(dsp),'g');
xlabel('fréquence en Hz');
ylabel('dB / (rad/sample)');
legend({'Périodogramme de Welch','DSP théorique'},'Location','southwest');
title('Périodogramme de Welch et DSP théorique de s(t)');

%% Avec bruit

for i = 1: length ( eb_n0 )
    error_cnt = 0;
    bit_cnt = 0;
    while error_cnt < 100
        %% bruit
        y_size = size(y);
        sigma2 = sigA2 * Eg ./ ( n_b * eb_n0 ) ; % Variance du bruit complexe en bande de base
        nl = sqrt ( sigma2(i) /2) * ( randn ( size (sl) ) + 1j* randn ( size (sl) ) ) ; % Generation du bruit blanc gaussien complexe
        yl = y + nl; % ajout du bruit

        %% Recepteur
        % démodulation
        fcos = 2*cos(2*pi*f0/fe*t);
        fsin = 2*cos(2*pi*f0/fe*t+pi/2);
        yI = yl .* fcos;
        yQ = yl .* fsin;
        yl_demod = yI + j*yQ;
        
        % filtre de reception
        ga = conj(flip(g)); %ga est l'inverse du filtre g; 
        rl = conv(yl_demod, ga); % application du filtre de reception

        %% Echantillonage à Ts
        rl_n = rl(T_filtre:Fse:(Ns-1)*Fse+T_filtre); % sous echantillonage, en prenant en compte le retard du aux filtres du canal et de reception

        %% association symbole/bit
        s_out = pskdemod(rl_n, M, pi/M, 'gray'); % on demodule le signal echantillone rl_n en qpsk (M = 4), on specifie la phase a l'origine pi/M (dans la constelation) et le codage gray
        s_out_bit = de2bi(s_out); % On genere la trame binaire associee a s_out
        
        error_cnt = error_cnt + biterr(s_out_bit(:), sb_bit(:)); % on incremente le nombre d'erreur
        bit_cnt = bit_cnt + Ns * n_b; % on compte le nombre de bits
    end
    TEB (i) = error_cnt / bit_cnt  % on calcule le TEB2
end


%% TEB

figure('Name', 'évolution du TEB en fonction de Eb/N0 en dB');
plot(eb_n0_dB, 10*log10(TEB), eb_n0_dB, 10*log10(Pb));
ylabel('TEB (dB)');
xlabel('Eb/N0 en dB');
title('évolution du TEB en fonction de Eb/N0 en dB');
legend({'TEB','TEB théorique'},'Location','southwest');







%% 4.1.4
eb_n0_dB = 10; % rapport signal sur bruit (Eb = energie par bits, N0 amplitude du bruit)
eb_n0 = 10.^( eb_n0_dB /10) ; % Liste des Eb/N0

TEB = zeros ( size ( eb_n0 ) ); % Tableau des TEB ( resultats )
Pb = qfunc ( sqrt (2*eb_n0 ) ) ; % Tableau des probabilites d'erreurs theoriques     %qfunc ( sqrt (1/2* eb_n0 ) ) ;

%% émetteur
sb = randi([0, M-1], Ns,1); % generation de nombre aleatoire entre 0 et M-1 (periode Te)
sb_bit = de2bi(sb); % on genere la tramme binaire associe a la tramme de symboles

%% association bit/symbole
ss = pskmod(sb, M, pi/M, 'gray'); % modulation en QPSK en mapping de gray
sigA2 = var(ss); % Variance des symboles
ss_surech = upsample(ss, Fse); % surechantillonage a Te

%% filtre de mise en forme
g = rcosdesign(0.5, 8, Fse); % cosinus sur eleve roll-off 0.5 et temps de propagation 4Ts
T_filtre = size(g); T_filtre = T_filtre(2); % on recupere la taille du filtre de canal
Eg = sum(g.^2); % Energie du filtre de mise en forme
sl = conv(ss_surech, g); % application du filtre au signal

%% Modulation
t = (0:1:length(sl)-1)';
TF = exp(j*2*pi*f0/fe*t);
s = real(sl.*TF);

%% filtre de canal
hl(1) = 1; % hl est un dirac 
y = conv(s, hl); % application du filtre du canal

%% bruit
y_size = size(y);
sigma2 = sigA2 * Eg ./ ( n_b * eb_n0 ) ; % Variance du bruit complexe en bande de base
nl = sqrt ( sigma2 /2) * ( randn ( size (sl) ) + 1j* randn ( size (sl) ) ) ; % Generation du bruit blanc gaussien complexe
yl = y + nl; % ajout du bruit

%% Recepteur
% démodulation
fcos = 2*cos(2*pi*f0/fe*t+pi/6);
fsin = 2*cos(2*pi*f0/fe*t+pi/2+pi/6); %
yI = yl .* fcos;
yQ = yl .* fsin;
yl_demod = yI + j*yQ;

% filtre de reception
ga = conj(flip(g)); %ga est l'inverse du filtre g; 
rl = conv(yl_demod, ga); % application du filtre de reception

%% Echantillonage à Ts
rl_n = rl(T_filtre:Fse:(Ns-1)*Fse+T_filtre); % sous echantillonage, en prenant en compte le retard du aux filtres du canal et de reception

%% association symbole/bit
s_out = pskdemod(rl_n, M, pi/M, 'gray'); % on demodule le signal echantillone rl_n en qpsk (M = 4), on specifie la phase a l'origine pi/M (dans la constelation) et le codage gray
s_out_bit = de2bi(s_out); % On genere la trame binaire associee a s_out

%% constelations
scatterplot(rl_n);



%% 4.1.5
eb_n0_dB = 10; % rapport signal sur bruit (Eb = energie par bits, N0 amplitude du bruit)
eb_n0 = 10.^( eb_n0_dB /10) ; % Liste des Eb/N0

TEB = zeros ( size ( eb_n0 ) ); % Tableau des TEB ( resultats )
Pb = qfunc ( sqrt (2*eb_n0 ) ) ; % Tableau des probabilites d'erreurs theoriques     %qfunc ( sqrt (1/2* eb_n0 ) ) ;

%% émetteur
sb = randi([0, M-1], Ns,1); % generation de nombre aleatoire entre 0 et M-1 (periode Te)
sb_bit = de2bi(sb); % on genere la tramme binaire associe a la tramme de symboles

%% association bit/symbole
ss = pskmod(sb, M, pi/M, 'gray'); % modulation en QPSK en mapping de gray
sigA2 = var(ss); % Variance des symboles
ss_surech = upsample(ss, Fse); % surechantillonage a Te

%% filtre de mise en forme
g = rcosdesign(0.5, 8, Fse); % cosinus sur eleve roll-off 0.5 et temps de propagation 4Ts
T_filtre = size(g); T_filtre = T_filtre(2); % on recupere la taille du filtre de canal
Eg = sum(g.^2); % Energie du filtre de mise en forme
sl = conv(ss_surech, g); % application du filtre au signal

%% Modulation
t = (0:1:length(sl)-1)';
TF = exp(j*2*pi*f0/fe*t);
s = real(sl.*TF);

%% filtre de canal
hl(1) = 1; % hl est un dirac 
y = conv(s, hl); % application du filtre du canal

%% bruit
y_size = size(y);
sigma2 = sigA2 * Eg ./ ( n_b * eb_n0 ) ; % Variance du bruit complexe en bande de base
nl = sqrt ( sigma2 /2) * ( randn ( size (sl) ) + 1j* randn ( size (sl) ) ) ; % Generation du bruit blanc gaussien complexe
yl = y + nl; % ajout du bruit

%% Recepteur
% démodulation
fcos = 2*cos(2*pi*(f0+10)/fe*t);
fsin = 2*cos(2*pi*(f0+10)/fe*t+pi/2);
yI = yl .* fcos;
yQ = yl .* fsin;
yl_demod = yI + j*yQ;

% filtre de reception
ga = conj(flip(g)); %ga est l'inverse du filtre g; 
rl = conv(yl_demod, ga); % application du filtre de reception

%% Echantillonage à Ts
rl_n = rl(T_filtre:Fse:(Ns-1)*Fse+T_filtre); % sous echantillonage, en prenant en compte le retard du aux filtres du canal et de reception

%% association symbole/bit
s_out = pskdemod(rl_n, M, pi/M, 'gray'); % on demodule le signal echantillone rl_n en qpsk (M = 4), on specifie la phase a l'origine pi/M (dans la constelation) et le codage gray
s_out_bit = de2bi(s_out); % On genere la trame binaire associee a s_out

%% constelations
scatterplot(rl_n);



%% 4.1.6
eb_n0_dB = 10; % rapport signal sur bruit (Eb = energie par bits, N0 amplitude du bruit)
eb_n0 = 10.^( eb_n0_dB /10) ; % Liste des Eb/N0

TEB = zeros ( size ( eb_n0 ) ); % Tableau des TEB ( resultats )
Pb = qfunc ( sqrt (2*eb_n0 ) ) ; % Tableau des probabilites d'erreurs theoriques     %qfunc ( sqrt (1/2* eb_n0 ) ) ;

%% émetteur
sb = randi([0, M-1], Ns,1); % generation de nombre aleatoire entre 0 et M-1 (periode Te)
sb_bit = de2bi(sb); % on genere la tramme binaire associe a la tramme de symboles

%% association bit/symbole
ss = pskmod(sb, M, pi/M, 'gray'); % modulation en QPSK en mapping de gray
sigA2 = var(ss); % Variance des symboles
ss_surech = upsample(ss, Fse); % surechantillonage a Te

%% filtre de mise en forme
g = rcosdesign(0.5, 8, Fse); % cosinus sur eleve roll-off 0.5 et temps de propagation 4Ts
T_filtre = size(g); T_filtre = T_filtre(2); % on recupere la taille du filtre de canal
Eg = sum(g.^2); % Energie du filtre de mise en forme
sl = conv(ss_surech, g); % application du filtre au signal

%% Modulation
t = (0:1:length(sl)-1)';
TF = exp(j*2*pi*f0/fe*t);
s = real(sl.*TF);

%% filtre de canal
hl(1) = 1; % hl est un dirac 
y = conv(s, hl); % application du filtre du canal

%% bruit
y_size = size(y);
sigma2 = sigA2 * Eg ./ ( n_b * eb_n0 ) ; % Variance du bruit complexe en bande de base
nl = sqrt ( sigma2 /2) * ( randn ( size (sl) ) + 1j* randn ( size (sl) ) ) ; % Generation du bruit blanc gaussien complexe
yl = y + nl; % ajout du bruit

%% Recepteur
% démodulation
fcos = 2*cos(2*pi*(f0+1000)/fe*t);
fsin = 2*cos(2*pi*(f0+1000)/fe*t+pi/2);
yI = yl .* fcos;
yQ = yl .* fsin;
yl_demod = yI + j*yQ;

% filtre de reception
ga = conj(flip(g)); %ga est l'inverse du filtre g; 
rl = conv(yl_demod, ga); % application du filtre de reception

%% Echantillonage à Ts
rl_n = rl(T_filtre:Fse:(Ns-1)*Fse+T_filtre); % sous echantillonage, en prenant en compte le retard du aux filtres du canal et de reception

%% association symbole/bit
s_out = pskdemod(rl_n, M, pi/M, 'gray'); % on demodule le signal echantillone rl_n en qpsk (M = 4), on specifie la phase a l'origine pi/M (dans la constelation) et le codage gray
s_out_bit = de2bi(s_out); % On genere la trame binaire associee a s_out

%% constelations
scatterplot(rl_n);