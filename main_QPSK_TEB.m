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

eb_n0_dB = 0:0.5:10; % rapport signal sur bruit (Eb = energie par bits, N0 amplitude du bruit)
eb_n0 = 10.^( eb_n0_dB /10) ; % Liste des Eb/N0
sigma2 = 0 ; % Variance du bruit complexe en bande de base

TEB = zeros ( size ( eb_n0 ) ); % Tableau des TEB ( resultats )
TEB2 = zeros ( size ( eb_n0 ) ); % Tableau des TEB ( resultats )
Pb = qfunc ( sqrt (2*eb_n0 ) ) ; % Tableau des probabilites d'erreurs theoriques

%% émetteur
sb = randi([0, M-1], Ns,1); % generation de nombre aleatoire entre 0 et M-1 (periode Te)
sb_bit = de2bi(sb); % on genere la tramme binaire associe a la tramme de symboles

%% association bit/symbole
ss = pskmod(sb, M, pi/M, 'gray'); % modulation en QPSK en mapping de gray
sigA2 = var(ss); % Variance des symboles
ss_surech = upsample(ss, Fse); % surechantillonage a Te

%% filtre de mise en forme
g = ones(1, 10); % reponse impultionnelle du filtre de mise en forme (1 quand 0<t<Ts) FSE = Ts/Te => nbr de coef a 1
T_filtre = size(g); T_filtre = T_filtre(2); % on recupere la taille du filtre de canal
Eg1 = sum(g.^2); % Energie du filtre de mise en forme
sl = conv(ss_surech, g); % application du filtre au signal

%% filtre de canal
hl(1) = 1; % hl est un dirac 
y = conv(sl, hl); % application du filtre du canal
        
%% itération de la chaine de transmission
for i = 1: length ( eb_n0 )
    error_cnt = 0;
    bit_cnt = 0;
    while error_cnt < 100
        
        %% bruit
        y_size = size(y);
        sigma2 = sigA2 * Eg1 ./ ( n_b * eb_n0 ) ; % Variance du bruit complexe en bande de base
        nl = sqrt ( sigma2(i) /2) * ( randn ( size (sl) ) + 1j* randn ( size (sl) ) ) ; % Generation du bruit blanc gaussien complexe
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
        
        error_cnt = error_cnt + biterr(s_out_bit(:), sb_bit(:)); % on incremente le nombre d'erreur
        bit_cnt = bit_cnt + Ns * n_b; % on compte le nombre de bits
        
    end
    TEB (i) = error_cnt / bit_cnt  % on calcule le TEB
end

%% afichage de la partie réel en temporel de QPSK avec filtre 111111
figure('Name', 'sl(t) et rl(t)');
n = linspace(0, 20*Fse-1, 20*Fse); % On genere les echantillons temporels
plot(n*Te, real(sl(1:20*Fse)), n*Te, real(rl(1:20*Fse))); % on trace la partie reel de sl(t) et rl(t) entre 0 et 10Ts-Te
title('partie reel de sl(t) et rl(t) pour t de 0 à 20Ts-Te'); 
xlabel('temps en secondes');
ylabel('amplitude');
legend({'sl(t)','rl(t)'},'Location','southwest');

[pw1, fpw1] = pwelch(sl,Nfft,0,Nfft,fe,'centered'); % diagramme de welch avec Nfft points et pas de recouvrement entre deux fenetres

figure('Name', 'évolution du TEB en fonction de Eb/N0 en dB');
plot(eb_n0_dB, 10*log10(TEB), eb_n0_dB, 10*log10(Pb));
ylabel('TEB (dB)');
xlabel('Eb/N0 en dB');
title('évolution du TEB en fonction de Eb/N0 en dB');
legend({'TEB avec g fenêtre','TEB théorique'},'Location','southwest');










%% avec un filtre de canal en cosinus sur eleve

%% filtre de mise en forme
g = rcosdesign(0.5, 8, Fse); % cosinus sur eleve roll-off 0.5 et temps de propagation 4Ts
T_filtre = size(g); T_filtre = T_filtre(2); % on recupere la taille du filtre de canal
Eg2 = sum(g.^2); % Energie du filtre de mise en forme
sl = conv(ss_surech, g); % application du filtre au signal


%% filtre de canal
hl(1) = 1; % hl est un dirac 
y = conv(sl, hl); % application du filtre du canal

for i = 1: length ( eb_n0 )
    error_cnt = 0;
    bit_cnt = 0;
    while error_cnt < 100
        
        %% bruit
        y_size = size(y);
        sigma2 = sigA2 * Eg2 ./ ( n_b * eb_n0 ) ; % Variance du bruit complexe en bande de base
        nl = sqrt ( sigma2(i) /2) * ( randn ( size (sl) ) + 1j* randn ( size (sl) ) ) ; % Generation du bruit blanc gaussien complexe
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
        
        error_cnt = error_cnt + biterr(s_out_bit(:), sb_bit(:)); % on incremente le nombre d'erreur
        bit_cnt = bit_cnt + Ns * n_b; % on compte le nombre de bits
        
    end
    TEB2 (i) = error_cnt / bit_cnt  % on calcule le TEB2
end

%% afichage de la partie réel en temporel de QPSK avec cos
figure('Name', 'sl(t) et rl(t)');
n = linspace(0, 20*Fse-1, 20*Fse); % On genere les echantillons temporels
plot(n*Te, real(sl(1:20*Fse)), n*Te, real(rl(1:20*Fse))); % on trace la partie reel de sl(t) et rl(t) entre 0 et 10Ts-Te
title('partie reel de sl(t) et rl(t) pour t de 0 à 20Ts-Te'); 
xlabel('temps en secondes');
ylabel('amplitude');
legend({'sl(t)','rl(t)'},'Location','southwest');

%% 3.3.2.1

figure('Name', 'évolution du TEB en fonction de Eb/N0 en dB');
plot(eb_n0_dB, 10*log10(TEB), eb_n0_dB, 10*log10(Pb), eb_n0_dB, 10*log10(TEB2));
ylabel('TEB (dB)');
xlabel('Eb/N0 en dB');
title('évolution du TEB en fonction de Eb/N0 en dB');
legend({'TEB avec g fenêtre','TEB théorique','TEB avec g en cosinus sur élevé'},'Location','southwest');



[pw2, fpw2] = pwelch(sl,Nfft,0,Nfft,fe,'centered'); % diagramme de welch avec Nfft points et pas de recouvrement entre deux fenetres
figure('Name', 'Périodogrammes de Welch');
plot(fpw1,10*log10(pw1),fpw2,10*log10(pw2*(Eg1/Eg2)),'g');
xlabel('fréquence en Hz');
ylabel('dB / (rad/sample)');
legend({'Avec g = fenêtre','Avec g = cosinus surélevé (puissance normalisée)'},'Location','southwest');
title('Périodogrammes de Welch');



