% Systeme representant le chariot pour le TP2, Automne 2022
% (c) Philippe Giguere 2022
clearvars
clf
myFontSize = 14; % Taille de la police de caractere pour les graphes

% Parametre du systeme
Pwifi   =  9.0;      % Puissance d'émission du routeur WiFi
d_init  =  5.0;      % Point de depart du robot
SRecepteur = 0.3;    % Bruit sur la puissance du récepteur radio, en mW.
nStep   = 400;       % Nombres de mesures/pas
dT      = 0.10;      % Intervalle de temps entre les mesures/pas
RatioBruitAccel = 0.2;

% Important! Il faut initialiser la matrice X de facon a avoir
% une dimension (2,1) et non (1,2). Sinon, le filtre EKF ne marchera pas.
xVrai = [d_init 0]'; % Etat reel, inconnu du filtre et du robot.

% Specifier les valeurs initiales des matrices.
% Ne pas oublier qu'ici, ce sont des covariances, pas des ecarts-types.

% Situation a)
X = [d_init 0]';
P = [0 0; 0 0];

% Situation b)
% X = [d_init 0]';
% P = [25 0; 0 0];

% Situation c)
% X = [8 0]';
% P = [0 0; 0 0];

% Situation d)
% X = [8 0]';
% P = [100 0; 0 1];

% Question 4.3
% X = [80 0]';
% P = [0 0; 0 3];

% Matrices de bruits
Cw = SRecepteur^2.0;  % Matrice de covariance sur les mesures

for iStep = 1:nStep
    % Simulation du systeme a chaque etape
    time = iStep*dT;
    
    % Commande d'accélération envoyee vers le systeme.
    U = 0.10*cos(0.04*pi*time);

    % ============== Debut de la simulation du deplacement reel ===========   
    % Le deplacement veritable du chariot selon les equations.
    % Je vous donne les equations, vous n'avez rien a changer ici.
    acceleration = U + RatioBruitAccel*U*randn;
    
    xVrai(1) = xVrai(1) + xVrai(2)*dT + acceleration/2*(dT^2); % Calcul du deplacement
    xVrai(2) = xVrai(2) + acceleration*dT;   % Calcul de la vitesse
    
    % Je simule pour vous la réponse du récepteur radio du chariot.
    z = Pwifi*exp(-xVrai(1)/8) + SRecepteur*randn;
    % =============== Fin de la simulation de deplacement reel ============   

    % ================ Debut de votre filtre E K F ou particule ==================
    % ATTENTION ATTENTION ATTENTION ATTENTION
    % Vous n'avez pas le droit d'utiliser xVrai dans votre filtre
    % car c'est la position et la vitesse reele du systeme, et 
    % elles vous sont inconnues.

    % ========= Calcul des matrices Jacobiennes =============
    F = [1 dT; 0 1];
    G = [dT^2/2; dT];
    H = [-(1/8.0)*Pwifi*exp(-X(1,1)/8.0) 0];

    % Le bruit sur la commande est variable. Il faut donc le calculer à 
    % chaque itération
    Cv = 0.20 * U;

    % ======= Propagation =========
    X = [X(1,1) + X(2,1) * dT + U(1,1) * dT^2 / 2; X(2,1) + U(1,1) * dT];
    P = F*P*F' + G*Cv*G';

    % ======= Mise-a-jour ========
    K = P*H' / (H*P*H' + Cw);
    r = z - Pwifi*exp(-X(1,1)/8.0);
    X = X + K*r;
    P = (ones(2) - K*H) * P;

    % ========= Fin des equations du filtre EKF ou particule =============
        
    % Cueillette des donnees pour les graphiques/statistiques
    AxVrai1(iStep)  = xVrai(1);
    AxVrai2(iStep)  = xVrai(2);
    AX1(iStep)      = X(1);
    AX2(iStep)      = X(2);
    AU(iStep)       = U;
    AZ(iStep)       = z;
    VarPosition(iStep) = P(1,1);
    ATime(iStep) = time;
    
    % Erreurs 
    if (time < 2)
        err_debut = mean(abs(AxVrai1-AX1));
    end
    
    % Pour voir votre filtre evoluer dans le temps
    clf
    h(1) = plot(ATime,AX1,'go');
    hold on;
    h(2) = plot(ATime,AxVrai1,'k-','LineWidth',2);
    h(3) = plot(ATime,-10*log(AZ/Pwifi),'r*'); % Ici on peut inverser le capteur, pour trouver la position correspondant a z.
    xlabel('Temps (s)');
    ylabel('Estime de position (m)');
    legend(h,{'EKF','Position Exacte','Mesure h_z^{-1}'});
    % ylim([0 20]);
    drawnow();
end

clf
subplot(3,1,1:2);
h(1) = plot(ATime,AX1,'go');
hold on;
h(2) = plot(ATime,AxVrai1,'k-','LineWidth',2);
h(3) = plot(ATime,-10*log(AZ/Pwifi),'r*'); % Ici on peut inverser le capteur, pour trouver la position correspondant a z.
xlabel('Temps (s)');
ylabel('Estime de position (m)');
legend(h,{'EKF','Position Exacte','Mesure h_z^{-1}'});
% ylim([0 20]);
text(3,15,sprintf('Erreur abs moy pour t < 2s = %.4f',err_debut));
text(3,16,sprintf('Erreur absolue moyenne = %.4f',mean(abs(AxVrai1-AX1))));

subplot(3,1,3);
h(1) = plot(ATime,VarPosition,'k-');
xlabel('Temps (s)');
ylabel('Variance P(1,1)');


