clearvars
clf
myFontSize = 14;

% Parametre du systeme
Pwifi   =  9.0;      % Puissance d'?mission du routeur WiFi
d_init  =  5.0;      % Point de depart du robot
SRecepteur = 0.3;    % Bruit sur la puissance du r?cepteur radio, en mW.
dT      = 0.10;      % Intervalle de temps entre les mesures/pas
RatioBruitAccel = 0.2;

% Important! Il faut initialiser la matrice X de facon a avoir
% une dimension (2,1) et non (1,2). Sinon, le filtre EKF ne marchera pas.
xVrai = [d_init 0]'; % Etat reel, inconnu du filtre et du robot.

% Matrices de bruits
Cw = SRecepteur^2.0;  % Matrice de covariance sur les mesures

nStep = 400; % Nombre d'iterations

% Nombre de particules C pris dans [4 10 40];
nParticules = 40;

% Je vais initialiser au hasard les particules, dans une boite de 5m par 5m
% centree sur (0,0)
X = ones(1,nParticules);

% situation a)
% X(1,:) = xVrai(1)*ones(1,nParticules); % position en x
% X(2,:) = xVrai(2)*ones(1,nParticules); % vitesse lineaire

% situation b)
X(1,:) = 3 + 9*rand(1,nParticules);
X(2,:) = zeros(1,nParticules);

w = ones(1,nParticules)/nParticules; % Poids egaux pour les particules

% ============== Boucle du filtre a particules =================
for iStep = 1:nStep

    % ======= Debut de la simulation du systeme reel =============
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
    % ======= Fin de la simulation du systeme reel =============
    
    % Le bruit sur la commande est variable. Il faut donc le calculer à 
    % chaque itération
    Cv = 0.20 * U;
    
    % ================ Debut du filtre a particules =============    
    for iParticule = 1:nParticules
        % Je propage chaque particule
        X(1,iParticule) = X(1,iParticule) + X(2,iParticule)*dT + (acceleration+Cv*randn)/2*(dT^2);
        X(2,iParticule) = X(2,iParticule) + (acceleration+Cv*randn)*dT;
        
        % Mise-a-jour des poids
        zhat = Pwifi*exp(-X(1,iParticule)/8.0);
        wnew = gauss(zhat-z,SRecepteur);
        w(iParticule) = w(iParticule).*wnew;

    end
    [X, w] = ParticleResampling(X,w,0.5);
    % ================ Fin du filtre a particules =============
    
    % Cueillette des donnees pour les graphiques/statistiques
    AxVrai1(iStep)  = xVrai(1);
    AX1(iStep)      = sum(X(1,:).*w); % Position pondérée
    AZ(iStep)       = z; 
    ATime(iStep)    = time;
    VarPosition(iStep) = std(X(1,:));
    
    % Pour voir le filtre evoluer dans le temps
    clf
    h(1) = plot(ATime,AX1,'go');
    hold on;
    h(2) = plot(ATime,AxVrai1,'k-','LineWidth',2);
    h(3) = plot(ATime,-10*log(AZ/Pwifi),'r*'); % Ici on peut inverser le capteur, pour trouver la position correspondant a z.
    xlabel('Temps (s)');
    ylabel('Estime de position (m)');
    legend(h,{'Filtre à particules','Position Exacte','Mesure h_z^{-1}'});
    ylim([0 25]);
    drawnow();
end

fprintf('C = %d\n', nParticules);
fprintf('Précision = %.4f\n', mean(VarPosition));

if (nParticules == 40),
    clf
    subplot(3,1,1:2);
    h(1) = plot(ATime,AX1,'go');
    hold on;
    h(2) = plot(ATime,AxVrai1,'k-','LineWidth',2);
    h(3) = plot(ATime,-10*log(AZ/Pwifi),'r*'); % Ici on peut inverser le capteur, pour trouver la position correspondant a z.
    xlabel('Temps (s)');
    ylabel('Estime de position (m)');
    legend(h,{'Filtre à particules','Position Exacte','Mesure h_z^{-1}'});
    ylim([0 25]);
    text(3,15,sprintf('Erreur absolue moyenne = %.4f',mean(abs(AxVrai1-AX1))));

    subplot(3,1,3);
    h(1) = plot(ATime,VarPosition,'k-');
    xlabel('Temps (s)');
    ylabel('Variance');
end