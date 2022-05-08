function [X w] = ParticleResampling(X,w,Ratio)
% ParticuleResampling   Pour effectuer le reechantillonnage des particules,
%                       si necessaire
% [XResampled WResampled] = ParticuleResampling(X,w,Ratio)
%   Entree
%        X : matrice de DxN, ou N est le nombre de particules et D le
%            nombre de variable d'etat
%        w : matrice de 1xN, ou N est le nombre de particules
%    Ratio : Ratio effectif en-deca duquel on reechantillonne. 
%
%  Sortie : X et w reechantillonnes
%  Ver 1.1

    % Normaliser les poids w;
    nParticules = size(w,2);
    Wnorm = sum(w);
    w = w./Wnorm;
    Copy = zeros(1,nParticules);
    
    % Resampling, pour combattre l'appauvrissement
    Neff = 1/sum(w.^2);
    
    if (isnan(Neff))
        Neff = 0;
        w = ones(1,nParticules)./nParticules;
    end

   
    % Verification si appauvrissement des particules
    if (Neff < Ratio*nParticules)
        % Effectuer un resampling.        
        display('Resampling');
        Q = cumsum(w);
        T = sort(rand(1,nParticules+1));
        T(nParticules+1) = 1;
        index = 1; jindex = 1;
        while (index <= nParticules)
            if (T(index) < Q(jindex))
                Copy(index) = jindex;
                index = index + 1;
            else
                jindex = jindex + 1;
            end
        end
        % Copie des particules, selon leur poids
        w = ones(1,nParticules)./nParticules;
        X = X(:,Copy);
    end
end