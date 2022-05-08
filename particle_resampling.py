import numpy as np
def particle_resampling(X, w, Ratio):
# ParticuleResampling   Pour effectuer le reechantillonnage des particules,
#                       si nécessaire
# [XResampled WResampled] = ParticuleResampling(X,w,Ratio)
#   Entree
#        X : matrice de DxN, où N est le nombre de particules et D le
#            nombre de variables d'état
#        w : matrice de 1xN, où N est le nombre de particules
#    Ratio : Ratio effectif en-deça duquel on réechantillonne.
#
#  Sorties : X et w réechantillonnés
#  Ver 1.1

    # Normaliser les poids w;
    nParticules = len(w)
    Wnorm = np.sum(w)
    w = w/Wnorm
    Copy = np.zeros(nParticules)

    # Resampling, pour combattre l'appauvrissement
    Neff = 1 / sum(w**2)

    if (np.isnan(Neff)):
        Neff = 0
        w = np.ones(nParticules)/nParticules

    # Verification si appauvrissement des particules
    if (Neff < Ratio*nParticules):
        # Effectuer un resampling.
        print('Resampling')
        Q = np.cumsum(w)
        T = np.sort(np.random.rand(nParticules+1))
        T[nParticules] = 1
        Copy = np.zeros(nParticules, dtype = np.uint)
        index = 0
        jindex = 0
        while (index < nParticules):
            if (T[index] < Q[jindex]):
                Copy[index] = jindex
                index = index + 1
            else:
                jindex = jindex + 1
        # Copie des particules, selon leur poids
        w = np.ones(nParticules)/nParticules
        X = X[:,Copy]
    return (X, w)

#X = np.random.randn(2,10)
#w = np.random.rand(10)
#(Yo, Yow) = particle_resampling(X,w,0.9)