# geostatistique-r

# Qu'est que ça prend pour utiliser les outils géostatistiques?
- Échantillhonage spatial
- Nombre de points suffisant
- Plan d'échantillonage approprié

# Variogramme
## Régles
- Longueru maximale: 1/2 d'un cêté ou 1/2 sqrt(aire)
- Nombre de lags? s'assurer d'avoir assez de paires de pounts a chaque lag
- répresentation des petites distances
- traitement de l'anisotropie: structure de variograme differénte selon les directions
-- variogram map

# Autocorrélation spatiale
- C de Geary - lié au variograme
-- if the sample is not different from 1, there's autocorrelation at that point
- I de Moran - présence d'autocorrélation "globale" ou par lag
-- if the sample is not different from 0, there's autocorrelation at that lag

# Fonction aléatoire
- Valeur moyenne et une fonction d'autocovariance/modéle de variogramme
- Étendue infinie
- Structurée spatialement, mais avec un aspect stochastique/aléatoire

## Simulation Gaussian non-conditionelle
- Pas des donées!
- Des nuggets vont s'aditioner au modéle spatial

# Estimation and kriging
- Estimate values at unsampled locations ~> interpolation

## Deterministic functions



