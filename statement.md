# Bonjour!

Ce programme est une simulation de l'épidémie de covid19 sur base d'un modèle différentiel SIR modifié.
Deux méthode d'intégration sont proposée :
* la méthode d'Euler (la plus simple)
* la méthode de Runge Kutta d'ordre 4 (la plus performante)

Cette simulation présente un scénario plutôt pessimiste, pour lequel le déconfinement provoque une deuxième vague suite à une remontée du taux de reproduction R au dessus de 1.  

@[Covid 19]({"stubs": ["test.py"], "command": "python3 test.py"})


# Exercice

Comparer les résultats des deux méthodes d'intégration avec différentes valeurs de subdivisions.
