En algebre on peut définir plusieurs types d'objets : les rationnels a/b, les polynomes P(X) = a + b*X +c*X^2 + ..., ou bien les matrices A ={ {a00,a01},{a10,a11}}. Chacun de ces objets obéit à des règles d'addition de muliplications, voire de division. Ces objets peuvent etre définis au-dessus d'autres objets, par exemple on peut faire des rationnels d'entier, mais aussi des rationnels de polynomes. Donc on a 3 types d'objets, qui se combinent différemment les uns avec les autres. Le but était d'aller le plus loin possible dans ces combinaisons.

Il y a des règles dans les combinaisons. Par exemple pour les rationnels d'entier, comme il y a une division avec reste pour les entiers, il y a l'algorithme du PGCD, et donc on peut simplifier les entiers. Idem pour les rationnels de polynomes de type T, si le type T a une division exacte on peut simplifier ces rationnels de polynome. Donc quand on écrit une class rationnel<T>, en fonction du type de T il y a des choses qu'on peut ou qu'on ne peut pas faire. J'ai fait donc des algorithmes de rationnel<T>, qui varient en fonction des algorithmes disponibles pour T.

Il y a donc une récurrence sur T lorsque j'écris une class rationnel<T>. J'ai géré cette récurrence. Pour exemple d'algorithme, on peut prendre a/b + c/d = (ad + bc)/cd. Donc quand on a deux objets (a,b) et (c,d), on obtient l'objet (ad+bc, cd). Ceci utilise l'addition et la multiplication pour le type de a. Ce genre d'algorithme utilise la multiplication et l'addition reliées au type de a (qui est le même que le type de b, de c, et de d). On écrit donc le même algorithhme, qui utilisera les algorithmes des oppérations du type a.

Les objets créés sont :
- les rationnels.
- les polynomes.
- les polynomes à n variables (en itératif et en récursif).
- les matrices.
- j'ai utilisé une class d'entiers infinis, et de réels à précision modifiable. Téléchargés sur internet.
- une class erreur qui contient un réel et l'erreur accumulée sur les opérations faites sur les réels.
- la diagonalisation de matrice, où on a les vecteurs propres en fonction de polynomes en lambda (les valeurs propres), sans connaitre les valeurs de lambda.
- les anneaux quotients : contient un élément, calculé modulo le quotient. par exemple les entiers modulo un nombre premier, ou un polynome modulo un autre polynome. Nécessite la division avec reste.
- les corps quotients : idem, mais avec une division exacte (calculée suivant le théorème de Bezout).
- les fichiers de récurrences sur les types T : si T est de type exacte ou approché (par exemple rationnel d'entier, ou réel). Ou si T a une division exacte (par exemple un rationnel d'autre chose), ou une division approchée (comme les entiers, ou les polynome<T> où T a une division exacte), ou sans aucune division.
- un fichier "unité" : donne l'élément nul ou l'élément 1, pour n'importe lequel des types au-dessus.
- un fichier "norme", pour savoir si l'objet est proche de 0. Par exemple pour minimiser les erreurs faites lors de l'inverse ou de la diagonalisation d'une matrice.
- des tests préparés dans le fichier "algebre.cpp", où se trouve le main.

Une partie du code prépare les algebres possibles (pré-compilation), et une autre exécute le code et fait les calculs. Une partie donc s'exécute lors de la compilation. Regarder le fichier "types.hpp".

Bonne lecture !
