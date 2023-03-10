En algebre on peut définir plusieurs types d'objets : les rationnels a/b, les polynomes P(X) = a + b*X +c*X^2 + ..., ou bien les matrices A ={ {a00,a01},{a10,a11}}. Chacun de ces objets obéit à des règles d'addition de muliplications, voire de division. Ces objets peuvent etre définis au-dessus d'autres objets, par exemple on peut faire des rationnels d'entier, mais aussi des rationnels de polynomes. Donc on a 3 types d'objets, qui se combinent différemment les uns avec les autres. Le but était d'aller le plus loin possible dans ces combinaisons.

Il y a des règles dans les combinaisons. Par exemple pour les rationnels d'entier, comme il y a une division avec reste pour les entiers, il y a l'algorithme du PGCD, et donc on peut simplifier les entiers. Idem pour les rationnels de polynomes de type T, si le type T a une division exacte on peut simplifier ces rationnels de polynome. Donc quand on écrit une classe rationnel<T>, en fonction du type de T il y a des choses qu'on peut ou qu'on ne peut pas faire. J'ai fait donc des algorithmes de rationnel<T>, qui varient en fonction des algorithmes disponibles pour T.

Il y a donc une récurrence sur T lorsque j'écris une classe rationnel<T>. J'ai géré cette récurrence. Pour exemple d'algorithme, on peut prendre a/b + c/d = (ad + bc)/cd. Donc quand on a deux objets a/b et c/d, on obtient l'objet (ad+bc)/(cd). Dans le programme, un rationnel consiste en la collection de deux objets (a,b), muni des opérations, qui décrit ce qu'on note en mathématiques a/b. Ceci utilise l'addition et la multiplication pour le type de a. Ce genre d'algorithme utilise la multiplication et l'addition reliées au type de a (qui est le même que le type de b, de c, et de d). On écrit donc le même algorithhme, qui utilisera les algorithmes des opérations du type a.

Les objets créés et fonctions implémentées sont :
- les rationnels<T>. Avec multiplication addition soustraction et division. Simplification si T a un algorithme de PGCD.
- les polynomes<T>. Avec addition et multiplication, et aussi PGCD et division avec reste dans le cas où T a une division exacte
- les polynomes<T> à n variables (deux versions : en itératif et en récursif). Multiplication et addition. Simplification P/Q si la fraction est exactement un polynome, dans le cas où T a une division exacte. Evaluation aussi : si on a P(X,Y,Z), on peut remplacer Z par un nombre, ou par un polynome Q(X,Y).
- les matrices<T>. Avec addition et multiplication, inverse de A et résolution AX=Y et noyau(A) si T a une division exacte. Determinant lent dans tous les cas, déterminant rapide si T a une division exacte, ou on peut calculer  déterminant<rationnel<T>> puis le simplifier si T a un PGCD.
- j'ai utilisé une classe d'entiers infinis, et de réels à précision modifiable. Trouvés sur internet.
- une classe erreur qui contient un réel et l'erreur accumulée sur les opérations faites sur les réels. Supporte l'addition la soustraction la multiplication et la division, et met à jour l'erreur cumulée. Est nul si la valeur est plus petite que l'erreur.
- la diagonalisation de matrice, où on a les vecteurs propres en fonction de polynomes en lambda (les valeurs propres), sans connaitre les valeurs de lambda.
- les anneaux quotients : contient un élément, calculé modulo le quotient. par exemple les entiers modulo un nombre premier, ou un polynome modulo un autre polynome. Nécessite la division avec reste.
- les corps quotients : idem, mais avec une division exacte (calculée suivant le théorème de Bezout).
- les fichiers de récurrences sur les types T : si T est de type exacte ou approché (par exemple rationnel d'entier, ou réel). Ou si T a une division exacte (par exemple un rationnel d'autre chose), ou une division approchée (comme les entiers, ou les polynome<T> où T a une division exacte), ou sans aucune division. Par exemple polynome<T> a une division approchée si T a une division exacte.
- un fichier "unité" : donne l'élément nul ou l'élément 1, pour n'importe lequel des types au-dessus. Prend un élément en parametre, car par exemple pour les polynomes à n variables il faut garder l'information du nombre de variables.
- un fichier "norme", pour savoir si l'objet est proche de 0. Par exemple pour minimiser les erreurs faites lors de l'inverse ou de la diagonalisation d'une matrice.
- des tests préparés dans le fichier "algebre.cpp", où se trouve le main.

Une partie du code prépare les algebres possibles (pré-compilation), et une autre exécute le code et fait les calculs. Une partie donc s'exécute lors de la compilation. Regarder le fichier "types.hpp".

Le but était que chaque objet peut être réutilisé dans la définition d'un nouvel objet "au-dessus" du premier. Tout ceci utilise le concept de template en C++. C'est-à-dire créer une classe qui prend une autre classe comme parametre. Cette implémentation se fait au moment de la compilation, et non au moment de l'utilisation du programme. Il faut donc tout préparer en amont.

C'était un exercice pour découvrir la notion de template, en essayant d'aller le plus loin possible dans ce projet.

Par exemple on peut noter les matrices de rationnels de polynômes à coefficients des rationnels de polynômes à n variables à coefficient rationnel d'entier (ou d'erreur de float_precision, pour le dernier). Le niveau maximal d'utilisation ! Ici pour calculer le polynôme caractéristique d'une matrice de type {{a,b},{c,d}}, ou du même genre la diagonalisation. (7 imbrications). a et b seraient des rationnels de polynomes à 4 variables, avec des coefficients rationnels ...

Bonne lecture !
