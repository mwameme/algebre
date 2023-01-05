utilisation : (et tout en bas l'explication)

Pour utiliser la class des polynome<T> sur une autre classe (la classe des coefficients), il y a besoin :
-que les coefficients aient les opérations + - * / (que la division soit parfaite : un corps).
-toujours la taille du vecteur coefficient doit être >=1 (au moins égale à 1)
-qu'il y ait l'opération =(true / false) pour mettre un coefficient valant 0 (l'élément neutre de l'addition) ou valant 1 (l'élément neutre de la multiplication)
les polynomes ont une multiplication, une addition, une soustraction, mais aussi une division (partielle) et son reste. On peut mettre polynome = true/false , auquel cas la taille des coefficients est 1, et le premier coefficient =true/false

Pour utiliser la class des matrice<T> sur une autre classe (la classe des coefficients), il y a besoin :
-les mêmes condition que pour les polynomes ...
les matrices peuvent être additionnées, soustraites, multipliées, mais il y a aussi une division (si determinant non-nul) exacte. On peut mettre une matrice = true/false (et alors la diagonale vaut = true/false , le reste = false).

Pour utiliser la class des rationnel sur une autre classe (la classe des coefficients), il y a besoin :
-de faire attention que le denominateur soit différent de 0
-qu'il y ait l'opération =(true / false) pour mettre un coefficient valant 0 (l'élément neutre de l'addition) ou valant 1 (l'élément neutre de la multiplication)
-qu'il y ait une division et un reste (de la division) partielle : utilisé pour réduire les fractions (diviser le numérateur et le dénominateur par le PGCD, calculé avec / et %)
les rationnels ont les opérations + - * et / (exacte). Le reste (%) est mis à 0 (la division est exacte). L'opératin =true/false met le numerateur = true/false, et le denominateur = true;

Quand on utilise les corps fini, attention que p soit un nombre premier, et le même des deux côtés (lors des opérations).
C'est un corps, donc il y a + - * /, et le reste (de la divisio) est mis à 0 (utilisé lors du pgcd).
le "= true/false" met x = 0 ou 1

Idem, le corps des nombres complexes est un corps ... vérifier avant division que le nombre est non-nul.

Toujours il y a la transformation d'un élément en booléen : si tous les coefficients sont nuls (polynome matrice corps-fini et complexe), ou si le numerateur est nul (fraction). C'est utile pour savoir si l'objet est nul ou non (pour les divisions par exemple).
Au lieu d'avoir un contructeur qui prend un booléen en parametre, j'utilise la copie d'un élément que je transforme avec = true/false. Ca permet de garder les caractéristiques lorsqu'il y a plusieurs éléments neutres (pour les matrices il faut aussi spécifier la taille ... le true/false ne donne pas cette information ... voilà pourquoi (et idem pour les corps finis)).

Au final toutes les opérations sont définies, avec la subtilité : la division est exacte (pour les corps), ou approchée (anneaux : polynomes ou entiers). Quand on a un anneau, on le transforme en corps avec rationnel<>
Il faut pouvoir mettre les éléments = true/false, et vérifier si un élément est nul (==? true/false)
Attention aux dépassements de la taille des entiers (j'utilise long long).


Explication :
partant d'une classe générique (notée T), si il y a de définis certaines opérations sur cette classe (+ - * / %), on peut générer une sur-classe qui prend T en paramètre. Par exemple, on peut définir des rationnel<T>. Comme la programmation est effectuée de manière générale (définir et utiliser les opérations sur ma_classe<T> , étant définies certaines opérations sur T) on peut composer ces classes. Par exemple pour calculer le polynome caractéristique d'une matrice<rationnel<entier>>, on utilise la classe rationne<polynome<rationnel<entier>>>. Comme il faut utiliser des opérations prédéfinies sur T, il y a la contrainte que ces opérations vérifient certaines conditions, et comme j'ai voulu itérer ces compositions de classes, il faut que les opérations sur ma_classe<T> vérifient aussi ces conditions.
Finalement il faut séparer en deux types de classes : les corps, avec une division exacte (b= (b/a) * a), et les anneaux (comme les polynomes ou les entiers) dont on peut faire des rationnels, qui eux sont des corps. Les polynomes et les matrices sont définies sur des corps (les coefficients doivent être des éléments de type corps). Normalement j'ai tout fait de manière composable, ce qui était le but de l'exercice. Vous pouvez créer vos propres anneaux (commutatifs je crois) et corps (idem), et vérifier que ça marche.