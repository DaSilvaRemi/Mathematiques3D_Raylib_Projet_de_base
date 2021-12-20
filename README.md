# Mathématiques pour la 3D
  
  E3FI - 1l
  21_E3FI_3I_SI5
  Da Silva Rémi et Bailleul Valentin
  Décembre 2021


## Chemin relatif du projet
  \projects\VS2019\raylib.sln

### En cas d'erreur d'include
Si une erreur d'include de fonctions occure, reciblez les projets :
Pour LES TROIS DOSSIERS 'core_basic_window' , 'core_basic_window_cpp' et 'raylib' :
    - Clic droit, "Propriétés", "Propriétés de configuration", "Répertoires VC++", "Général", "Répertoire Include"
    - Selectionnez la flèche vers le bas qui s'affiche au survol de la section Répertoire Include
    - Cliquez sur "<Modifier...>",
    - Ouvrez l'explorateur de fichiers de Windows, allez dans "C:\Program Files (x86)\Windows Kits\10\Include"
    - Copiez le lien "C:\Program Files (x86)\Windows Kits\10\Include" (selon votre configuration, il peut être légèrement modifié)
    - Revenez sur la fenêtre pop-up Visual Studio
    - Cliquez sur l'icône jaune représentant un dossier en haut à droite de la fenêtre
    - Collez le lien C:\ , puis OK, Appliquer, OK
    - Ensuite, Clic droit sur les 3 dossiers, "Recibler les projets", puis sélectionnez la version qui se trouve dans votre explorateur d'où vous avez trouvé le lien "C:\Program Files (x86)\Windows Kits\10\Include"
    - Puis "OK"

Une fois ceci fait POUR LES TROIS DOSSIERS, faites F5 pour compiler et lancer le jeu

# Description technique

## Remarques particulières
Les limitations de notre projet :
* L’intersection avec la RoundedBox  ne fonctionne uniquement que si celle-ci n’est pas tournée par un quaternion
* Méthode d’intersection Segment-Box à revoir (non utilisée dans notre projet final, celle-ci n’est pas considérée comme un bug à part entière)
* Les machines sur lesquelles nous travaillons ont un rendu avec très peu de fps


## Bugs connus :
* Parfois, dans la scène de jeu, la balle traverse les Quads muraux si elle rentre en collision avec un autre objet dans la même frame


## Voies d’améliorations :
* Nous voudrions avoir une RoundedBox qui puisse avoir une collision adéquate lorsque celle-ci est tournée sur elle-même
* Ajouter des segments-détecteurs de collision, à la sphère, afin de détecter plus précisément les collisions arrivant à différents endroits de cette sphère


<hr />


# Répartition des tâches
## Rémi Da Silva : 
* Réalisation des méthodes de dessin
* Réalisation des méthodes d’intersection en collaboration avec Valentin
* Réalisation des méthodes de conversion de référentiels et conversion de coordonnées
* Création des structures mathématiques
* Structuration du code en plusieurs fichiers
* Débogage des méthodes de dessin
* Création de la scène de jeu


## Valentin Bailleul : 
* Réalisation de la V2 de la RoundedBox
* Ecriture de l’intersection avec la RoundedBox
* Résolution et simplification des divers formules de calculs d’intersections
* Réalisation des méthodes d’intersection en collaboration avec Rémi
* Débogage pointu des dessins, des intersections avec un segment
* Amélioration de la scène de jeu


## Courte vidéo du projet
[Lien google drive vers la vidéo](https://drive.google.com/file/d/1Tb73aJxSwJBbJnyB7VOgOq8hKEaePGPa/view?usp=sharing)