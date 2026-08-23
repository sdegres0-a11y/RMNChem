# RMNChem — analyse de formule brute et spectre ¹H simulé

Add-in pour Casio Graph 90+E (fx-CG 50). On saisit les indices d'une formule
brute C\_a H\_b O\_c N\_d, l'add-in calcule le degré d'insaturation, cherche la
molécule dans une base embarquée, trace sa formule développée et affiche le
spectre RMN ¹H correspondant avec déplacements, intégrations et multiplicités.

Base embarquée : 46 molécules, dont 9 formules à isomères multiples
(C8H10, C2H6O, C3H8O, C4H10O, C3H6O, C4H8O, C2H4O2, C3H6O2, C2H5NO).
Quand la formule n'est pas tabulée, l'add-in bascule sur un mode secours qui
donne le DI exact, les fragments compatibles et l'échelle des zones de
déplacement.

---

## Sept points du cahier des charges à corriger

Avant de compiler, quelques précisions : plusieurs éléments de la demande
initiale ne correspondent pas à ce que fait réellement le fxSDK.

**1. La commande de compilation n'est pas `fxsdk build-fx`.**
Celle-là produit un `.g1a`, le format des fx-9860G monochromes.
Pour la fx-CG 50, donc la Graph 90+E, c'est `fxsdk build-cg`, qui produit
un `.g3a`.

**2. Le fichier `addin.json` n'existe pas.**
Le fxSDK n'a jamais utilisé ce nom. La configuration se répartit entre
`project.cfg` (métadonnées lues par l'outil en ligne de commande) et l'appel
`generate_g3a()` du `CMakeLists.txt`, qui est celui qui fixe réellement le nom
affiché et les icônes au moment de la fabrication du `.g3a`.
CMake est le système de build par défaut depuis la version 2.3.
Les deux fichiers sont fournis ici.

**3. L'icône ne fait pas 64×48 mais 92×64.**
Il en faut deux, un état non sélectionné et un état sélectionné. Elles sont
stockées non compressées en RGB565 dans l'en-tête du `.g3a`.

**4. La cible « moins de 50 Ko » est tenable, mais de justesse.**
L'en-tête `.g3a` occupe à lui seul 0x7000 octets, soit 28 672, dont
2 × 92 × 64 × 2 = 23 552 pour les icônes. Il reste donc environ 21 Ko pour le
code et les données. Les données chimiques pèsent ~3,5 Ko sur SH4 ; avec `-Os`
et `--gc-sections` (déjà dans le `CMakeLists.txt`) le reste passe sans
difficulté. Ne pas ajouter de calcul flottant, c'est ce qui ferait exploser le
budget en tirant la bibliothèque mathématique.

**5. L'écran ne fait pas 384×216 mais 396×224.**
Le 384×216 correspond à la zone que l'OS Casio laisse au Basic. Un add-in gint
prend la main sur tout l'écran. Le code utilise `DWIDTH` et `DHEIGHT`, il
s'adapte donc tout seul si tu le recompiles pour une autre machine.

**6. MicroPython n'a pas d'API clavier.**
Le module `casioplot` fourni par Casio expose le dessin, pas la lecture des
touches. La variante MicroPython livrée dans `micropython/` saisit donc les
indices au prompt avec `input()` au lieu des touches de fonction. C'est la
raison pour laquelle la version C reste la bonne : l'interface interactive
demandée n'est pas réalisable en Python sur cette machine.

**7. `dprint` ne sait pas formater les flottants.**
Tout est en entiers dans le code : les déplacements chimiques sont stockés en
dixièmes de ppm (72 pour 7,2) et convertis à l'affichage, virgule française
comprise.

---

## Installation de la chaîne de compilation

Sous Linux, ou sous WSL si tu es sur Windows. Compter une bonne heure la
première fois : il faut compiler un cross-compilateur complet.

Dépendances système d'abord.
Sur Debian, Ubuntu et WSL : `sudo apt install python3-pil libusb-dev`,
plus `udisks2` en option.
Sur Arch et Manjaro : `sudo pacman -S python-pillow libusb sdl2 ncurses patch`.

Ensuite GiteaPC, l'outil qui automatise le reste. Son script d'installation est
sur la forge de Planète Casio, à `Lephenixnoir/GiteaPC` — récupère la commande
depuis le README du dépôt plutôt que de recopier une URL qui peut avoir bougé.
Redémarre le terminal après, le script modifie le `PATH`.

Puis, dans cet ordre :

```sh
giteapc install Lephenixnoir/fxsdk Lephenixnoir/sh-elf-binutils \
                Lephenixnoir/sh-elf-gcc Lephenixnoir/sh-elf-gdb
giteapc install Lephenixnoir/OpenLibm Vhex-Kernel-Core/fxlibc
giteapc install Lephenixnoir/sh-elf-gcc
giteapc install Lephenixnoir/gint
```

Le double passage sur `sh-elf-gcc` n'est pas une erreur : au premier passage
la bibliothèque C++ ne peut pas être construite faute de bibliothèque C, et la
réinstallation reprend là où elle s'était arrêtée une fois OpenLibm et fxlibc
en place.

Deux cas particuliers.
Si `udisksctl` n'existe pas sur ta machine, ajoute `:noudisks2` après
`Lephenixnoir/fxsdk`. Et si binutils et GCC sont déjà installés dans une
version postérieure au fxSDK 2.9, `:any` après les deux évite une
recompilation d'une trentaine de minutes.

---

## Compilation

Le plus sûr est de laisser le fxSDK générer un projet neuf, puis d'y déposer
les fichiers d'ici. Ça garantit que `project.cfg` et l'environnement CMake sont
au bon format pour ta version installée.

```sh
fxsdk new RMNChem
cd RMNChem
```

Copie ensuite `src/main.c`, `src/chem_db.h`, `assets-cg/icon-uns.png`,
`assets-cg/icon-sel.png` et le `CMakeLists.txt` fourni par-dessus les fichiers
générés. Puis :

```sh
fxsdk build-cg
```

Si tu as des options de configuration à passer, lance une première fois avec
`-c` pour configurer, puis sans pour construire.
Le résultat est `RMNChem.g3a` à la racine du projet.

Vérifie la taille : `ls -l RMNChem.g3a`. Attends-toi à un fichier entre 32 et
40 Ko. Au-delà de 50, regarde d'abord si un `printf` avec `%f` ne s'est pas
glissé quelque part.

---

## Transfert vers la calculatrice

Deux méthodes.

**En clé USB.** Branche la Graph 90+E, elle propose un mode de connexion :
choisis la mémoire de stockage. Elle apparaît comme un volume amovible. Copie
`RMNChem.g3a` à la racine, pas dans un sous-dossier. Démonte proprement, puis
débranche. L'add-in apparaît dans le menu principal.

**Avec fxlink**, l'outil de communication du fxSDK :

```sh
fxsdk send-cg RMNChem.g3a
```

Il monte le volume, copie, démonte. Pratique quand on itère.

---

## Utilisation

**Écran de saisie.** Les flèches haut/bas changent de champ, gauche/droite
décrémentent et incrémentent, les chiffres saisissent directement. La formule
brute et le degré d'insaturation se mettent à jour en direct, avec le nombre
d'isomères trouvés. Les raccourcis F1 à F4 chargent des formules d'exemple,
F5 remet à zéro, F6 ouvre l'aide-mémoire des zones.

**Écran résultat.** En haut la formule, le nom et le DI. À gauche la formule
développée, groupes caractéristiques en rouge, avec le cercle aromatique pour
les cycles benzéniques. À droite le tableau des signaux : déplacement,
intégration, multiplicité, attribution. En bas le spectre de 12 à 0 ppm avec la
courbe d'intégration en orange.

Les protons labiles sont marqués d'un astérisque, avec le rappel qu'ils
disparaissent à l'échange par D₂O.

F1 et F2 naviguent entre isomères quand il y en a plusieurs.

**Mode secours.** Pour une formule absente de la base : DI exact, liste des
fragments compatibles, nombre d'hydrogènes à placer, et l'échelle des six zones
de déplacement usuelles.

**Formules impossibles.** Si le DI ressort demi-entier, l'add-in le signale
plutôt que d'arrondir : ça veut dire ion ou radical, ou une erreur de saisie.
Idem si b dépasse 2a + 2 + d.

---

## Dépannage

**Écran noir au lancement, il faut retirer les piles.** Le nom interne de
l'add-in est en cause sur certaines versions d'OS. Change `INTERNAL` dans
`project.cfg` (`@RMNCHM2` par exemple) et recompile.

**L'add-in n'apparaît pas dans le menu.** Le `.g3a` doit être à la racine de la
mémoire de stockage. Vérifie aussi que le fichier n'a pas été renommé avec une
extension en double par le gestionnaire de fichiers.

**Erreur CMake sur `find_package(Gint 2.9 REQUIRED)`.** gint n'est pas installé
ou l'est dans une version antérieure. Relance `giteapc install -u
Lephenixnoir/gint`.

**Erreur sur `generate_g3a`.** Les chemins d'icônes sont relatifs à la racine
du projet. Vérifie que `assets-cg/` existe et contient bien deux PNG de 92×64.

---

## Variante MicroPython

`micropython/rmnchem.py` tourne sans rien installer : copie le fichier à la
racine de la calculatrice vue en clé USB, puis ouvre-le depuis le menu Python.

C'est une version réduite. Base de 16 molécules au lieu de 46, tracé plus lent
parce que `casioplot` n'expose que `set_pixel` et oblige à réimplémenter le
tracé de segments, et surtout saisie au prompt faute d'API clavier. La formule
développée n'est pas dessinée, seulement le spectre et le tableau.

Utile pour dépanner ou pour vérifier un calcul de DI. Pour le reste, la version
C fait tout mieux.

---

## Organisation des fichiers

```
RMNChem/
├── CMakeLists.txt          build fxSDK, cible fx-CG 50
├── project.cfg             métadonnées de l'add-in
├── src/
│   ├── main.c              interface, tracé, machine à états
│   └── chem_db.h           structures compactes et base des 46 molécules
├── assets-cg/
│   ├── icon-uns.png        icône menu, non sélectionnée (92×64)
│   └── icon-sel.png        icône menu, sélectionnée (92×64)
├── tools/
│   ├── check_db.c          validateur de la base, à compiler avec gcc
│   └── sim/                simulateur de rendu sur PC
└── micropython/
    └── rmnchem.py          variante sans compilation
```

---

## Vérification de la base

`tools/check_db.c` recompte tout : nombre de C, O et N contre les indices
déclarés, somme des hydrogènes portés par les atomes contre b, somme des
intégrations du spectre contre b, valences de chaque atome, connexité du
squelette par union-find, validité des index de liaison, parité du degré
d'insaturation, doublons de noms.

```sh
gcc -std=c11 -Wall -Wextra -Wno-unused-const-variable \
    -Isrc tools/check_db.c -o /tmp/check_db && /tmp/check_db
```

Sortie attendue : `OK (0 erreur)` sur les 46 molécules.

Lance-le après chaque ajout. C'est lui qui a rattrapé les erreurs pendant le
développement, notamment sur le groupe nitro où l'oxygène anionique a bien une
valence 1 dans la forme N⁺(=O)–O⁻.

Le simulateur de `tools/sim/` rejoue les primitives de dessin sur PC et sort un
PNG par écran, en signalant les débordements de texte. Il évite de flasher la
calculatrice à chaque essai de mise en page.

---

## Ajouter une molécule

Dans `chem_db.h`, trois tableaux et une entrée.

Les atomes en premier, coordonnées relatives en pixels, deux entiers signés sur
8 bits. Les macros `ZX(i)` et `ZY(i)` donnent les positions d'une chaîne en
zigzag, `HEX6` et ses variantes placent un cycle à six sommets. Le champ `nH`
porte les hydrogènes, `show` indique s'il faut écrire l'étiquette, `hl` colore
en rouge les groupes caractéristiques.

Les liaisons ensuite : deux index d'atomes, un ordre, un drapeau aromatique.

Les signaux enfin : déplacement en dixièmes de ppm, intégration, multiplicité
prise dans `M_S` à `M_M`, drapeau labile, et un index dans le tableau `ATTR[]`
des attributions partagées. Réutilise une chaîne existante quand c'est
possible, chaque nouvelle coûte sa longueur en Flash.

Puis l'entrée `MOL(...)` dans `MOLECULES[]`, qui calcule les tailles toute
seule par `sizeof`.

Si la formule existe déjà, ajoute-la au groupe d'isomères correspondant pour
qu'elle apparaisse dans la pagination.

Recompile le validateur avant de recompiler l'add-in.
