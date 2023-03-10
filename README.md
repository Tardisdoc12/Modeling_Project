# Modeling_Project <br />
Il faut avoir un compte github pour effectuer ce tutoriel. <br />
<br />
## Tutoriel d'utilisation de github : <br />
### Tutoriel Premier Téléchargement <br />
Pour le premier téléchargement, il faut télécharger le projet sur l'ordinateur. <br />
Tout d'abord, placez-vous dans le fichier où vous souhaitez avoir le fichier. Puis ouvrez un terminal. <br />
Il faut écrire : ```git clone url_du_Projet``` <br />
l'url est copiable depuis le bouton ```<> Code```. <br />
Normalement, le dossier devrait apparaitre avec tous les fichiers et sous-dossier.<br />
<br />
### Mettre à jour le projet. <br />
Lorsque vous voulez mettre à jour le fichier sur votre ordinateur, il n'y a pas besoin de retélécharger le projet en entier. <br />
Tout d'abord, il faut ouvrir un terminal et se placer dans le dossier du projet. <br />
Il suffit ensuite d'écrire : ```git pull``` et de se connecter via le terminal à son compte. <br />
<br />
### Effectuer des modifications sur un fichier. <br />
Evidemment, le but est que tout le monde puisse modifier le projet à sa guise! <br />
Lorsque les modifications ont été sauvegardé sur votre ordinateur, ouvrez un terminal dans le dossier du Projet. <br />
Dans le terminal, écrivez les commandes suivantes : <br />
```git add *``` pour ajouter toutes les modifications sur les fichiers du projets ou ```git add NameOfTheFile``` pour ajouter les modification sur 1 fichiers <br />
```git commit -m"la modification ou la correction apporter"```<br />
```git push```<br />
Il faudra ensuite se connecter sur votre compte pour finaliser l'opération. Voir la section ci dessous pour ça.<br />
<br />
### Se connecter pour les modifications sur le terminal.<br />
Si on essaie de se connecter lors d'un ```git push``` avec  le mot de passe utilisateur, cela nous renvoie une erreur.<br />
La solution est de créer un token pour cela, qui correspond à un mot de passe random mieux protéger!<br />
Pour le créer, il suffit d'aller sur notre profil dans ```settings```. <br />
Une fois dedans, aller dans ```developper settings```, puis dans ```Personnal access tokens``` et dans ```Tokens```.<br />
Cliquer sur ```Generate new token```, prendre la version ```classic```.<br />
Renter son mot de passe utilisateur.<br />
Il faut changer expiration par ```No expiration```. Cocher ensuite les cases souhaiter.<br />
Le token est ensuite générer. Il faut le sauvegarder à un endroit car pour tous les changements ce sera le même token.<br />
Il n'y a donc pas besoin d'en générer un à chaque modifications.<br />
<br />
## Organisation des fichiers (sujet à modification) : <br />
Pour chaque "projet", il y a un dossier correspondant. <br />
A l'intérieure de chaque dossier se trouve :<br/>
- Les codes en Fortran <br />
- les codes en C++ <br />
- Les tests si besoin <br />
- Un README.md qui expliquera le projet et problèmes rencontrés si besoin. (un historique)<br />
