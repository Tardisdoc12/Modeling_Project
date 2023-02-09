#include<iostream>
#include<cmath>
#include<fstream>
#include<vector>
#include<sstream>
using namespace std;
/*
       ******************************************************
       *                                                    *
       *                  PARAMETERS                        *
       *                                                    *
       *   maximal number of atoms          : nmax          *
       *   max. number of neighbours        : nvmax         *
       *   time steps number                : npas          *
       *   periodicity for writing          : npw           *
       *   temperature  (Kelvin)            : temp          *
       *   composition                      : cmet          *
       *   atoms number of metal 1          : nat1          *
       *   chemical potential difference    : dmu           *
       *     dmu= mu_B - mu_A if we change A in B           *
       *                                                    *
       ******************************************************
*/

//fonction pour ouvrir des fichiers-----------------------------------------------------------------------------
bool open_file(string filename){
  ifstream fichier_entree(filename);
  if(!fichier_entree){
    cout<<"ERROR file of the entry is not ready"<<endl;
    return false;
  }
  else{
    fichier_entree.close();
    return true;
  }
}
//---------------------------------------------------------------------------------------------------------------

//fonction pour split des chaines de caractères------------------------------------------------------------------
void split(string &chaine, char delimiteur, vector<string> &elements){
   stringstream ss(chaine);
   string sousChaine;
   while (getline(ss, sousChaine, delimiteur))
   {
     elements.push_back(sousChaine);
   }
}

vector<string> split( string &chaine, char delimiteur) {
 vector<string> elements;
 split(chaine, delimiteur, elements);
 return elements;
}
//---------------------------------------------------------------------------------------------------------------


//fonction pour construire les différents paramètres nécessaires:------------------------------------------------
void construction_parameters(vector<float>& elements,string filename){
  if(!open_file(filename)){
    cout<<"ERROR";
  }
  else{

    ifstream in(filename);
    string _r;
    getline(in,_r);
    string line;
    while(getline(in,line)){
      vector<string> inLine=split(line,' ');
      elements.push_back(stof(inLine[0]));
    }
    in.close();
  }
}
//---------------------------------------------------------------------------------------------------------------

//Fonction random------------------------------------------------------------------------------------------------
int randomf(){
  srand (time(NULL));
  return rand()%10000+1;
}
//---------------------------------------------------------------------------------------------------------------

//traduction de la subroutine subnomi----------------------------------------------------------------------------
string subnomi(int mnomi){
    string filename;
    int nbr_chiffre;
    string mnomi_str=to_string(mnomi);
    string filename_inter="";
    for(size_t i=0;i<mnomi_str.length();i++){
      filename_inter=filename_inter+"0";
    }
    filename="p"+filename_inter+".xyz";
    return filename;
}
//----------------------------------------------------------------------------------------------------------------

//Calcul de l'énergie---------------------------------------------------------------------------------------------
float energy(int natot){
  float energy_total;
  for(int i=0;i<natot;i++){
    float energy_i;
  }
  return energy_total;
}
//----------------------------------------------------------------------------------------------------------------
int main(){

  //initialisation avec les fichiers-----------------------------------------------------------------------------
  string filename="./in.dat";
  vector<float> elements;
  construction_parameters(elements,filename);
  float npas,npw,temp,dmu,idmumax,ddmu,imax,jmax,kmax,iconf,irand;
  npas=elements[0];npw=elements[1];temp=elements[2];dmu=elements[3];idmumax=elements[4];ddmu=elements[5];
  imax=elements[6];jmax=elements[7];kmax=elements[8];iconf=elements[9];irand=elements[10];float npeq=20*npw;
  cout<<"starting number ="<<imax<<endl;
  //fin de l'initialisation--------------------------------------------------------------------------------------
  cout<<subnomi(250)<<endl;
  //"Ceci est un changement"


  return 0;
}
