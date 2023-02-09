#include<iostream>
#include<cmath>
#include<fstream>
#include<string>
#include<vector>
#include<sstream>
#include<typeinfo>
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
       *   cut distance                     : dcut          *
       *     dmu= mu_B - mu_A if we change A in B           *
       *                                                    *
       ******************************************************
*/

const double dcut=sqrt(2)+0.01;
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



vector<string> treat_line(string line, vector<string> &sous_maille){
  sous_maille=split(line,' ');
  cout<<"j"<<sous_maille[0]<<endl<<"a"<<sous_maille[1]<<endl<<"b"<<sous_maille[2]<<"c"<<endl;
  vector<string> true_maille;
  for(int i=0;i<sous_maille.size();i++){
    if(sous_maille[i]!=" "){
      true_maille.push_back(sous_maille[i]);
    }
  }
  return true_maille;
}
//----------------------------------------------------------------------------------------------------------------

//creation de la structure (atomes de toutes les mailles)
vector<vector<string>> maille_crea(string filename){
  vector<vector<string>>  maille; //dans chaque emplacement du vecteur de vecteur, il y aura : type atome,x,y,z
  if (open_file(filename)){
    ifstream file(filename);
    string line;
    getline(file, line);
    getline(file, line);
    vector<string> sous_maille;
    while(getline(file,line)){
      vector<string> true_sous_maille=treat_line(line,sous_maille);
      maille.push_back(true_sous_maille);
    }
    file.close();
    return maille;
  }
  else{
    return maille;
  }

}

/*double distance(vector<vector<string>> atome,int i,int j){
  double dist=sqrt(pow(stod(atome[i][1])-stod(atome[j][1]),2)+pow(stod(atome[i][2])-stod(atome[j][2]),2)+pow(stod(atome[i][3])-stod(atome[j][3]),2));
  return dist;
  }

void voisin(vector<vector<string>> maille){
    vector<int> nvois=vector<int>(maille.size(),0);
    vector<vector<int>> ivois=vector<vector<int>>(maille.size(),vector<int>(0,0));
    for(int i=0;i<maille.size();i++){
      for(int j=i+1;j<maille.size();j++){
        if(distance(maille,i,j)<dcut){
          nvois[i]=nvois[i]+1;
          nvois[j]=nvois[j]+1;
          ivois[i].push_back(j);
          ivois[j].push_back(i);
        }
      }
    }
}*/



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
  vector<vector<string>> a;
  a=maille_crea("./init.dat");
  for(int i=0;i<a[3].size();i++){
    cout<<a[3][i]<<endl;
  }
  cout<<typeid(a[3][0]).name()<<endl;


  return 0;
}
