#include<iostream>
#include<cmath>
#include<fstream>
#include<string>
#include<vector>
#include<sstream>
#include<typeinfo>
#include<regex>
#include<algorithm>
#include<iterator>
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

const double dcut=1/sqrt(2)+0.01;
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


vector<string> split2(const string com){
  vector<string> result{};
  const std::regex re("\\s+");
  std::transform(std::sregex_token_iterator(com.begin(), com.end(), re,-1), {}, std::back_inserter(result),[](std::string s) { return std::regex_replace(s, std::regex("\""), ""); });
  result.erase(result.begin());
  return result;
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
    filename="p"+mnomi_str+".xyz";
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

//creation de la structure (atomes de toutes les mailles)---------------------------------------------------------
vector<vector<string>> maille_crea(string filename,vector<float>& boxe){
  vector<vector<string>>  maille; //dans chaque emplacement du vecteur de vecteur, il y aura : type atome,x,y,z
  if (open_file(filename)){
    ifstream file(filename);
    string line;
    getline(file, line);
    getline(file, line);
    vector<string> boxe_file=split2(line);
    for(int h=0;h<boxe_file.size();h++){
      boxe.push_back(stof(boxe_file[h]));
    }
    while(getline(file,line)){
      vector<string> true_sous_maille=split2(line);
      maille.push_back(true_sous_maille);
    }
    file.close();
    return maille;
  }
  else{
    return maille;
  }

}
//-----------------------------------------------------------------------------------------------------------------

//Function to calculate each neighbor------------------------------------------------------------------------------
double distance(float xij,float yij, float zij){
  float xij2=pow(xij,2);
  float yij2=pow(yij,2);
  float zij2=pow(zij,2);
  double dist=sqrt(xij2+yij2+zij2);
  return dist;
}

void voisin(vector<vector<string>> maille,vector<int> &nvois,vector<vector<int>> &ivois,vector<float>& boxe){
  for(int i=0;i<maille.size()-1;i++){
    for(int j=i+1;j<maille.size();j++){
      //preparation of the limit conditions :--------------------------------------------
      float xij=stof(maille[j][1])-stof(maille[i][1]);
      float yij=stof(maille[j][1])-stof(maille[i][1]);
      float zij=stof(maille[j][1])-stof(maille[i][1]);
      float boxe_x=boxe[0];
      float boxe_y=boxe[1];
      float boxe_z=boxe[2];
      //Limit conditions :----------------------------------------------------------------
      if(abs(xij+boxe_x)<xij){
        xij=xij+boxe_x;
      }
      if(abs(xij-boxe_x)<xij){
        xij=xij-boxe_x;
      }
      if(abs(yij+boxe_y)<yij){
        yij=yij+boxe_y;
      }
      if(abs(yij-boxe_y)<yij){
        yij=yij-boxe_y;
      }
      if(abs(zij+boxe_z)<zij){
        zij=zij+boxe_z;
      }
      if(abs(zij-boxe_z)<zij){
        zij=zij-boxe_z;
      }
      //calucl and verification if i is neighbor of j:
      if(distance(xij,yij,zij)<dcut){
        nvois[i]=nvois[i]+1;
        nvois[j]=nvois[j]+1;
        ivois[i].push_back(j);
        ivois[j].push_back(i);
      }
    }
  }
}
//-------------------------------------------------------------------------------------------------------------------

//Main---------------------------------------------------------------------------------------------------------------
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
  vector<vector<string>> maille;
  vector<float> boxe;
  maille=maille_crea("./init.dat",boxe);
  vector<int> nvois=vector<int>(maille.size(),0);
  vector<vector<int>> ivois=vector<vector<int>>(maille.size(),vector<int>(0,0));
  voisin(maille,nvois,ivois,boxe);
  for(int i=0;i<ivois[0].size();i++){
    cout<<ivois[0][i]<<endl;
  }
  return 0;
}
