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
#include <sys/stat.h>
#include <sys/types.h>
#include <iomanip>
#include <ctime>
using namespace std;

#ifndef __nano_h__
#define __nano_h__

namespace nano{
  //test of opening the file-----------------
  bool open_file(string filename);

  //split a string by spaces-----------------
  vector<string> split(const string com);

  //construction of all parameters (Needed for starting calcul)-----------
  void initialize(string filename);

  //randomizer-------------------------------
  int randomf(int range);

  //create the name of the file of the output:
  string name_file(int mnomi);

  //Function to calculate distance-----------
  double distance(float xij,float yij, float zij);

  //Function to write some parameters to plot courb in the future :
  void write_courbe(string filename);

  //Classe of nanoparticle-------------------
  class Nanoparticle{
  private:
    float apar; //=-0.0135;
    float ecoh; //=0.55/12.;
    float x,y,z;
    string kind;
  public:
    //creator---------------------------------
    Nanoparticle();

    //setting the parameters in the maille (kind,x,y,z) :
    void set_param(string line);

    //return position in maille :--------------
    vector<float> position();

    //change the kind of the nanoparticule-----
    void change_type(string new_type);

    //return the kind of the particle:---------
    string getKind();

    //Destructor------------------------------
    ~Nanoparticle();
  };

  //Classe boxe (limit conditions)------------
  class Boxe{
    private :
    float x,y,z;
    public :

    //Creator of the boxe---------------------
    Boxe();

    //set the dimension-----------------------
    void setDimension(string line);

    //return the dimension of the boxe--------
    vector<float> getDimension();

    //Destructor------------------------------
    ~Boxe();

  };

  //Classe Maille-----------------------------
  class Maille{
  private:
    vector<Nanoparticle> maille;
    vector<int> nvois;
    vector<vector<int>> ivois;
    Boxe boxe;
    float N_atoms;
    vector<float> energy_atoms;
  public:
    //creation de la structure ----------------
    Maille(string filename);

    //calcul of the energy of the maille :-----
    float energy();

    //create vectors ivois and nvois :---------
    void voisin();

    //return of ivois
    vector<vector<int>> getIVois();

    //change the kind of particle i :----------
    void changeParticle(int site,string new_kind);

    //Function to write all parameters in a file :
    void write_parameters(string filename,int pas);

    //return the entire maille-----------------
    vector<Nanoparticle> getMaille();

    //return a nanoparticle of the maille------
    Nanoparticle getParticle(int position);

    //return of nvois
    vector<int> getNVois();

    //return the number of atom in the maille--
    float getNumberOfAtoms();

    //return the boxe of the maille :
    Boxe getBoxe();

    //Destructor-------------------------------
    ~Maille();
  };

  //change an atom of the maille in the impurity :
  int mc_exchange(Maille& maille,int ipas);

  //Monte Carlo of a maille (maybe to put in class maille)---------
  float Monte_Carlo(Maille& maille);

}
#endif
