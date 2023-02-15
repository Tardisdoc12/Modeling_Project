#include "nano.h"
using namespace std;

string dateTime(){
  int rc;
  time_t temp;
  struct tm *timeptr;
  temp = time(NULL);
  timeptr = localtime(&temp);
  char s[100];
  strftime(s,sizeof(s),"%b-%d-%r", timeptr);
  string str=s;
  return str;
}

string creates_directory(string name){
  string date=dateTime();
  date.pop_back();
  date.pop_back();
  date.pop_back();
  int status;
  string command="mkdir -p ./"+name+"-"+date;
  const char * com = command.c_str();
  status = system(com); // Creating a directory
  if (status == -1){
    cerr << "Error : Directory not created"<< endl;
  }
  return name+"-"+date;
}




namespace nano{
  namespace{
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
    float dcut=1/sqrt(2)+0.001;
    float energy_sum;
    float c2sum;
    int rejected=0;
    float temp,dmu,ddmu;
    int npas,npw,imax,jmax,kmax,npeq,idmumax,iconf,irand;
    vector<float> energy_atoms;
    float cbol=8.62e-5;
    bool verification_initialized=false;
    string impurity;
    string base;
    float V;
    float Tau;
  }


//functions of nano :  //nano:://---------------------------

bool open_file(string filename){
  ifstream fichier_entree(filename);
  if(!fichier_entree){
    ofstream fichier(filename);
    fichier.close();
    return false;
  }
  else{
    fichier_entree.close();
    return true;
  }
}

//construction of all parameters (Needed for starting calcul)
void initialize(string filename){
  vector<string> elements;
  if(!open_file(filename)){
    cerr<<"ERROR cannot open file "<<filename;
  }
  else{
    ifstream in(filename);
    string _r;
    getline(in,_r);
    string line;
    cout<<"Charging the paramters... \n Please wait..."<<endl;
    while(getline(in,line)){
      vector<string> inLine=split(line);
      elements.push_back(inLine[0]);
    }
    in.close();
  }
  npas=stoi(elements[0]);npw=stoi(elements[1]);temp=stof(elements[2]);
  dmu=stof(elements[3]);idmumax=stoi(elements[4]);ddmu=stof(elements[5]);
  imax=stoi(elements[6]);jmax=stoi(elements[7]);kmax=stoi(elements[8]);iconf=stoi(elements[9]);
  irand=stoi(elements[10]);npeq=20*npw;
  verification_initialized=true;
  cout<<"Parameters loaded.\n Loading of V and Tau..."<<endl;
  //properties of the system :
  //potential V
  if(!open_file("../properties/systems.dat")){
    cerr<<"IMPORTANT ERROR : File \"systems.dat\" not found!";
  }
  else{
    impurity=elements[12];base=elements[11];
    ifstream file("../properties/systems.dat");
    string line;
    string system=base+"-"+impurity;
    vector<string> datas;
    while(getline(file,line)){
      datas=split(line);
      if(datas[0]==system){
        cout<<"V ="<<datas[1]<<endl;
        V=stof(datas[1]);
      }
    }
    file.close();
  }
  //Tau
  if(!open_file("../properties/elements.dat")){
    cerr<<"IMPORTANT ERROR : File \"elements.dat\" not find! Please verify the existence of the file!";
  }
  else{
    ifstream elem("../properties/elements.dat");
    vector<string> element;
    string line;
    while(getline(elem,line)){
      element=split(line);
      if((element[0]==base)||(element[0]==impurity)){
        Tau=Tau+stof(element[2]);
      }
    }
    Tau=Tau/12.;
    cout<<"Tau vaut ="<<Tau<<endl;
    elem.close();
  }
  cout<<"End of loading V and Tau."<<endl;
}

vector<string> split(const string com){
  vector<string> result;
  const regex re("\\s+");
  transform(sregex_token_iterator(com.begin(), com.end(), re,-1), {}, back_inserter(result),[](string s) { return regex_replace(s, regex("\""), ""); });
  return result;
}

int randomf(int range){
  srand (time(NULL));
  return rand()%range;
}

string name_file(int mnomi){
  string filename;
  string mnomi_str=to_string(mnomi);
  filename="p"+mnomi_str+".dat";
  return filename;
}

double distance(float xij,float yij, float zij){
  float xij2=pow(xij,2);
  float yij2=pow(yij,2);
  float zij2=pow(zij,2);
  double dist=sqrt(xij2+yij2+zij2);
  return dist;
}

int mc_exchange(Maille& maille,int ipas){
  //our initial energy--------------------------------------------------------------------------------
  int n_impu=0;
  int Natom_tot=maille.getNumberOfAtoms();
  int Natom_base=maille.getNumberOfAtoms();
  //we are keeping data about energy
  vector<float> energy_atom0;
  copy(energy_atoms.begin(), energy_atoms.end(), energy_atom0.begin());

  //-------------------------------------------------------------------------------------------------
  for(int times=0;times<Natom_tot;times++){
    //random atom to pick
    int rand_atom=randomf(Natom_tot-1);
    string save_type=maille.getParticle(rand_atom).getKind();
    maille.changeParticle(rand_atom,impurity);
    Natom_base--;
    n_impu--;
    //calculate the new energy
    float ener_new=maille.energy();
    //calculate (thanks to boltzman term compared to a random number)
    float de = ener_new+Natom_base*dmu;
    if(de<=0){
      continue;
    }
    float boltzman = exp(-de/(cbol*temp));
    float p=randomf(RAND_MAX)/(float) RAND_MAX;
    if(p>boltzman){
      continue;
    }
    //we wait for an equilibrium of the system than start to count how many impurities were rejected
    if(ipas>npeq){
      rejected++;
    }
    //we return to the previous state
    copy(energy_atom0.begin(), energy_atom0.end(), energy_atoms.begin());
    maille.getParticle(rand_atom).change_type(save_type);
  }
  return Natom_base;
}

float Monte_Carlo(Maille& maille){
  string Directory_name="Simulation";
  string dir=creates_directory(Directory_name);
  maille.voisin();
  if(!verification_initialized){
    cerr<<"Error : parameters not initialized! Please use the function nano::initialize(string filename) to initialize them!";
  }
  for(int pas=0;pas<npas;pas++){
    int N_atom1=mc_exchange(maille,pas);
    if(pas>npeq){
      //Calculs about the energy of the exchange----
      float energy_total = maille.energy();
      energy_sum=energy_sum+energy_total;
      c2sum=c2sum+(1-N_atom1/maille.getNumberOfAtoms());
      //Writing into the file for some configurations
      string filename=name_file(pas);
      string path="./"+dir+"/"+filename;
      maille.write_parameters(path,pas);
    }
  }
  return energy_sum;
}
//end of functions in nano------------------------------------



//class Boxe{
  Boxe::Boxe(){
    x=0;y=0;z=0;
  }

  void Boxe::setDimension(string line){
    vector<string> dim=split(line);
    dim.erase(dim.begin());
    x=stof(dim[0]);y=stof(dim[1]);z=stof(dim[2]);
  }

  vector<float> Boxe::getDimension(){
    vector<float> dim={x,y,z};
    return dim;
  }

  Boxe::~Boxe(){

  }
//};


//class Nanoparticle{
  Nanoparticle::Nanoparticle(){
    x=0;y=0;z=0;
    kind="Null";
    apar=0;ecoh=0;
  }

  void Nanoparticle::set_param(string line){
    vector<string> true_sous_maille=split(line);
    true_sous_maille.erase(true_sous_maille.begin());
    kind=true_sous_maille[0];
    x=stof(true_sous_maille[1]);y=stof(true_sous_maille[2]);z=stof(true_sous_maille[3]);
    if(!open_file("../properties/elements.dat")){
      cerr<<"IMPORTANT ERROR : File \"elements.dat\" not find! Please verify the existence of the file!";
    }
    else{
      ifstream file("../properties/elements.dat");
      string line;
      vector<string> datas;
      while(getline(file,line)){
        datas=split(line);
        if(datas[0]==kind){
          apar=stof(datas[1]);
          ecoh=stof(datas[2]);
        }
      }
      file.close();
    }
  }

  vector<float> Nanoparticle::position(){
    vector<float> position={x,y,z};
    return position;
  }

  string Nanoparticle::getKind(){
    return kind;
  }

  void Nanoparticle::change_type(string new_type){
    kind=new_type;
  }

  Nanoparticle::~Nanoparticle(){

  }
//};


//class nano::Maille{
  Maille::Maille(string filename){
    if (open_file(filename)){
      ifstream file(filename);
      string line;
      getline(file, line);
      N_atoms=stof(line);
      nvois=vector<int>(N_atoms,0);
      ivois=vector<vector<int>>(N_atoms,vector<int>(0,0));
      energy_atoms=vector<float>(N_atoms,0);
      getline(file, line);
      boxe.setDimension(line);
      while(getline(file,line)){
        Nanoparticle particle;
        particle.set_param(line);
        maille.push_back(particle);
      }
      file.close();
    }
    else{
      cerr<<"failed to open the file "<<filename;
    }

  }

  float Maille::energy(){
    float energy_total=0.;
    for(int i=0;i<N_atoms;i++){
      if(maille[i].getKind()==impurity){
        for(int j=0;j<nvois[i];j++){
          int voisin_k=ivois[i][j];
          if(maille[voisin_k].getKind()==impurity){
            energy_atoms[i]=energy_atoms[i]+V;
          }
        }
        energy_atoms[i]=energy_atoms[i]+nvois[i]*(Tau-V);
      }
      energy_total=energy_total+energy_atoms[i];
    }
    return energy_total;
  }

  void Maille::voisin(){
    for(int i=0;i<N_atoms-1;i++){
      for(int j=i+1;j<N_atoms;j++){
        //preparation of the limit conditions :--------------------------------------------
        float xij=maille[j].position()[0]-maille[i].position()[0];
        float yij=maille[j].position()[1]-maille[i].position()[1];
        float zij=maille[j].position()[2]-maille[i].position()[2];
        float boxe_x=boxe.getDimension()[0];
        float boxe_y=boxe.getDimension()[1];
        float boxe_z=boxe.getDimension()[2];
        //Limit conditions :----------------------------------------------------------------
        if(abs(xij+boxe_x)<abs(xij)){
          xij=xij+boxe_x;
        }
        if(abs(xij-boxe_x)<abs(xij)){
          xij=xij-boxe_x;
        }
        if(abs(yij+boxe_y)<abs(yij)){
          yij=yij+boxe_y;
        }
        if(abs(yij-boxe_y)<abs(yij)){
          yij=yij-boxe_y;
        }
        if(abs(zij+boxe_z)<abs(zij)){
          zij=zij+boxe_z;
        }
        if(abs(zij-boxe_z)<abs(zij)){
          zij=zij-boxe_z;
        }
        //calcul and verification if i is neighbor of j:
        if(distance(xij,yij,zij)<dcut){
          nvois[i]=nvois[i]+1;
          nvois[j]=nvois[j]+1;
          ivois[i].push_back(j);
          ivois[j].push_back(i);
        }
      }
    }
  }

  void Maille::write_parameters(string filename,int pas){
    if(pas%npw==0){
      ofstream params(filename);
      params<<N_atoms<<endl;
      params<<boxe.getDimension()[0]<<" "<<boxe.getDimension()[1]<<" "<<boxe.getDimension()[2]<<endl;
      for(int indi=0;indi<maille.size();indi++){
        params<<maille[indi].getKind()<<" "<<maille[indi].position()[0]<<" "<<maille[indi].position()[1]<<" "<<maille[indi].position()[2]<<endl;
      }
      params.close();
    }
  }

  vector<vector<int>> Maille::getIVois(){
    return ivois;
  }

  vector<Nanoparticle> Maille::getMaille(){
    return maille;
  }

  void Maille::changeParticle(int site,string new_kind){
    maille[site].change_type(new_kind);
  }

  Nanoparticle Maille::getParticle(int position){
    return maille[position];
  }

  float Maille::getNumberOfAtoms(){
    return N_atoms;
  }

  vector<int> Maille::getNVois(){
    return nvois;
  }

  Boxe Maille::getBoxe(){
    return boxe;
  }

  Maille::~Maille(){

  }
//};
};
