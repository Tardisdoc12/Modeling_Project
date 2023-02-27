#include "nano.h"
using namespace std;

string dateTime(){ //function to find the date and to convert it to a string
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

string creates_directory(string name){ //create a directory with the given name
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

namespace nano{ //namespace to ease the use for people
  namespace{ //anonymous namespace to protect the different parameters
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
    vector<float> c2sumvect(5,0.);//vector of concentration of all 5 type of particle (surface,heart..)
    vector<string> type_nano{"6","7","8","9","12"};//vector containing the type of the atoms depending on their positions
    int rejected=0;//number of impurities rejected when
    float temp,dmu,ddmu;
    int npas,npw,imax,jmax,kmax,npeq,idmumax,iconf,irand;
    vector<float> energy_atoms;
    float cbol=8.62*pow(10,-5);
    //float cbol=1.2*pow(10,-23);
    bool verification_initialized=false;
    string base;
    string impurity;
    float V;float Tau;
  }


//functions of nano :  //nano:://---------------------------

bool open_file(string filename){ //test the opening with ifstream
  ifstream fichier_entree(filename);
  if(!fichier_entree){
    ofstream fichier(filename); //if the opening ifstream we open it with ofstream
    fichier.close();
    return false;
  }
  else{
    fichier_entree.close();
    return true;
  }
}

//construction of all parameters (Needed for starting calcul)
void initialize(string filename){ //initialize all the parameters, it is mandatory!!
  vector<string> elements;
  if(!open_file(filename)){
    cerr<<"ERROR cannot open file "<<filename; //if there is no file send error
  }
  else{
    ifstream in(filename);
    string _r;
    getline(in,_r);
    string line;
    cout<<"Charging the paramters... \n Please wait..."<<endl; //inform the user that the charging is starting
    while(getline(in,line)){
      vector<string> inLine=split(line); //reading the file and we will attribute to all parameters
      elements.push_back(inLine[0]);
    }
    in.close();
  }
  npas=stoi(elements[0]);npw=stoi(elements[1]);temp=stof(elements[2]);
  dmu=stof(elements[3]);idmumax=stoi(elements[4]);ddmu=stof(elements[5]);
  imax=stoi(elements[6]);jmax=stoi(elements[7]);kmax=stoi(elements[8]);iconf=stoi(elements[9]);
  irand=stoi(elements[10]);npeq=20*npw;
  verification_initialized=true; //to let all function runs after!
  cout<<"Parameters loaded.\n Loading of V and Tau..."<<endl;
  //properties of the system :
  impurity=elements[12];base=elements[11];
  /*system_1.setParameters(base,impurity);
  if(impurity2!="Null"){
    system_2.setParameters(base,impurity2);
  }*/
  string fileName="../properties/elements.dat";
  if(open_file(fileName)){
    ifstream Element(fileName);
    string Line;
    vector<string> element;
    while(getline(Element,Line)){
      element=split(Line);
      if(element[0]==base){
        Tau=Tau+stof(element[2]);
      }
      else if(element[0]==impurity){
        Tau=Tau-stof(element[2]);
      }
    }
    Tau=Tau/12.;
  }
  //For system :----------------------------------------
  string filename2="../properties/systems.dat";
  if(open_file(filename2)){
    ifstream file(filename2);
    //for system:
    string line;
    while(getline(file,line)){
      vector<string> sys=split(line);
      if((base+"-"+impurity==sys[0])||(impurity+"-"+base==sys[0])){
        V=stof(sys[1]);
      }
    }
  }
  cout<<"End of loading V and Tau."<<endl;
}
//end of the initialization

vector<string> split(const string com){ //split a string by space
  vector<string> result;
  const regex re("\\s+"); //the space character then split the string
  transform(sregex_token_iterator(com.begin(), com.end(), re,-1), {}, back_inserter(result),[](string s) { return regex_replace(s, regex("\""), ""); });
  return result;
}

int randomf(int range){ //return a random number
  return rand()%range;
}

string name_file(int mnomi){ //create a new file name
  string filename;
  string mnomi_str=to_string(mnomi);
  filename="p"+mnomi_str+".xyz";
  return filename;
}

double distance(float xij,float yij, float zij){
  float xij2=pow(xij,2);
  float yij2=pow(yij,2);
  float zij2=pow(zij,2);
  double dist=sqrt(xij2+yij2+zij2); //define the distance
  return dist;
}

float energy(Maille& maille){
  float energy_total=0;
  for(int i=0;i<maille.getNumberOfAtoms();i++){
    if(maille.getParticleKind(i)==impurity){
      for(int j=0;j<maille.getNVois()[i];j++){
        int voisin_k=maille.getIVois()[i][j];
        if(maille.getParticleKind(voisin_k)==impurity){
          energy_atoms[i]=energy_atoms[i]+V;
        }
      }
      energy_atoms[i]=energy_atoms[i]+maille.getNVois()[i]*(Tau-V);
    }
    energy_total=energy_total+energy_atoms[i];
  }
  clear_vect(energy_atoms);
  return energy_total;
}

//exchange particles and calcul the new energy
int mc_exchange(Maille& maille,int ipas){
  srand (time(NULL));
  //our initial energy--------------------------------------------------------------------------------
  int n_impu=0;
  int Var=0;
  int Natom_tot=maille.getNumberOfAtoms();
  int Natom_base=maille.getNumberOfAtoms();
  float ener_0=0;
  //we are going to keep data about energy in this vector:
  vector<float> energy_atom0=vector<float>(maille.getNumberOfAtoms(),0);
  //-------------------------------------------------------------------------------------------------
  for(int times=0;times<Natom_tot;times++){
    //saving before exchange:
    copy(energy_atoms.begin(), energy_atoms.end(), energy_atom0.begin());
    //random atom to pick
    int rand_atom=randomf(Natom_tot-1);
    string save_type=maille.getParticleKind(rand_atom);
    ener_0=energy(maille);
    if(save_type!=impurity){
      Var=-1;
      n_impu++;
      maille.changeParticle(rand_atom,impurity);
    }
    else{
      Var=1;
      n_impu--;
      maille.changeParticle(rand_atom,base);
    }
    //---------------------------------------
    //modification of the kind of atom random
    //cout<<save_type<<" devient "<<maille.getParticleKind(rand_atom)<<" pour l'indice "<<rand_atom<<endl;
    //calculate the new energy
    float ener_new=energy(maille);
    //calculate (thanks to boltzman term compared to a random number)
    float de = ener_new;
    //cout<<"l'energie nouvelle vaut ="<<de<<endl;
    //cout<<"l'ancienne energie vaut ="<<ener_0<<endl;
    de-=ener_0;
    de+=Var*dmu;
    //cout<<"DE="<<de<<endl;
    if(de<=0){
      continue;
    }
    float boltzman = exp(-de/(cbol*temp));
    float p=(randomf(RAND_MAX))/(float) RAND_MAX;
    //cout<<"Boltzman vaut ="<<boltzman<<endl;
    //cout<<"p="<<p<<endl;
    if(p<boltzman){
      continue;
    }
    //we wait for an equilibrium of the system than start to count how many impurities were rejected
    if(ipas>npeq){
      rejected++;
    }
    //we return to the previous state
    copy(energy_atom0.begin(), energy_atom0.end(), energy_atoms.begin());
    Var=0;
    if(save_type==base){
      n_impu--;
    }
    else{
      n_impu++;
    }
    maille.changeParticle(rand_atom,save_type);
  }
  return Natom_base;
}

 void concen_neighbor(Maille maille){
  vector<int> Nimpu_vect(5,0);
  vector<int> Nbase_vect(5,0);
  for(int i=0;i<405;i++){
    for(int j=0;j<type_nano.size();j++){
      if(maille.getParticleType(i)==type_nano[j]){
        cout<<maille.getParticleType(i)<<endl;
        if(maille.getParticleKind(i)==impurity){
          Nimpu_vect[j]+=Nimpu_vect[j];
        }
        else{
          Nbase_vect[j]+=Nbase_vect[j];
        }
      }
    }
  }
  cout<<"base"<<Nbase_vect[2]<<endl;
  for(int j=0;j<type_nano.size();j++){
    c2sumvect[j]=(float)Nimpu_vect[j]/(float)(Nimpu_vect[j]+Nbase_vect[j]);
  }
}

void clear_vect(vector<float> vect){
  for(int i=0;i<vect.size();i++){
    vect[i]=0.;
  }
}


float Monte_Carlo(Maille& maille){ //make the monte carlo algorithm
  cout<<"je fais le monte"<<endl;
  for(int pas=0;pas<npas;pas++){ //start the algorithm
    cout<<"on est au pas="<<pas<<endl;
    int N_atom1=mc_exchange(maille,pas); //exchange atoms
    if(pas>npeq){
      //Calculs about the energy of the exchange----
      float energy_total = energy(maille); //calcul the energy
      energy_sum=energy_sum+energy_total; //energy final
      //Writing into the file for some configurations
    }
    if(pas==npas-1){
      cout<<"fini"<<endl;
      concen_neighbor(maille);
      cout<<"calculation of the concentrations..."<<endl;
    }
  }
  return energy_sum;
}
//end of functions in nano------------------------------------

//function to write parameters and to do all loops of monte carlo :
void DoMonteCarlo(Maille& maille){
  energy_atoms=vector<float>(maille.getNumberOfAtoms(),0); //initialize vectors to ease the program
  string Directory_name="Simulation_TOh_405_";
  string dir=creates_directory(Directory_name);
  create_concen_files(dir);//creates the concentration files
  maille.voisin();//create the ivois and nvois of the maille
  cout<<"voisin_done"<<endl;
  if(!verification_initialized){ //verify that all parameters are initialize and throw an error if not!
    cerr<<"Error : parameters not initialized! Please use the function nano::initialize(string filename) to initialize them!";
  }
  Maille maille_work=maille;
  for(int pasMu=0;pasMu<idmumax;pasMu++){
    //Maille maille_work=maille;
    cout<<"On exectute le monte carlo pour la boucle dmu="<<dmu<<endl;
    float energy=Monte_Carlo(maille_work);
    cout<<"On a fini le monte carlo pour dmu="<<dmu<<endl;
    string filename=name_file(pasMu);
    string path="./"+dir+"/"+filename;
    cout<<"On écrit les paramètres!"<<endl;
    maille_work.write_parameters(path,pasMu);//write the configuration in a file which is in the directory!
    write_all(dir);//write in all the concentration files
    cout<<"writing concentrations..."<<endl;
    clear_vect(c2sumvect);//clear the concentration vector so there is not overwritting
    dmu=dmu+ddmu;
    cout<<endl;
  }

}

void writeConcen(string path,float concen){
  ofstream file(path,ios::in|ios::ate|ios::out);
  file<<dmu<<"\t"<<concen<<endl;
  file.close();
}

void write_all(string directory){
  for(int j=0;j<type_nano.size();j++){
    writeConcen("./"+directory+"/concentration"+type_nano[j]+".dat",c2sumvect[j]);
  }
}

void create_concen_files(string directory){
  for(int j=0;j<type_nano.size();j++){
    string path2="./"+directory+"/concentration"+type_nano[j]+".dat";
    if(open_file(path2)){
      cerr<<"ERROR : FILE already exist please be sure of the name fileERROR : FILE already exist please be sure of the name file";
    }
  }
}


//class Boxe{
  /*Boxe::Boxe(){ //new class (new kind of variables like int, double etc...)
    x=0;y=0;z=0;
  }

  void Boxe::setDimension(string line){ //set the dimension of the boxe
    vector<string> dim=split(line);
    dim.erase(dim.begin());
    x=stof(dim[0]);y=stof(dim[1]);z=stof(dim[2]);
  }

  vector<float> Boxe::getDimension(){ //return the dimension (not constant but must be in futur)
    vector<float> dim={x,y,z};
    return dim;
  }

  Boxe::~Boxe(){ //destructor (useless here but can be usefull if we make inherit a new class)!

  }*/
//};


//class Nanoparticle{ //class Nanoparticle
  Nanoparticle::Nanoparticle(){ //creator for a variables of this kind
    x=0;y=0;z=0;
    kind="Null";
    apar=0;ecoh=0;
  }

  void Nanoparticle::set_param(string line){ // will set all parameters of the nanoparticle (need for the maille) the string will give the kind, the position (x,y,z)
    vector<string> true_sous_maille=split(line); //split the string to have the kind,x,y,z separate in a vector
    true_sous_maille.erase(true_sous_maille.begin()); // erase the first term that will be a space
    kind=true_sous_maille[0]; //attribute the different variables
    type_vois=true_sous_maille[4];
    x=stof(true_sous_maille[1]);y=stof(true_sous_maille[2]);z=stof(true_sous_maille[3]);
    if(!open_file("../properties/elements.dat")){ //will test to open the file for the properties of the kind of the nanoparticles
      cerr<<"IMPORTANT ERROR : File \"elements.dat\" not find! Please verify the existence of the file!";
    }
    else{ //if the file exist, lets go! we will read it.
      ifstream file("../properties/elements.dat");
      string line;
      vector<string> datas;
      while(getline(file,line)){ //we will need to read each lines
        datas=split(line);
        if(datas[0]==kind){ //if we find the good line
          apar=stof(datas[1]); //attribute the parameters of the nanoparticles
          ecoh=stof(datas[2]);
        }
      }
      file.close();
    }
  }

  vector<float> Nanoparticle::position(){ //will return the position! (again => must be const) (I or Yasmina will do it)
    vector<float> position={x,y,z};
    return position;
  }

  const string Nanoparticle::getKind(){//will return the kind (const)
    return kind;
  }

  const string Nanoparticle::getType_vois(){//will return the type depending on the number of neighbours
    return type_vois;
  }

  void Nanoparticle::change_type(string new_type){ //will change the type of the variable
    kind=new_type; //change the kind of the particle
    //find again the new properties of the particle!
    if(!open_file("../properties/elements.dat")){ //will test to open the file for the properties of the kind of the nanoparticles
      cerr<<"IMPORTANT ERROR : File \"elements.dat\" not find! Please verify the existence of the file!";
    }
    else{ //if the file exist, lets go! we will read it.
      ifstream file("../properties/elements.dat");
      string line;
      vector<string> datas;
      while(getline(file,line)){ //we will need to read each lines
        datas=split(line);
        if(datas[0]==kind){ //if we find the good line
          apar=stof(datas[1]); //attribute the parameters of the nanoparticles
          ecoh=stof(datas[2]);
        }
      }
      file.close();
    }
  }

  Nanoparticle::~Nanoparticle(){ //destructor of the class (useful for heritage)

  }
//};


//class nano::Maille{
  Maille::Maille(string filename){ //create the maille with the file initial
    if (open_file(filename)){//try to open the file
      ifstream file(filename);
      string line;
      getline(file, line);
      N_atoms=stof(line); //attribute the number of the file
      nvois=vector<int>(N_atoms,0);
      ivois=vector<vector<int>>(N_atoms,vector<int>(0,0));
      getline(file, line);
      //boxe.setDimension(line); //initialize the boxe of the maille
      while(getline(file,line)){ //read the file line by line
        Nanoparticle particle;
        particle.set_param(line);
        maille.push_back(particle); //push the particle corresponding to the line
      }
      file.close();
    }
    else{ //if the file cannot be open return an error
      cerr<<"failed to open the file "<<filename<<" in order to create the system (maille)";
    }

  }

  void Maille::voisin(){ //will make nvois and ivois
    for(int i=0;i<N_atoms-1;i++){
      for(int j=i+1;j<N_atoms;j++){
        //preparation of the limit conditions :--------------------------------------------
        float xij=maille[j].position()[0]-maille[i].position()[0];
        float yij=maille[j].position()[1]-maille[i].position()[1];
        float zij=maille[j].position()[2]-maille[i].position()[2];
        /*float boxe_x=boxe.getDimension()[0];
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
        }*/
        //calcul and verification if i is neighbor of j:
        if(distance(xij,yij,zij)<dcut){
          //if it is the case put the information in the vector ivois and nvois
          nvois[i]=nvois[i]+1;
          nvois[j]=nvois[j]+1;
          ivois[i].push_back(j);
          ivois[j].push_back(i);
        }
      }
    }
  }

  void Maille::write_parameters(string filename,int pas){ //will wirte parameters in the file
    float maille_parameter=3.62;
    //cout<<"on a écrit dans un fichier la configuration!"<<endl;
    ofstream params(filename);
    params<<N_atoms<<endl;
    //params<<boxe.getDimension()[0]<<" "<<boxe.getDimension()[1]<<" "<<boxe.getDimension()[2]<<endl;
    for(int indi=0;indi<maille.size();indi++){
      params<<maille[indi].getKind()<<" "<<maille[indi].position()[0]*maille_parameter<<" "<<maille[indi].position()[1]*maille_parameter<<" "<<maille[indi].position()[2]*maille_parameter<<" "<<maille[indi].getType_vois()<<endl;
    }
    params.close();
  }

  //FUNCTIONS THAT MUST BE CONST
  const vector<vector<int>> Maille::getIVois(){ //will return ivois
    return ivois;
  }

  const vector<Nanoparticle> Maille::getMaille(){ //will return maille
    return maille;
  }

  const Nanoparticle Maille::getParticle(int position){ //will return particle of the site position
    return maille[position];
  }

  const float Maille::getNumberOfAtoms(){ //return the number of atoms in the maille
    return N_atoms;
  }

  const vector<int> Maille::getNVois(){ //return number of voisin
    return nvois;
  }

  /*const Boxe Maille::getBoxe(){ //return the boxe of the maille
    return boxe;
  }*/

  const string Maille::getParticleKind(int sitePosition){
    return maille[sitePosition].getKind();
  }

  const string Maille::getParticleType(int sitePosition){
    return maille[sitePosition].getType_vois();
  }
  //END of the const functions

  void Maille::changeParticle(int site,string new_kind){ //will change the particle of the site i to the new kind
    maille[site].change_type(new_kind);
  }

  Maille Maille::operator=(const Maille& source){
    // Construit une copie temporaire de la source :
    Maille Temp(source);
    // Échange le contenu de cette copie avec l'objet courant :
    swap(Temp,*this);
    // Renvoie l'objet courant (modifié) et détruit les données
    // de la variable temporaire (contenant les anciennes données) :
    return *this;
  }

  Maille::~Maille(){ //destructor of the maille
    }
  }
//};
