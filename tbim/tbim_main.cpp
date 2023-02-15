#include "nano.h"
using namespace std;
using namespace nano;

//Main---------------------------------------------------------------------------------------------------------------
int main(){
  Maille maille("./init.dat"); //on créé la maille qui est notre system
  initialize("./in.dat"); //on initialise les paramètres
  Monte_Carlo(maille);
  return 0;
}
