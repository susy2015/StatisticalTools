#include "TheLimits.h"
#include "SusyScan.h"
#include "GeneratorMasses.h"
#include "simpleGenMasses.h"

#include <fstream>
#include <iostream>
#include <cmath>

void TheLimits::Fill(int argc, char** argv)
{
   for (int i = 1; i<argc; ++i)
   {
     add( new SusyScan(argv[i]) );
   }
}

void TheLimits::Fill(const std::string filelist, const std::string directory)
{
   std::ifstream masses_file;
   masses_file.open(filelist.c_str());
   std::string filename;
   while (1) {
//      GeneratorMasses * p = new GeneratorMasses;
      masses_file >> filename;
      if (!masses_file.good()) break;
      std::string file = directory.c_str() + filename;
      add( new SusyScan(file));
   }
   //std::cout << filelist << ": "<< _scan.size() <<std::endl;
   masses_file.close();
}

void TheLimits::match()
{
  for (std::vector<SusyScan*>::iterator it=_scan.begin();it!=_scan.end();++it){
    bool match = false;
    for (std::vector<GeneratorMasses*>::const_iterator gt=_masses.begin();gt!=_masses.end();++gt){
      if ((*it)->Mzero==(*gt)->Mzero && 
          (*it)->Mhalf==(*gt)->Mhalf && 
	  (*it)->Azero==(*gt)->Azero &&
	  (*it)->TanBeta==(*gt)->TanBeta &&
	  (*it)->Mu==(*gt)->Mu ) {
        (*it)->M1 = (*gt)->M1;
        (*it)->M2 = (*gt)->M2;
        (*it)->M3 = (*gt)->M3;
        (*it)->MGL = (*gt)->MGL;
        (*it)->MUL = (*gt)->MUL;
        (*it)->MB1 = (*gt)->MB1;
        (*it)->MSN = (*gt)->MSN;
        (*it)->MNTAU = (*gt)->MNTAU;
        (*it)->MZ1 = (*gt)->MZ1;
        (*it)->MW1 = (*gt)->MW1;
        (*it)->MHL = (*gt)->MHL;
        (*it)->MUR = (*gt)->MUR;
        (*it)->MB2 = (*gt)->MB2;
        (*it)->MEL = (*gt)->MEL;
        (*it)->MTAU1 = (*gt)->MTAU1;
        (*it)->MZ2 = (*gt)->MZ2;
        (*it)->MW2 = (*gt)->MW2;
        (*it)->MHH = (*gt)->MHH;
        (*it)->MDL = (*gt)->MDL;
        (*it)->MT1 = (*gt)->MT1;
        (*it)->MER = (*gt)->MER;
        (*it)->MTAU2 = (*gt)->MTAU2;
        (*it)->MZ3 = (*gt)->MZ3;
        (*it)->MHA = (*gt)->MHA;
        (*it)->MDR = (*gt)->MDR;
        (*it)->MT2 = (*gt)->MT2;
        (*it)->MZ4 = (*gt)->MZ4;
        (*it)->MHp = (*gt)->MHp;
	match = true;
      }	  
    }
    //if (!match) std::cout << "No match for M0="<<(*it)->Mzero
    //                      << ", M12="<<(*it)->Mhalf<<std::endl;
  }
}

void TheLimits::matchSimple()
{
  for (std::vector<SusyScan*>::iterator it=_scan.begin();it!=_scan.end();++it){
    bool match = false;
    for (std::vector<simpleGenMasses*>::const_iterator gt=_simpleGenMasses.begin();gt!=_simpleGenMasses.end();++gt){
      if ((*it)->Mzero==(*gt)->Mzero && 
          (*it)->Mhalf==(*gt)->Mhalf ){
        (*it)->Mgluino = (*gt)->Mgluino;
        (*it)->Msquark = (*gt)->Msquark;
        (*it)->Mlsp = (*gt)->Mlsp;
	match = true;
      }	  
    }
  }
}

void TheLimits::FillGeneratorMasses(std::string file)
{
   std::ifstream masses_file;
   masses_file.open(file.c_str());
   while (1) {
      GeneratorMasses * p = new GeneratorMasses;
      masses_file >> p->Mzero
                  >> p->Mhalf
                  >> p->TanBeta
                  >> p->Mu  
                  >> p->Azero
                  >> p->Mtop  
                  >> p->muQ	  
                  >> p->Q	  
                  >> p->M1	  
                  >> p->M2
                  >> p->M3	  
                  >> p->MGL	  
                  >> p->MUL	  
                  >> p->MB1	  
                  >> p->MSN	  
                  >> p->MNTAU	  
                  >> p->MZ1	  
                  >> p->MW1	  
                  >> p->MHL	  
                  >> p->MUR	  
                  >> p->MB2	  
                  >> p->MEL	  
                  >> p->MTAU1	  
                  >> p->MZ2	  
                  >> p->MW2	  
                  >> p->MHH	  
                  >> p->MDL	  
                  >> p->MT1	  
                  >> p->MER	  
                  >> p->MTAU2	  
                  >> p->MZ3	  
                  >> p->MHA	  
                  >> p->MDR	  
                  >> p->MT2	  
                  >> p->MZ4	  
                  >> p->MHp;

      if (!masses_file.good()) break;
      //std::cout << p->Mzero<<", "<<p->Mhalf<<", "<<p->TanBeta<<std::endl;
      if (fabs(p->Mu)!=1.) {
         std::cerr << "check lines near m0=" << p->Mzero << ", m1/2=" << p->Mhalf << std::endl;
         break;
      }	
      _masses.push_back( p );
   }
}

void TheLimits::FillSimpleGenMasses(std::string file)
{
   std::ifstream masses_file;
   masses_file.open(file.c_str());
   while (1) {
      simpleGenMasses * p = new simpleGenMasses;
      masses_file >> p->Mzero
                  >> p->Mhalf
                  >> p->Mgluino
                  >> p->Mlsp
                  >> p->Msquark;

      if (!masses_file.good()) break;
      _simpleGenMasses.push_back( p );
   }
}

void TheLimits::FillEmptyPointsByInterpolationInM0M12()
{
  //first find out where to expect points
  std::cout<< "start: TheLimits::FillEmptyPointsByInterpolationInM0M12()" <<std::endl;
  double gridy=9999, miny=9999, maxy=0, gridx=9999, minx=9999, maxx=0;
  for (std::vector<SusyScan*>::const_iterator it=_scan.begin(); it!=_scan.end(); ++it){
    if ((*it)->Mzero<minx) minx=(*it)->Mzero;
    if ((*it)->Mzero>maxx) maxx=(*it)->Mzero;
    if ((*it)->Mhalf<miny) miny=(*it)->Mhalf;
    if ((*it)->Mhalf>maxy) maxy=(*it)->Mhalf;
    for (std::vector<SusyScan*>::const_iterator zt=it+1; zt!=_scan.end(); ++zt){
      if ( fabs((*it)->Mzero - (*zt)->Mzero) < gridx && 
           (*it)->Mhalf==(*zt)->Mhalf) gridx = fabs((*it)->Mzero - (*zt)->Mzero);
      if ( fabs((*it)->Mhalf - (*zt)->Mhalf) < gridy && 
           (*it)->Mzero==(*zt)->Mzero) gridy = fabs((*it)->Mhalf - (*zt)->Mhalf);
    }
  } 
  //Now, interpolate
  for (std::vector<SusyScan*>::const_iterator it=_scan.begin(); it!=_scan.end(); ++it){
     //find next neighbor in x and y for it:
    std::vector<SusyScan*>::const_iterator nextx=_scan.end(), nexty=_scan.end();
    double dx=9999, dy=9999; 
    for (std::vector<SusyScan*>::const_iterator zt=_scan.begin(); zt!=_scan.end(); ++zt){
      if (it==zt) continue;
      if ( fabs((*it)->Mzero - (*zt)->Mzero) < dx && 
           (*it)->Mzero < (*zt)->Mzero && (*it)->Mhalf==(*zt)->Mhalf) {
        dx = fabs((*it)->Mzero - (*zt)->Mzero);
	nextx = zt;
      }	
      if ( fabs((*it)->Mhalf - (*zt)->Mhalf) < dy && 
           (*it)->Mhalf < (*zt)->Mhalf && (*it)->Mzero==(*zt)->Mzero) {
        dy = fabs((*it)->Mhalf - (*zt)->Mhalf);
	nexty = zt;
      }
      if (dx==gridx && dy==gridy) break;	
    }
    //interpolate in x:
    if (dx!=gridx && nextx!=_scan.end()){
        //std::cout << "m0 = " <<(*it)->Mzero  << ", m12="<< (*it)->Mhalf << std::endl;
       double dist = (*nextx)->Mzero - (*it)->Mzero;
       for (double x=gridx; x<dist; x+=gridx ){
         //std::cout << "adding m0 = " <<x << ", m12="<< (*it)->Mhalf << std::endl;
	 _scan.push_back( new SusyScan( ( (**it * (x/dist)) + (**nextx * (1.-x/dist)) )));
       }	 
    }	

  }
  std::cout<< "done: TheLimits::FillEmptyPointsByInterpolationInM0M12()" <<std::endl;
}








