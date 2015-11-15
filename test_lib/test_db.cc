/* -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

 * File Name : test_db.cc

 * Purpose :

 * Creation Date : 03-10-2012

 * Last Modified : Thu Oct  4 11:58:46 2012

 * Created By : shultz

 _._._._._._._._._._._._._._._._._._._._._.*/


#include "semble/semble_meta.h"
#include "semble/semble_key_val_db.h"
#include "hadron/hadron_npart_npt_corr.h"
#include "io/adat_xmlio.h"
#include "io/key_val_db.h"
#include <string>
#include <exception>
#include <vector>
#include <iostream>
#include "ensem/ensem.h"
#include "AllConfStoreDB.h"


using namespace SEMBLE;

typedef SembleExtendedKeyHadronNPartIrrep_t K;
typedef SembleMassOverlapData_t D;
typedef ADATIO::SerialDBKey<K> SK;
typedef ADATIO::SerialDBData<D> SD;


#define CFG 1
#define SIZE 50
#define BASE 0
#define STEP 2
#define NKEYS 5
// i think there are ~700 in the file

int main(void)
{
  std::string filename("someXML.xml");
  ADATXML::XMLReader xml(filename);
  ADATXML::Array<Hadron::KeyHadronNPartNPtCorr_t> NPoints;

  try
  {
    read(xml,"/Keys",NPoints);
  }
  catch(const std::string &e)
  {
    std::cerr << "ERR: " << __func__ << " " << e << std::endl;
    exit(1);
  }

  std::vector<Hadron::KeyHadronNPartIrrep_t> irrep_keys;

  for(int i = 0; i < NKEYS; i++)
  {
    irrep_keys.push_back(NPoints[i].npoint[1].irrep);   
    //   std::cout << irrep_keys.back() << std::endl;    
  }

  std::string particle_id("pion");

  std::vector<K> ExtKeys;
  std::vector<D> ExtData; 
  std::vector<Hadron::KeyHadronNPartIrrep_t>::const_iterator it;

  int i = BASE;
  for(it = irrep_keys.begin(); it != irrep_keys.end(); ++it , i+=STEP)
  {
    //   std::cout << i << " " << *it;
    ExtKeys.push_back(K(particle_id,*it));
    ENSEM::EnsemReal E,Z;
    E.resize(SIZE);
    Z.resize(SIZE);
    E = toScalar(double(i));
    Z = toScalar(double(STEP)) * E;
    ExtData.push_back(D(E,Z));

    //  std::cout << i << " E " << toScalar(ENSEM::mean(ExtData.back().E())) << std::endl;
    //  std::cout << i << " Z " << toScalar(ENSEM::mean(ExtData.back().Z())) << std::endl;
  }

  std::vector<int> cfg;
  for(int i = 0; i < CFG; ++i)
    cfg.push_back(i);

  FILEDB::AllConfStoreDB< SK , SD > db(cfg);
  std::string dbname("foobar.sdb");

  if(db.open(dbname, O_RDWR | O_TRUNC | O_CREAT, 0664) != 0)
  {
    std::cerr << __func__ << ": error opening dbase= " << dbname << std::endl;
    exit(1);
  }

  std::cout << "\n\n"
    << "*****************************************************************************\n"
    << "*****************************INSERTING DATA**********************************\n"
    << "*****************************************************************************\n"
    << std::endl;


  for(int i = 0; i < irrep_keys.size(); ++i)
  {
    SK key; 
    SD data;
    key.key() = ExtKeys[i];
    data.data() = ExtData[i];
    if(db.insert(key,std::vector<SD>(1,data)) == 0)
    {
      std::cout << __func__ <<": inserted data \n DATA: " << toScalar(ENSEM::mean(data.data().E())) 
        << " " << toScalar(ENSEM::mean(data.data().Z())) << "\n KEY: " << key.key();
    }
    else
    {
      std::cerr << __func__ << ": could not insert key \n" << key.key() << std::endl;
      exit(1);
    }
  }

  db.close();


  std::cout << "\n\n"
    << "*****************************************************************************\n"
    << "********************FINISHED WRITING DATABASE********************************\n"
    << "*****************************************************************************\n"
    << std::endl;


  if(db.open(dbname, O_RDONLY , 0400) != 0)
  {
    std::cerr << __func__ << ": error opening database " << dbname << std::endl;
    exit(1);
  }

  std::vector<SK> keys;
  db.keys(keys);

  if(keys.size() != irrep_keys.size())
  {
    std::cerr << __func__ << ": error recording proper number of keys " << std::endl;
    std::cerr << " put in " << irrep_keys.size() << " keys and got out " << keys.size() << " keys \n";
    exit(1);
  }

  std::cout << "\n\n"
    << "*****************************************************************************\n"
    << "**********************FINISHED READING KEYS**********************************\n"
    << "*****************************************************************************\n"
    << std::endl;

  std::cout << "\n\n Keys: \n" << std::endl;

  for(int i = 0; i < keys.size(); ++i)
    std::cout << keys[i].key() << std::endl;


  for(int i = 0; i < ExtKeys.size(); ++i)
  {
    std::vector<SD> dat;
    db.get(ExtKeys[i] , dat);
    if(!!!(ExtData[i] == dat.begin()->data()))
    {
      std::cout << i << std::endl;
      std::cout << "E_in.size() =  " << ExtData[i].E().size() << "\nE_out.size() = " << dat.begin()->data().E().size() << std::endl;
      std::cout << "Z_in.size() =  " << ExtData[i].Z().size() << "\nZ_out.size() = " << dat.begin()->data().Z().size() << std::endl;
      std::cout << "E_in = " << toScalar(ENSEM::mean(ExtData[i].E())) << " +/- " <<  toScalar(ENSEM::variance(ExtData[i].E())) << std::endl;
      std::cout << "E_out = " << toScalar(ENSEM::mean(dat.begin()->data().E())) << " +/- " <<  toScalar(ENSEM::variance(dat.begin()->data().E())) << std::endl;
      std::cout << "Z_in = " << toScalar(ENSEM::mean(ExtData[i].Z())) << " +/- " <<  toScalar(ENSEM::variance(ExtData[i].Z())) << std::endl;
      std::cout << "Z_out = " << toScalar(ENSEM::mean(dat.begin()->data().Z())) << " +/- " <<  toScalar(ENSEM::variance(dat.begin()->data().Z())) << std::endl;
      std::cerr << __func__ << ": couldnt get back the right data for key \n " << ExtKeys[i] << std::endl;
      exit(1);
    }
  }

  std::cout << "\n\n"
    << "*****************************************************************************\n"
    << "*************************TEST SUCCESSFUL*************************************\n"
    << "*****************************************************************************\n"
    << std::endl;



  return 0;

}




