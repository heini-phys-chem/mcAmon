#define HAVE_EIGEN


#include <stdio.h>                 
                                   
#include <iostream>                
#include <sstream>
#include <vector>                  
#include <math.h>                  
#include <string.h>                
#include <fstream>                 
#include <sstream>                 
#include <random>                  
#include <getopt.h>                
                                   
#include <openbabel/babelconfig.h> 
#include <openbabel/base.h>        
#include <openbabel/mol.h>         
#include <openbabel/obconversion.h>
#include <openbabel/forcefield.h>  
#include <openbabel/obutil.h>      

#include "utils.hpp"


int main(int argc, char *argv[]) {

  // remove xyz files
  std::remove("*.xyz");

  // get cmd line options
  McAmon::Option opts = McAmon::get_options(argc, argv);

  // Read in Molecules
  OpenBabel::OBMol target = McAmon::readfile(opts.target);
  OpenBabel::OBMol ligand = McAmon::readfile(opts.ligand);

  // Define temp molecule to move
  OpenBabel::OBMol mol;
  
  // define vector for center of mass (cem)
  OpenBabel::vector3 com_ligand = McAmon::get_com(ligand);
  OpenBabel::vector3 com_target = McAmon::get_com(target);
  OpenBabel::vector3 vec;

  // move molecules to com and 2*com away
  McAmon::move_molecule(target, 2*com_target);
  McAmon::move_molecule(ligand, -com_ligand);

  // define rotation
  double theta = 90.0;
  OpenBabel::vector3 rot;
  rot.Set(0.0, 1.0, 0.0);

  // rotate molecule 90 degrees
  McAmon::rotate_molecule(ligand, rot, theta);
//  rot.Set(0.0, 1.0, 0.0);
//  McAmon::rotate_molecule(ligand, rot, theta);
//  rot.Set(1.0, 0.0, 0.0);
//  McAmon::rotate_molecule(ligand, rot, theta);

  // concatenate target and ligand
  OpenBabel::OBMol target_ligand = target;
  target_ligand += ligand;

  // get start and end ID of ligand in concatenated molecule
  int start_id   = target_ligand.NumAtoms() - ligand.NumAtoms() + 1;
  int end_id     = target_ligand.NumAtoms() + 1; 


  com_target = McAmon::get_com(target);

  // vector between both com's (target and ligand)
  vec = com_ligand - com_target;
  //McAmon::rotate_molecule(ligand, rot, theta, start_id, end_id);

  // write initial target-ligand complex
  McAmon::write_xyz(target_ligand, "test_mid.xyz");

  // array containing the scaling factors
  double arr[3] = { 1.5, 1.0, 0.5 };
  std::string f_out = "test_";

  // loop over scaling array, mvoe molecules, and write them to an xyz file
  for (int i = 0; i < 3; i++) {
    mol = target_ligand;
    McAmon::move_molecule(mol, arr[i]*vec, start_id, end_id);
    std::stringstream ss;
    ss << f_out << i << ".xyz";
    McAmon::write_xyz(mol, ss.str());
  }
}
