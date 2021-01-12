#define HAVE_EIGEN


#include <stdio.h>                 
                                   
#include <iostream>                
//#include <sstream>
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
  OpenBabel::OBConversion conv;
  conv.SetInAndOutFormats("xyz", "xyz");

  mol = target;

  McAmon::center_molecule(mol);
  McAmon::rotate_to_yz_plane(mol);

  OpenBabel::vector3 tmp;
  tmp.Set(0.0, -1.0, 0.0);
  McAmon::move_molecule(mol, tmp);
  //OpenBabel::vector3 ex;
  //ex.Set(1.0, 0.0, 0.0);
  //double angle180  = 3.14159;

  //McAmon::rotate_molecule(mol, ex, angle180);

  target = mol;
  std::ofstream ofs_target("test_target.xyz");
  conv.Write(&target, &ofs_target);
  ofs_target.close();

  mol = ligand;

  McAmon::center_molecule(mol);
  McAmon::rotate_to_yz_plane(mol);

  OpenBabel::vector3 ex;
  ex.Set(1.0, 0.0, 0.0);
  double angle  = 2.44346;

  McAmon::rotate_molecule(mol, ex, angle);

  ligand = mol;
  std::ofstream ofs_ligand("test_ligand.xyz");
  conv.Write(&ligand, &ofs_ligand);
  ofs_ligand.close();


  OpenBabel::OBMol target_ligand = target;
  target_ligand += ligand;
  // get start and end ID of ligand in concatenated molecule
  int start_id   = target_ligand.NumAtoms() - ligand.NumAtoms() + 1;
  int end_id     = target_ligand.NumAtoms() + 1; 

  std::string target_name = opts.target;
  std::string f_out = McAmon::get_fout(target_name);
  // move molecules to com and 2*com away
  //McAmon::move_molecule(target, ex);
  //McAmon::move_molecule(ligand, -com_ligand);

  // array containing the scaling factors
  double arr[3] = { 2.0, 3.0, 4.0 };
  OpenBabel::vector3 vec;
  vec.Set(1.0, 0.0, 0.0);

  // write initial target-ligand complex
  std::ofstream ofs_move(f_out);

  // loop over scaling array, mvoe molecules, and write them to an xyz file
  for (int i = 0; i < 3; i++) {
    mol = target_ligand;
    McAmon::move_molecule(mol, arr[i]*vec, start_id, end_id);
    conv.Write(&mol, &ofs_move);
  }

  ofs_move.close();


//
//  // define rotation
//  double theta = 1.5708; // 90 degrees in radian
//  OpenBabel::vector3 rot;
//  rot.Set(1.0, 1.0, 1.0);
//
//  // rotate molecule 90 degrees
//  McAmon::rotate_molecule(ligand, rot, theta);
//
//  // concatenate target and ligand
//  OpenBabel::OBMol target_ligand = target;
//  target_ligand += ligand;
//
//  // minimize molecule
//  McAmon::minimize_molecule(target_ligand);
//
//  OpenBabel::OBConversion conv;
//  conv.SetInAndOutFormats("xyz", "xyz");
//  std::cout << f_out << std::endl;
//
//
//
//  com_target = McAmon::get_com(target);
//
//  // vector between both com's (target and ligand)
//  vec = com_ligand - com_target;
//
//  // array containing the scaling factors
//  double arr[3] = { -0.5, 0.5, 1.0 };
//
//  // write initial target-ligand complex
//  std::ofstream ofs_move(f_out);
//  conv.Write(&target_ligand, &ofs_move);
//
//  // loop over scaling array, mvoe molecules, and write them to an xyz file
//  for (int i = 0; i < 3; i++) {
//    mol = target_ligand;
//    McAmon::move_molecule(mol, -1*arr[i]*vec, start_id, end_id);
//    conv.Write(&mol, &ofs_move);
//  }
//
//  ofs_move.close();
}
