#ifndef SRC_UTILS_HPP
#define SRC_UTILS_HPP
#define HAVE_EIGEN   


namespace McAmon {

// Container for command-line options
struct Option {

    std::string target;
    std::string ligand;

};


// Option parser
Option get_options (int argc, char **argv) {

    int c;

    Option opts;

    while (1) {
        static struct option long_options[] = {
            {"target",        required_argument, 0, 'a'},
            {"ligand",        required_argument, 0, 'b'},
            //{0, 0, 0, 0}
        };

        int option_index = 0;

        c = getopt_long (argc, argv, "a:b:", long_options, &option_index);

        if (c == -1)
            break;

        switch (c) {

            case 0:
                if (long_options[option_index].flag != 0) break;
            break;

            case 'a':
              opts.target = optarg;
              break;

            case 'b':
              opts.ligand = optarg;
              break;

            }
    }


  /* Print any remaining command line arguments (not options). */
  if (optind < argc)
    {
      printf ("non-option ARGV-elements: ");
      while (optind < argc)
        printf ("%s ", argv[optind++]);
      putchar ('\n');
    }


    if (opts.ligand.compare("") == 0) {
        printf("ERROR: No ligand xyz-file specified.\n");

        exit(0);
    }

    if (opts.target.compare("") == 0) {
        printf("ERROR: No target xyz-file specified.\n");

        exit(0);
    }

    return opts;


} // Namespace McAmon 


OpenBabel::OBMol readfile(std::string filename) {

    OpenBabel::OBMol mol;

    OpenBabel::OBConversion conv;
    conv.SetInAndOutFormats("xyz", "xyz");
    std::ifstream ifs;

    ifs.open(filename.c_str());
    conv.Read(&mol, &ifs);
    ifs.close();

    return mol;
}

OpenBabel::vector3 crossProduct(OpenBabel::vector3 v1, OpenBabel::vector3 v2) {

  OpenBabel::vector3 v;
  double X, Y, Z;

  X = v1.GetY() * v2.GetZ() - v1.GetZ() * v2.GetY();
  Y = v1.GetZ() * v2.GetX() - v1.GetX() * v2.GetZ();
  Z = v1.GetX() * v2.GetY() - v1.GetY() * v2.GetX();

  v.Set(X, Y, Z);

  return v;

}

double dotProduct(OpenBabel::vector3 v1, OpenBabel::vector3 v2) {

  double dot;

  dot = v1.GetX()*v2.GetX() + v1.GetY()*v2.GetY() + v1.GetX()*v2.GetZ();

  return dot;

}

double norm(OpenBabel::vector3 v) {
  double norm;

  norm = sqrt(v.GetX()*v.GetX() + v.GetY()*v.GetY() + v.GetZ()*v.GetZ()); 
  
  return norm;

}




OpenBabel::vector3 get_com(OpenBabel::OBMol &mol) {

  OpenBabel::vector3 com;
  com.Set(0.0, 0.0, 0.0);
  OpenBabel::OBAtom *atom;

  OpenBabel::vector3 temp;

  for (unsigned int i = 1; i < mol.NumAtoms() + 1; i++) {
    atom = mol.GetAtom(i);
    temp = atom->GetVector();
    com += temp;
  }

  com /= mol.NumAtoms();

  return com;

}

void move_molecule(OpenBabel::OBMol &mol, OpenBabel::vector3 move, int start_id = 1, int end_id = -1) {

  if (end_id == -1) end_id = mol.NumAtoms() + 1;

  OpenBabel::vector3 temp;
  OpenBabel::OBAtom *atom;

  for (int i = start_id; i < end_id; i++) {

    atom = mol.GetAtom(i);
    temp = atom->GetVector();
    temp += move;
    atom->SetVector(temp);
  }

}

OpenBabel::vector3 rotate(const OpenBabel::vector3 &V, const OpenBabel::vector3 &J, const double T) {

  double x = V.x();
  double y = V.y();
  double z = V.z();

  double u = J.x();
  double v = J.y();
  double w = J.z();

  double norm = std::sqrt(u*u + v*v + w*w);
  double inv_norm_sqrt = 1.0 / (norm*norm);
  double sint = std::sin(T);
  double cost = std::cos(T);

  double a = (u * (u*x + v*y + w*z) + (x * (v*v + w*w) - u * (v*y + w*z)) * cost + norm * (-w*y + v*z) * sint) * inv_norm_sqrt;
  double b = (v * (u*x + v*y + w*z) + (y * (u*u + w*w) - v * (u*x + w*z)) * cost + norm * ( w*x - u*z) * sint) * inv_norm_sqrt;
  double c = (w * (u*x + v*y + w*z) + (z * (u*u + v*v) - w * (u*x + v*y)) * cost + norm * (-v*x + u*y) * sint) * inv_norm_sqrt;

  OpenBabel::vector3 rotated;
  rotated.Set(a, b, c);

  return rotated;

}

void rotate_molecule(OpenBabel::OBMol &mol, OpenBabel::vector3 direction, double theta, int start_id = 1, int end_id = -1) {

  if (end_id == -1) end_id = mol.NumAtoms() + 1;

  OpenBabel::OBAtom *atom;
  OpenBabel::vector3 temp;

  for (int i = start_id; i < end_id; i++) {
    atom = mol.GetAtom(i);
    temp = atom->GetVector();
    temp = rotate(temp, direction, theta);
    
    atom->SetVector(temp);
  }

}


void rotate_to_yz_plane(OpenBabel::OBMol &mol) {

  OpenBabel::vector3 v1, v2, v, ex;                                
                                                                   
  OpenBabel::OBAtom *atom1, *atom2, *atom3;                        
                                                                   
  atom1 = mol.GetAtom(1);                                          
  atom2 = mol.GetAtom(2);                                          
  atom3 = mol.GetAtom(3);                                          
                                                                   
  v1 = atom1->GetVector() - atom2->GetVector();                    
  v2 = atom3->GetVector() - atom2->GetVector();                    
                                                                   
  //v = -1*McAmon::crossProduct(v1, v2);                             
  v = -1*crossProduct(v1, v2);                             
  double norm_v, norm_ex, theta, dot;                              
                                                                   
  ex.Set(1.0, 0.0, 0.0);                                           
  //dot = McAmon::dotProduct(v, ex);                                 
  dot = dotProduct(v, ex);                                 
                                                                   
  //norm_v  = McAmon::norm(v);                                       
  norm_v  = norm(v);                                       
  //norm_ex = McAmon::norm(ex);                                      
  norm_ex = norm(ex);                                      
  theta = acos(dot / (norm_v*norm_ex));                            
                                                                   
  //McAmon::rotate_molecule(mol, McAmon::crossProduct(v, ex), theta);
  McAmon::rotate_molecule(mol, crossProduct(v, ex), theta);

}

void write_xyz(OpenBabel::OBMol mol, std::string filename) {

  OpenBabel::OBConversion conv;
  conv.SetInAndOutFormats("xyz", "xyz");

  std::ofstream ofs(filename);

  conv.Write(&mol, &ofs);
  ofs.close();
}

std::string get_fout(std::string name) {

  std::string delimiter = "/";
  size_t pos = 0;
  std::string token;

  while ((pos = name.find(delimiter)) != std::string::npos) {
    token = name.substr(0, pos);
    //std::cout << token << std::endl;
    name.erase(0, pos + delimiter.length());

  }
  //std::cout << name << std::endl;

  return name;
}

void minimize_molecule(OpenBabel::OBMol &mol) {
  OpenBabel::OBForceField* pFF = OpenBabel::OBForceField::FindForceField("MMFF94");

  pFF->Setup(mol);
  double e = pFF->Energy();

  const int steps = 200;
  const double crit = 5.0e-4;

  pFF->ConjugateGradientsInitialize(steps, crit);

  bool done = true;

  while (done) {
    done = pFF->ConjugateGradientsTakeNSteps(1);
    if (pFF->DetectExplosion()) {
      std::cerr << "explosion has occured!" << std::endl;
      exit(1);
    } else {
      pFF->GetCoordinates(mol);
    }
  }

  e = pFF->Energy();
  mol.SetEnergy(e);

}

void center_molecule(OpenBabel::OBMol &mol) {
  OpenBabel::vector3 v;
  OpenBabel::OBAtom *atom;

  atom = mol.GetAtom(2);

  v = -1*atom->GetVector();

  move_molecule(mol, v);

}


}

#endif

