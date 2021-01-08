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

void write_xyz(OpenBabel::OBMol mol, std::string filename) {

  OpenBabel::OBConversion conv;
  conv.SetInAndOutFormats("xyz", "xyz");

  std::ofstream ofs(filename);

  conv.Write(&mol, &ofs);
  ofs.close();
}

}

#endif

