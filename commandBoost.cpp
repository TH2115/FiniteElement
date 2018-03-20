 /* command line arguments using boost program options.
 */
#include <iostream>
#include <cstdlib>
#include "getNumbers.h"
using namespace std;

#include <boost/program_options.hpp>

// An alias to reduce typing
namespace po = boost::program_options;

void getVal (int argc, char* argv[]){

        po::options_description opts(
		"Test case value definitions.");
    opts.add_options()
        ("a", po::value<double>()->default_value(0.0),
				 "Geometric parameter describing the height of the plate, a*x*x + b*x + h1.")
        ("h1",  po::value<double>()->default_value(1.0),
				 "Height of the left side of the plate [m].")
        ("h2",  po::value<double>()->default_value(1.0),
				 "Height of the right side of the plate [m].")
        ("L",  po::value<double>()->default_value(2.0),
				 "Length of the plate [m].")
        ("t_p",  po::value<double>()->default_value(0.2),
				 "Thickness of plate [m].")
        ("Kxx",  po::value<double>()->default_value(250.0),
				 "Thermal conductivity [W/mk].")
        ("Kyy",  po::value<double>()->default_value(250.0),
				 "Thermal conductivity [W/mk].")
        ("Kxy",  po::value<double>()->default_value(0.0),
				 "Thermal conductivity [W/mk].")
        ("Nelx",  po::value<int>()->default_value(10),
				 "Number of elements in the x direction.")
        ("Nely",  po::value<int>()->default_value(5),
				 "Number of elements in the y direction.")
        ("T_BC",  po::value<char>()->default_value('L'),
				 "Side in which the temperature is applied.")
        ("q_BC",  po::value<char>()->default_value('R'),
				 "Side in which the heat flux is applied.")
        ("help",       "Print help message.");

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, opts), vm);
    po::notify(vm);

    if (vm.count("help")) {
        cout << "Test case value definition." << endl;
        cout << opts << endl;
        return 0;
    }

    const double a = vm["a"].as<double>();
    const double h1 = vm["h1"].as<double>();
    const double h2 = vm["h2"].as<double>();
    const double L = vm["L"].as<double>();
    const double t_p = vm["t_p"].as<double>();
    const double Kxx = vm["Kxx"].as<double>();
    const double Kyy = vm["Kyy"].as<double>();
    const double Kxy = vm["Kxy"].as<double>();
    const int Nelx = vm["Nelx"].as<int>();
    const int Nely = vm["Nely"].as<int>();
    const char T_BC = vm["T_BC"].as<char>();
    const char q_BC = vm["q_BC"].as<char>();

    /////// creating matrix D /////////////
        double* D = new double[2*2];
        D[0] = k_xx;
        D[1] = k_xy;
        D[2] = k_xy;
        D[3] = k_yy;
    //// creatimg D_test matrix to test if positive definite
        double* D_test = new double[2*2];
        bool pass;
        D_test[0] = k_xx;
        D_test[1] = k_xy;
        D_test[2] = k_xy;
        D_test[3] = k_yy;
        pass = checkPositiveDefinite(2,  D_test);

        if (pass == false){
            cout << "D matrix is not positive definite!" << endl;
        }


}



