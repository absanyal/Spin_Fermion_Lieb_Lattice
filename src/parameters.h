#ifndef PARAMETERS_H
#define PARAMETERS_H
#include <cassert>

using std::ifstream;
using std::string;

typedef std::complex<double> cd;

class parameters
{
    // private:
    //     /* data */
public:
    int Lx, Ly, numw, Mx, My, nx, ny, seed, TCA, c_Lx, c_Ly, num_basis, k_only;
    double minw, maxw, JH, JH_d, G, t, tx, ty,
        td, d_theta_x, d_theta_y, d_phi_x, d_phi_y, T, T_stop,
        window, sweeps, filling, measure_after, measure_every,
        lambda, lambda_R, B, fixed_ky, Delta_c, Delta_s;
    bool plot_up_down_separately;

    double matchstring(string, string);
    void load(string);
};

double parameters::matchstring(string file, string match)
{
    string test;
    string line;
    ifstream readFile(file);
    double amount;
    bool pass = false;
    while (std::getline(readFile, line))
    {
        std::istringstream iss(line);
        if (std::getline(iss, test, '=') && pass == false)
        {
            // ---------------------------------
            if (iss >> amount && test == match)
            {
                // cout << amount << endl;
                pass = true;
            }
            else
            {
                pass = false;
            }
            // ---------------------------------
            if (pass)
                break;
        }
    }
    if (pass == false)
    {
        // string errorout = match;
        // errorout += "= argument is missing in the input file!";
        // throw std::invalid_argument(errorout);
        cout << match << "= argument is missing in the input file!";
    }
    cout << match << " = " << amount << endl;
    return amount;
}

void parameters::load(string inputfile)
{
    seed = matchstring(inputfile, "seed");
    k_only = matchstring(inputfile, "k_only");
    B = matchstring(inputfile, "B");
    Delta_c = matchstring(inputfile, "Delta_c");
    Delta_s = matchstring(inputfile, "Delta_s");
    fixed_ky = matchstring(inputfile, "fixed_ky");
    Lx = matchstring(inputfile, "Lx");
    Ly = matchstring(inputfile, "Ly");
    TCA = matchstring(inputfile, "TCA");
    c_Lx = matchstring(inputfile, "Lx");
    c_Ly = matchstring(inputfile, "Ly");
    JH = matchstring(inputfile, "JH");
    JH_d = matchstring(inputfile, "JH_d");
    t = matchstring(inputfile, "t");
    lambda = matchstring(inputfile, "lambda");
    lambda_R = matchstring(inputfile, "lambda_R");
    // ty = matchstring(inputfile, "ty");
    // td = matchstring(inputfile, "td");
    d_theta_x = matchstring(inputfile, "d_theta_x");
    d_theta_y = matchstring(inputfile, "d_theta_y");
    d_phi_x = matchstring(inputfile, "d_phi_x");
    d_phi_y = matchstring(inputfile, "d_phi_y");
    numw = matchstring(inputfile, "numw");
    minw = matchstring(inputfile, "minw");
    maxw = matchstring(inputfile, "maxw");
    G = matchstring(inputfile, "G");
    plot_up_down_separately = matchstring(inputfile, "plot_up_down_separately");
    T = matchstring(inputfile, "temperature_start");
    T_stop = matchstring(inputfile, "temperature_stop");
    window = matchstring(inputfile, "window");
    sweeps = matchstring(inputfile, "sweeps");
    measure_after = matchstring(inputfile, "measure_after");
    measure_every = matchstring(inputfile, "measure_every");
    filling = matchstring(inputfile, "filling");
}

#endif