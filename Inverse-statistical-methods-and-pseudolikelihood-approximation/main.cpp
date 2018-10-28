//
//  main.cpp
//  Inverse-statistical-methods-and-pseudolikelihood-approximation
//
//  Created by Samuel Bosch on 11/10/18.
//  Copyright © 2018 Samuel Bosch. All rights reserved.
//

#include <iostream>
#include <vector>
#include <random>
#include <array>
#include <fstream>
#include <iomanip>
#include <string>

using namespace std;


int main(int argc, const char * argv[]) {
    random_device rd;  // only used once to initialise (seed) engine
    mt19937 rng(rd()); // random-number engine used (Mersenne-Twister in this case)
    int min_spin = 0;  // min value spin can take
    int max_spin = 21; // max value spin can take
    int N = 25;        //Number of atoms in molecule
    double T = 1;      // This is the temperature of the experiment
    // double K_b = 1.38064852e-23; // Boltzman constant
    double K_b = 1;
    uniform_int_distribution<int> random_spin(min_spin,max_spin); // definition of random function
    uniform_int_distribution<int> random_atom(0,N-1); // definition of random function

    array<std::array<double, 525>, 525> J;
    ifstream myfile;
    myfile.open("/Users/samuelbosch/OneDrive/Faks/EPFL_M1/Computer_simulation_of_physical_systems/Project/Inverse-statistical-methods-and-pseudolikelihood-approximation/J.txt");
    if(!myfile){ //Testing if file is actually open.
        cout << "Error opening output file" << endl;
        return -1;
    }else{
        cout << "File 'J.txt' is open" << "\n\n";
    }
    for (int i = 0; i < 525; i++){
        for(int j = 0; j < 525; ++j){
            myfile >> J[i][j]; // Here we read the file number by number
        }
    }

    // Random initialization of spins
    cout << "Initial spin configuration:\n";
    vector<int> v(N);    // declares a vector of integers
    for(int i=0; i<N; i++){
        auto random_integer = random_spin(rng);
        cout << random_integer << ' ';
        v[i] = random_integer;
    }

    // Calculation of the energy using Pott's model (H = -J*sum(Kronecker_delta(i,j))
    double E = 0.0;
    for(int i=0; i<N; i++){
        for(int j=i+1; j<N; j++){
            if (abs(J[21*i+v[i]][21*j+v[j]]) > 0){
                E = E - J[21*i+v[i]][21*j+v[j]]; // Is this the right formula for calculating the energy?
            }
        }
    }
    cout << "\nEnergy(initial)" << " = " << E << "\n\n";

    //cout << "\nStarting random iterations...\n";
    int max_number_of_interations = 10000;
    int iteration_number = 0;
    for(; iteration_number<max_number_of_interations; iteration_number++){
        auto atom_number = random_atom(rng); //Pick random atom for changing the spin
        int old_spin = v[atom_number]; //Saving old spin in case we still want to use it
        double E_old = E;
        v[atom_number] = random_spin(rng); //Pick random new spin for selected atom
        if(v[atom_number] == old_spin){
            continue; // If the spin didn't change, we do nothing
        }
        // Calculation of the new energy (new, more efficient algorithm)
        for(int i=0; i<N; i++){
            if (atom_number==i){
                continue;
            }
            if (abs(J[21*i+v[i]][21*atom_number+v[atom_number]]) > 0){
                E = E - J[21*i+v[i]][21*atom_number+v[atom_number]]; // Here we add the new energy
            }
            if (abs(J[21*i+v[i]][21*atom_number+old_spin]) > 0){
                E = E + J[21*i+v[i]][21*atom_number+old_spin]; // Here we substract the old energy
            }
        }

        if(E > E_old){
            double uniform_random_0_1 = (double)rand()/(double)RAND_MAX;
            double prob = exp(-(E-E_old)/(K_b*T));
            if (prob < uniform_random_0_1){
                // With some small probability, we revert the change (Metropolis algorithm)
                E = E_old;
                v[atom_number] = old_spin;
            }
        }
    }

    cout << "Final spin configuration:\n";
    for(int i=0; i<N; i++){
        cout << v[i] << ' ';
    }
    cout << endl;
    cout << "Energy(" << iteration_number << ")" << " = " << E << "\n";
    cout << '\n';
    
    return 0;
}
