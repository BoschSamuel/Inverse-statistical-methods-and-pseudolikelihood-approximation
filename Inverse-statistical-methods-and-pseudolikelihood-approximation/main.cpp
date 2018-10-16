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

using namespace std;


int main(int argc, const char * argv[]) {
    random_device rd;     // only used once to initialise (seed) engine
    mt19937 rng(rd());    // random-number engine used (Mersenne-Twister in this case)
    int min_spin = 0;  // min value spin can take
    int max_spin = 21; // max value spin can take
    int N = 50; //Number of atoms in molecule
    double T = 300; // This is the temperature of the experiment
    // double K_b = 1.38064852e-23; // Boltzman constant
    double K_b = 0.008;
    uniform_int_distribution<int> random_spin(min_spin,max_spin); // definition of random function
    uniform_int_distribution<int> random_atom(0,N-1); // definition of random function
    
    // Random initialization of spins
    //cout << "Initial random distribution:\n";
    vector<int> v(N);    // declares a vector of integers
    for(int i=0; i<N; i++){
        auto random_integer = random_spin(rng);
        cout << random_integer << ' ';
        v[i] = random_integer;
    }
    
    // Calculation of the energy using Pott's model (H = -J*sum(Kronecker_delta(i,j))
    //cout << "\nTotal energy = ";
    int E = 0;
    for(int i=0; i<N; i++){
        for(int j=i+1; j<N; j++){
            if(v[i] == v[j]){
                E = E-1;
            }
        }
    }
    cout << E << "\n\n";

    //cout << "\nStarting random iterations...\n";
    int max_number_of_interations = 10000;
    for(int iteration_number=0; iteration_number<max_number_of_interations; iteration_number++){
        auto atom_number = random_atom(rng); //Pick random atom for changing the spin
        int old_spin = v[atom_number]; //Saving old spin in case we still want to use it
        //cout << "Atom is = " << atom_number << " Spin was = " << old_spin;
        int E_old = E;
        v[atom_number] = random_spin(rng); //Pick random new spin for selected atom
        //cout << " and now is = " << v[atom_number] << '\n';
        if(v[atom_number] == old_spin){
            continue; // If the spin didn't change, we do nothing
        }
        // Calculation of the new energy
        E = 0;
        for(int i=0; i<N; i++){
            for(int j=i+1; j<N; j++){
                if(v[i] == v[j]){
                    E = E-1;
                }
            }
        }
        //cout << "Energy was = " << E_old << " and now is = " << E << '\n';
        if(E > E_old){
            //cout << "Energy want to increase\n";
            double uniform_random_0_1 = (double)rand()/(double)RAND_MAX;
            double prob = exp(-(E-E_old)/(K_b*T));
            if (prob < uniform_random_0_1){
                // With some small probability, we revert the change (Metropolis algorithm)
                E = E_old;
                v[atom_number] = old_spin;
            }
            //cout << prob << ' ' << uniform_random_0_1 << '\n';
        }else{
            //cout << "Energy decreased\n";
        }
        //cout << "Energy now is = " << E << "\n\n";
        cout << E << '\n';
    }
    
    
    cout << "Final spin configuration:\n";
    for(int i=0; i<N; i++){
        cout << v[i] << ' ';
    }
    
    
    cout << '\n';
    return 0;
}
