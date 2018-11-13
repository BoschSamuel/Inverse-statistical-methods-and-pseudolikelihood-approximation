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


int main(int argc, const char * argv[]){
    random_device rd;  // only used once to initialise (seed) engine
    mt19937 rng(rd()); // random-number engine used (Mersenne-Twister in this case)
    int min_spin = 0;  // min value spin can take
    int max_spin = 20; // max value spin can take
    int N = 25;        //Number of atoms in molecule
    double T = 1;  // This is the temperature of the experiment
    // double K_b = 1.38064852e-23; // Boltzman constant
    double K_b = 1;
    uniform_int_distribution<int> random_spin(min_spin,max_spin); // definition of random function
    uniform_int_distribution<int> random_atom(0,N-1); // definition of random function

    array<std::array<double, 525>, 525> J;
    ifstream myfile;
    myfile.open("/Users/samuelbosch/OneDrive/Faks/EPFL_M1/Computer_simulation/Project/Inverse-statistical-methods-and-pseudolikelihood-approximation/J.txt");
    if(!myfile){ //Testing if file is actually open.
        cout << "Error opening file" << endl;
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
            E = E + J[21*i+v[i]][21*j+v[j]]; // Is this the right formula for calculating the energy?
        }
    }
    cout << "\nEnergy(initial)" << " = " << E << "\n\n";

    //cout << "\nStarting random iterations...\n";
    int energy_change_counter = 0;
    int max_number_of_interations = 100000;
    vector<double> Energy(max_number_of_interations);
    int iteration_number = 0;
    for(; iteration_number<max_number_of_interations; iteration_number++){
        Energy[energy_change_counter] = E;
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
            E = E + J[21*i+v[i]][21*atom_number+v[atom_number]]; // Here we add the new energy
            E = E - J[21*i+v[i]][21*atom_number+old_spin]; // Here we substract the old energy
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
        energy_change_counter++;
    }

    cout << "Final spin configuration:\n";
    for(int i=0; i<N; i++){
        cout << v[i] << ' ';
    }
    cout << '\n';
    for(int iteration_number=0; iteration_number<energy_change_counter; iteration_number++){
        // cout << Energy[iteration_number] << '\n';
    }
    
    cout << endl;
    //cout << "Energy(" << iteration_number << ")" << " = " << E << "\n";
    cout << '\n';

    double mean = accumulate(begin(Energy), end(Energy), 0.0) / Energy.size();
    //cout << "Average Energy = " << mean << '\n';
    
    
    // Writing the Energy vs step number into a .txt file
    // The specific path was need, as it is otherwise saved in the xcode hidden folder
    ofstream energy_file;
    energy_file.open("/Users/samuelbosch/OneDrive/Faks/EPFL_M1/Computer_simulation/Project/Inverse-statistical-methods-and-pseudolikelihood-approximation/Energy_vs_time.txt");
    if (myfile.is_open()) { cout << "File 'Energy_vs_time.txt' is open\n\n"; }
    for(int i=0; i<energy_change_counter; i++){
        energy_file << Energy[i] << '\n';
    }
    energy_file.close();
    
    
    
// Autocorrelation function
    double autocorrelation_fraction = 0.03; // Through what fraction of the data do you want the autocorr. function to go?
    vector<double> autocorrelation((int)(autocorrelation_fraction*energy_change_counter));
    for(int i=0; i<int(autocorrelation_fraction*energy_change_counter); i++){
        autocorrelation[i] = 0;
        for(int j=0; j<energy_change_counter-i; j++){
            autocorrelation[i] += Energy[j]*Energy[j+i] - mean*mean;
        }
        if (energy_change_counter-i>0){
            autocorrelation[i] /= (energy_change_counter-i);
        }
    }
    double normalisation_factor = autocorrelation[0];
    for(int i=0; i<int(autocorrelation_fraction*energy_change_counter); i++){
        autocorrelation[i] /= normalisation_factor;
    }

    // activate the following lines if you want to get the autocorrelation printed
    // cout << "\n\n\n\n" << "autocorrelation:" << '\n';
    // for(int i=0; i<int(autocorrelation_fraction*energy_change_counter); i++){
        // cout << autocorrelation[i] << '\n';
    // }
    
    
    // The blocking method analysis
    int n = 10; //number of blocks
    vector<double> block_averages(n);
    double sum = 0;
    int k = int(energy_change_counter/n);
    int j = 0;
    for (int i=0; i<energy_change_counter; i++){
        sum += Energy[i];
        if (i==k){
            block_averages[j] = sum/(int(energy_change_counter/n));
            sum = 0;
            k += int(energy_change_counter/n);
            j++;
        }
    }
    cout << "Blocking method analysis with " << n << " blocks:\n";
    for (int i=0; i<n; i++){
        cout << block_averages[i] << '\n';
    }
    
    
    return 0;
}
