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
    int min = 0;
    int max = 21;
    int N = 20;
    uniform_int_distribution<int> uni(min,max); // guaranteed unbiased
    
    cout << "Initial random distribution:\n";
    vector<int> v(N);    // declares a vector of integers
    for(int i=0; i<N; i++){
        auto random_integer = uni(rng);
        cout << random_integer << ' ';
        v[i] = random_integer;
    }
    
    
    cout << "\n\nTotal energy = ";
    int E = 0;
    for(int i=0; i<N; i++){
        for(int j=i+1; j<N; j++){
            if(v[i] == v[j]){
                E = E-1;
            }
        }
    }
    cout << E;
    
    
    
    
    cout << '\n';
    return 0;
}
