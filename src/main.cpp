#include <iostream>
#include "physics_params.h"
int main() {
    reactor::physicsparameters p = reactor::createDefaultPhysics();
    std::cout <<"beta  =" << p.beta << "\n";
    std::cout<<"Pnom   ="<<p.Pnom<<"\n";
    std::cout << "Hello, reactor simulator\n";
    return 0;
}