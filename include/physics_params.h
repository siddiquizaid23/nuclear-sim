#pragma once

namespace reactor {

struct physicsparameters {
    double beta;                    // delayed neutron fraction
    double lambda_delayed;          // precursor decay constant (1/s)
    double prompt_generation_time;  // neutron generation time (s)
    double alphaf;                  // fuel temperature coefficient (1/K)
    double gamma_decay;             // decay heat fraction
    double lambda_decay;            // decay heat decay constant (1/s)
    double Tin;                     // coolant inlet temperature (K)
    double Cf;                      // fuel thermal capacitance (J/K)
    double Cclad;                   // cladding thermal capacitance (J/K)
    double Cc;                      // coolant thermal capacitance (J/K)
    double hAgap;                   // gap conductance x area (W/K)
    double hAconvnom;               // nominal convective coefficient x area (W/K)
    double mnom;                    // nominal mass flow rate (kg/s)
    double Cp;                      // coolant specific heat (J/kg/K)
    double Pnom;                    // nominal thermal power (W)
    double flowexp;                 // flow exponent in heat transfer
    double Tref;                    // reference fuel temperature (K), set at startup
};

physicsparameters createDefaultPhysics();

}  // namespace reactor