#pragma once
#include "physics_params.h"

namespace reactor {

struct State {
    double P;      // relative power (1.0 = steady state)
    double C;      // delayed neutron precursors
    double Tf;     // average fuel temperature (K)
    double D;      // decay heat power (relative)
    double Tclad;  // average cladding temperature (K)
    double Tc;     // average coolant temperature (K)
};

struct derivatives {
    double dpdt;
    double dcdt;
    double dtfdt;
    double dddt;
    double dtcladdt;
    double dtcdt;
};

// Steady state at full power. Also sets params.Tref.
State init_critical(physicsparameters& params);

// Rates of change of every state variable.
derivatives compute_deri(const State& state, double rho_user, double m_dot,
                         const physicsparameters& params);

State rk4_step (const State& state , double rho_user , double m_dot, double dt,
const physicsparameters& params);

}  // namespace reactor