#include "physics_params.h"

namespace reactor {

physicsparameters createDefaultPhysics() {
    physicsparameters params;
    params.beta = 0.0065;
    params.lambda_delayed = 0.09;
    params.prompt_generation_time = 0.001;
    params.alphaf = -5.0e-6;
    params.gamma_decay = 0.06;
    params.lambda_decay = 0.1;
    params.Tin = 290.0;
    params.Cf = 2.65e6;
    params.Cclad = 1.15e6;
    params.Cc = 3.1e7;
    params.hAgap = 5.3e5;
    params.hAconvnom = 1.77e6;
    params.mnom = 505.0;
    params.Cp = 4200.0;
    params.Pnom = 100.0e6;
    params.flowexp = 0.8;
    params.Tref = 0.0;  // set later by init_critical
    return params;
}

}  // namespace reactor