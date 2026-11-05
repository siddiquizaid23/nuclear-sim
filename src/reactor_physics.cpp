#include "physics_params.h"
#include"reactor_model.h"
#include<cmath>
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
    params.Tin = 290.0;
    params.Tref = 0.0;  // set later by init_critical
    return params;
}
State init_critical(physicsparameters& params){
  State state{};
  double Q = (1.0 + params.gamma_decay) * params.Pnom;
  state.P = 1.0;
  state.C = params.beta / (params.lambda_delayed * params.prompt_generation_time);
  state.Tc = params.Tin + Q / (params.mnom * params.Cp);
  state.Tclad = state.Tc + Q / params.hAconvnom;
  state.Tf = state.Tclad + Q / params.hAgap;
  state.D = params.gamma_decay;
  params.Tref = state.Tf;
  return state;
}

derivatives compute_deri(const State& state, double rho_user,double m_dot, const physicsparameters& params ){

    derivatives d;

    double rho = rho_user + params.alphaf *(state.Tf - params.Tref);

    double flow_factor = std::pow(m_dot / params.mnom , params.flowexp);
    double hAconv = flow_factor * params.hAconvnom;


    d.dpdt = ((rho - params.beta) / params.prompt_generation_time) * state.P
             + params.lambda_delayed * state.C;

    d.dcdt = (params.beta / params.prompt_generation_time) * state.P
             - params.lambda_delayed * state.C;

    d.dtfdt = ((state.P + state.D) * params.Pnom
               - params.hAgap * (state.Tf - state.Tclad)) / params.Cf;

    d.dddt = params.lambda_decay * (state.P * params.gamma_decay - state.D);

    d.dtcladdt = (params.hAgap * (state.Tf - state.Tclad)
                  - hAconv * (state.Tclad - state.Tc)) / params.Cclad;

    d.dtcdt = (hAconv * (state.Tclad - state.Tc)
               + m_dot * params.Cp * (params.Tin - state.Tc)) / params.Cc;

    return d;
}

State offset(const State& s, const derivatives& d, double h) {
    State out;
    out.P      = s.P      + h * d.dpdt;
    out.C      = s.C      + h * d.dcdt;
    out.Tf     = s.Tf     + h * d.dtfdt;
    out.D      = s.D      + h * d.dddt;
    out.Tclad  = s.Tclad  + h * d.dtcladdt;
    out.Tc     = s.Tc     + h * d.dtcdt;
    return out;

}

State rk4_step(const State& state, double rho_user, double m_dot , double dt, const physicsparameters& params){
    derivatives k1 = compute_deri(state ,rho_user,m_dot,params);
    derivatives k2 = compute_deri(offset(state,k1 ,dt/ 2),rho_user,m_dot,params);
    derivatives k3 = compute_deri(offset(state,k2,dt/2),rho_user,m_dot,params);
    derivatives k4 = compute_deri(offset(state,k3,dt),rho_user,m_dot,params);

    derivatives avg;
      avg.dpdt = (k1.dpdt + 2*k2.dpdt + 2*k3.dpdt + k4.dpdt) / 6.0;
      avg.dcdt = (k1.dcdt +2*k2.dcdt+2*k3.dcdt+k4.dcdt)/6.0;
      avg.dtfdt = (k1.dtfdt + 2*k2.dtfdt+2*k3.dtfdt+k4.dtfdt)/6.0;
      avg.dddt = (k1.dddt + 2*k2.dddt + 2*k3.dddt+k4.dddt)/6.0;
      avg.dtcladdt =(k1.dtcladdt+ 2*k2.dtcladdt+2*k3.dtcladdt+k4.dtcladdt)/6.0;
      avg.dtcdt =(k1.dtcdt + 2*k2.dtcdt+2*k3.dtcdt+k4.dtcdt)/6.0;

      return offset(state,avg,dt);
}
}  // namespace reactor