#include <cstdio>
#include "reactor_model.h"
using namespace reactor;

// Run the reactor from steady state with a constant reactivity insertion.
void run(const char* label, double rho, double duration_s) {
    physicsparameters p = createDefaultPhysics();
    State s = init_critical(p);

    const double dt = 0.01;
    const int steps = static_cast<int>(duration_s / dt);
    const int print_every = static_cast<int>(30.0 / dt);  // every 30 s

    std::printf("\n%s (rho = %+.4f)\n", label, rho);
    std::printf("%8s %10s %10s %10s\n", "t (s)", "P", "Tf (K)", "Tc (K)");

    for (int i = 0; i <= steps; ++i) {
        if (i % print_every == 0) {
            std::printf("%8.1f %10.4f %10.2f %10.2f\n", i * dt, s.P, s.Tf, s.Tc);
        }
        s = rk4_step(s, rho, p.mnom, dt, p);
    }
}

int main() {
    run("Positive insertion", +0.001, 600.0);
    run("Negative insertion", -0.001, 600.0);

    return 0;
}