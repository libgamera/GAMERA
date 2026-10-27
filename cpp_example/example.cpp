#include <iostream>
#include "Radiation.h"
#include "Utils.h"

int main() {
    double t = 2.7;
    double energy_density = 0.25 * eV_to_erg;

    Radiation radiation;

    radiation.AddThermalTargetPhotons(t, energy_density);
    auto photons = radiation.GetTargetPhotons();
    for (const auto& elems: photons) {
        for (double val: elems) {
            std::cout << val << ' ';
        }
        std::cout << '\n';
    }

    return 0;
}
