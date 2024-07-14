#include "pintr/pintr.hpp"
#include "pintr/utils/args.hpp"

#include <filesystem>
#include <spdlog/spdlog.h>

int main(int argc, char **argv) {
    // Call the setup function to configure logging
    pintr::setup();

    // Create a database instance
    std::filesystem::path databasedir;
    bool download_missing = false;

    for (int i = 1; i < argc; ++i) {
        bool found = pintr::args::parse_download_missing(i, argc, argv, download_missing);
        if (!found) {
            pintr::args::parse_database(i, argc, argv, databasedir);
        }
    }

    pintr::Database database(download_missing, false, databasedir);

    // Create a basis
    auto basis = pintr::BasisAtomCreator<double>()
                     .set_species("Rb")
                     .restrict_quantum_number_n(58, 62)
                     .restrict_quantum_number_l(0, 3)
                     .create(database);

    SPDLOG_INFO("Number of basis states: {}", basis->get_number_of_states());

    // Create a system
    auto system = pintr::SystemAtom<double>(basis);
    system.set_electric_field({0, 0, 1});
    system.set_magnetic_field({0, 0, 1});
    system.enable_diamagnetism(true);

    return 0;
}
