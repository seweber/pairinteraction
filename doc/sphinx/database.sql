-- Version 1.0

-- Unit System: Atomic Units
--
-- Base units:
-- - length: Bohr radius
-- - time: reduced Planck constant divided by Hartree energy
-- - temperature: Hartree energy divided by Boltzmann constant
-- - mass: electron mass
-- - current: elementary charge times Hartree energy divided by reduced Planck constant
--
-- Fundamental constants in this unit system:
-- - reduced Planck constant:   1
-- - Coulomb constant:          1
-- - Bohr magneton:             0.5
-- - elementary charge:         1
-- - electron mass:             1
--

BEGIN TRANSACTION;

-- Create a table of states
-- The energies are defined such that the mean ionization threshold is zero.
CREATE TABLE states (
    id INTEGER PRIMARY KEY,
    energy REAL NOT NULL, -- the energy of the state
    parity INTEGER NOT NULL, -- the parity of the state
    n INTEGER NOT NULL, -- the principal quantum number of the Rydberg electron (is zero if multi-channel)
    f REAL NOT NULL, -- the total momentum quantum number (equals exp_j if hyperfine splitting neglected)
    exp_nu REAL NOT NULL, -- the expectation value of effective principal quantum number of the Rydberg electron
    exp_l REAL NOT NULL, -- the expectation value of the orbital quantum number of the Rydberg electron
    exp_j REAL NOT NULL, -- the expectation value of the total angular quantum number of the Rydberg electron
    exp_s REAL NOT NULL, -- the expectation value of the total spin quantum number of all valence electrons
    std_nu REAL NOT NULL, -- the standard deviation of the effective principal quantum number of the Rydberg electron
    std_l REAL NOT NULL, -- the standard deviation of the orbital quantum number of the Rydberg electron
    std_j REAL NOT NULL, -- the standard deviation of the total angular quantum number of the Rydberg electron
    std_s REAL NOT NULL -- the standard deviation of the total spin quantum number of all valence electrons
);

CREATE TRIGGER check_energy
BEFORE INSERT ON states
BEGIN
    SELECT RAISE(FAIL, "id and energy must be inserted in increasing order")
    FROM states
    WHERE id > NEW.id OR energy > NEW.energy;
END;

-- Create a table of wigner 3j symbols multiplied by the phase factor (-1)^(f-m)
CREATE TABLE wigner (
    f_initial REAL, -- the total momentum quantum number of the initial state
    f_final REAL, -- the total momentum quantum number of the final state
    m_initial REAL,
    m_final REAL,
    kappa INTEGER,
    q INTEGER,
    val REAL NOT NULL,
    PRIMARY KEY (f_initial, f_final, m_initial, m_final, kappa, q)
);

-- Create a table of matrix elements of the reduced electric dipole operator e*r*sqrt(4*pi/3)*Y_1
CREATE TABLE matrix_elements_d (
    id_initial INTEGER,
    id_final INTEGER,
    val REAL NOT NULL,
    PRIMARY KEY (id_initial, id_final)
);

-- Create a table of matrix elements of the reduced electric quadrupole operator e*r^2*sqrt(4*pi/5)*Y_2
CREATE TABLE matrix_elements_q (
    id_initial INTEGER,
    id_final INTEGER,
    val REAL NOT NULL,
    PRIMARY KEY (id_initial, id_final)
);

-- Create a table of matrix elements of the reduced electric octupole operator e*r^3*sqrt(4*pi/7)*Y_3
CREATE TABLE matrix_elements_o (
    id_initial INTEGER,
    id_final INTEGER,
    val REAL NOT NULL,
    PRIMARY KEY (id_initial, id_final)
);

-- Create a table of matrix elements of the reduced magnetic dipole operator
-- \mu = - \mu_B / \hbar * (g_L * L + g_S * S) + \mu_N / \hbar g_I I,
-- taking into account all spins: Rydberg electron, core electron, and nucleus
CREATE TABLE matrix_elements_mu (
    id_initial INTEGER,
    id_final INTEGER,
    val REAL NOT NULL,
    PRIMARY KEY (id_initial, id_final)
);

-- Create a table of matrix elements of the reduced operator e*r^2*sqrt(4*pi)*Y_0, occurring in
-- the operator of the diamagnetic interaction
CREATE TABLE matrix_elements_dia (
    id_initial INTEGER,
    id_final INTEGER,
    val REAL NOT NULL,
    PRIMARY KEY (id_initial, id_final)
);

-- Create indices to speed up queries
CREATE INDEX idx_states_energy ON states (energy);
CREATE INDEX idx_wigner_f_initial ON wigner (f_initial);
CREATE INDEX idx_wigner_f_final ON wigner (f_final);
CREATE INDEX idx_matrix_elements_d_id_initial ON matrix_elements_d (id_initial);
CREATE INDEX idx_matrix_elements_d_id_final ON matrix_elements_d (id_final);
CREATE INDEX idx_matrix_elements_q_id_initial ON matrix_elements_q (id_initial);
CREATE INDEX idx_matrix_elements_q_id_final ON matrix_elements_q (id_final);
CREATE INDEX idx_matrix_elements_o_id_initial ON matrix_elements_o (id_initial);
CREATE INDEX idx_matrix_elements_o_id_final ON matrix_elements_o (id_final);
CREATE INDEX idx_matrix_elements_dia_id_initial ON matrix_elements_dia (id_initial);
CREATE INDEX idx_matrix_elements_dia_id_final ON matrix_elements_dia (id_final);
CREATE INDEX idx_matrix_elements_mu_id_initial ON matrix_elements_mu (id_initial);
CREATE INDEX idx_matrix_elements_mu_id_final ON matrix_elements_mu (id_final);

COMMIT;
