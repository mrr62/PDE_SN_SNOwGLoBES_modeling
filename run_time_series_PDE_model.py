from PDE_model import MySupernovaModel
from data_PDE_loader import load_neutrino_data
from snewpy.neutrino import Flavor
from astropy import units as u
import numpy as np

# Load luminosity data (in 10^51 erg/s)
time_nu_e, lum_nu_e = load_neutrino_data("/home/mrr62/miniconda3/lib/python3.11/site-packages/snewpy/PDE_SN_files/IP_ENL_QCD.csv", unit=u.erg/u.s, scale_factor=1e51)
time_nu_ebar, lum_nu_ebar = load_neutrino_data("/home/mrr62/miniconda3/lib/python3.11/site-packages/snewpy/PDE_SN_files/IP_EA-NL_QCD.csv", unit=u.erg/u.s, scale_factor=1e51)
time_nu_x, lum_nu_x = load_neutrino_data("/home/mrr62/miniconda3/lib/python3.11/site-packages/snewpy/PDE_SN_files/IP_m_tNL_QCD.csv", unit=u.erg/u.s, scale_factor=1e51)

# Load mean energy data (in MeV)
time_nu_e_energy, E_nu_e = load_neutrino_data("/home/mrr62/miniconda3/lib/python3.11/site-packages/snewpy/PDE_SN_files/IP_ENE_QCD.csv", unit=u.MeV)  #nu_e = electron neutrino
time_nu_ebar_energy, E_nu_ebar = load_neutrino_data("/home/mrr62/miniconda3/lib/python3.11/site-packages/snewpy/PDE_SN_files/IP_EA-NE_QCD.csv", unit=u.MeV) #nu_ebar = electron anti-neutrino
time_nu_x_energy, E_nu_x = load_neutrino_data("/home/mrr62/miniconda3/lib/python3.11/site-packages/snewpy/PDE_SN_files/IP_m_tNE_QCD.csv", unit=u.MeV) #nu_x = mu/tau neutrino

# Ensure time arrays are the same
if not (np.array_equal(time_nu_e, time_nu_e_energy) and 
        np.array_equal(time_nu_ebar, time_nu_ebar_energy) and
        np.array_equal(time_nu_x, time_nu_x_energy)):
    raise ValueError("Time arrays for luminosity and mean energy do not match:(")

# Organize data into dictionaries
luminosity = {
    Flavor.NU_E: lum_nu_e,
    Flavor.NU_E_BAR: lum_nu_ebar,
    Flavor.NU_X: lum_nu_x,
    Flavor.NU_X_BAR: lum_nu_x
}

mean_energy = {
    Flavor.NU_E: E_nu_e,
    Flavor.NU_E_BAR: E_nu_ebar,
    Flavor.NU_X: E_nu_x,
    Flavor.NU_X_BAR: E_nu_x
}

# Define pinching parameter alpha (2.23 for all flavors)
alpha_value = 2.23  # Single constant value for pinching parameter
alpha = {
    Flavor.NU_E: np.full(len(time_nu_e), alpha_value),
    Flavor.NU_E_BAR: np.full(len(time_nu_ebar), alpha_value),
    Flavor.NU_X: np.full(len(time_nu_x), alpha_value),
    Flavor.NU_X_BAR: np.full(len(time_nu_x), alpha_value)
}  #made into array to match time array


# Initialize and run model
metadata = {"progenitor_mass": "50 Msol", "eos": "LS220"} 
my_model = MySupernovaModel(time_nu_e, luminosity, mean_energy, alpha, metadata)

# Energy grid and flavors
E_eval = np.linspace(1, 100, 100) * u.MeV
flavors = [Flavor.NU_E, Flavor.NU_E_BAR, Flavor.NU_X, Flavor.NU_X_BAR]
t_grid = time_nu_e

# Output file
output_filename = "PDE_SN_model_timeseries.dat"

# Save full time series in SNOwGLoBES format
with open(output_filename, "w") as f:
    for t in t_grid:
        spectra = my_model.get_initial_spectra(t, E_eval, flavors=flavors)
        f.write(f"# Time = {t.to_value(u.s):.6f} s\n")
        f.write("# Energy (MeV) " + " ".join([f"{flavor}" for flavor in flavors]) + "\n")
        for i in range(len(E_eval)):
            f.write(f"{E_eval[i].value:.6f} " +
                    " ".join(f"{spectra[flavor][i]:.6e}" for flavor in flavors) + "\n")

print(f"Saved full time series neutrino flux data to {output_filename}")
