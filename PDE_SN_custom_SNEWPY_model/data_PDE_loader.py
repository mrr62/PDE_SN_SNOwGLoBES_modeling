import pandas as pd
from astropy import units as u

def load_neutrino_data(file_path, unit=1, scale_factor=1):
    """
    Load neutrino data from a CSV file and apply unit corrections.

    Parameters
    ----------
    file_path : str
        Path to the CSV file.
    unit : astropy.Unit
        Physical unit of the values (e.g., erg/s for luminosity, MeV for energy).
    scale_factor : float, optional
        Factor to multiply values (e.g., 1e51 for luminosity in 10^51 erg/s).

    Returns
    -------
    time : astropy.Quantity
        Time values in seconds.
    values : astropy.Quantity
        Corresponding values in the specified unit.
    """
    df = pd.read_csv(file_path, names=["time", "value"], skiprows=1)  # Assuming no header

    # Convert time from milliseconds (ms) to seconds (s)
    time = df["time"].values * 0.001 * u.s #attach units so array (hopefully) works
    #time = time  # Convert to seconds

    # Apply scale factor to values (e.g., 10^51 for luminosity)
    values = df["value"].values * unit * scale_factor

    return time, values
