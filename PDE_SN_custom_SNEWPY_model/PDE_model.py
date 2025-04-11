import numpy as np
from astropy import units as u
from scipy.special import loggamma
from snewpy.models.base import SupernovaModel
from snewpy.neutrino import Flavor

class MySupernovaModel(SupernovaModel):
    """Custom Supernova Model based on luminosity and mean energy data."""

    def __init__(self, time, luminosity, mean_energy, alpha, metadata):
        """
        Initialize the supernova model.

        Parameters
        ----------
        time : ndarray of astropy.Quantity
            Time points where the model flux is defined.
        luminosity : dict of astropy.Quantity
            Dictionary of luminosities (erg/s) for each neutrino flavor.
        mean_energy : dict of astropy.Quantity
            Dictionary of mean energies (MeV) for each neutrino flavor.
        alpha : dict of float
            Pinching parameter for each neutrino flavor.
        metadata : dict
            Dictionary of model parameters.
        """
        self.luminosity = luminosity
        self.mean_energy = mean_energy
        self.alpha = alpha

        super().__init__(time, metadata)

    def get_initial_spectra(self, t, E, flavors=Flavor):
        """
        Compute the initial neutrino spectra at the source.

        Parameters
        ----------
        t : astropy.Quantity
            Time to evaluate initial spectra.
        E : astropy.Quantity or ndarray of astropy.Quantity
            Energies to evaluate the initial spectra.
        flavors: iterable of snewpy.neutrino.Flavor
            Return spectra for these flavors only.

        Returns
        -------
        initial_spectra : dict
            Dictionary of model spectra, keyed by neutrino flavor.
        """
        t = u.Quantity(t, ndmin=1)
        E = u.Quantity(E, ndmin=1)    # Convert energy to erg
        E_values = E                  # Convert to a unitless array for calculations
        E_values = np.expand_dims(E_values, axis=0)  # Reshape for broadcasting

        initial_spectra = {}

        # Ensure input time is in the same unit as the model time grid
        t = t

        for flavor in flavors:
            # Interpolate L, E_a, and alpha over time
            L_values = np.interp(t, self.time, self.luminosity[flavor])
    
            # Convert mean energy (originally in MeV) to erg
            mean_energy_vals = self.mean_energy[flavor]  # This should be an astropy Quantity array
            Ea_values = np.interp(t, self.time, mean_energy_vals.to(u.MeV).value) * u.MeV.to(u.erg)



            a_values = np.interp(t, self.time, self.alpha[flavor])

            # Reshape to match the energy dimension
            L_values = np.expand_dims(L_values, axis=1)
            Ea_values = np.expand_dims(Ea_values, axis=1)
            a_values = np.expand_dims(a_values, axis=1)

            # Ensure E_values is unitless
            #E_values = E.to_value(u.erg)  # Convert to erg and strip units
            #E_values = np.expand_dims(E_values, axis=0)  # Reshape for broadcasting

            # Ensure L_values, Ea_values, and E_values are dimensionless
            # Strip units and get dimensionless float arrays
            L_norm = L_values.to_value(u.erg / u.s)
            Ea_norm = Ea_values #this one might already be a float?? idk
            E_norm = E.to_value(u.erg)




            # Convert everything to numpy float arrays
           # L_norm = L_norm.astype(float)
           # Ea_norm = Ea_norm.astype(float)
           # E_norm = E_norm.astype(float)
           # a_values = a_values.astype(float)

            # Compute the spectrum using the pinched thermal distribution
            spectrum = np.exp(
            np.log(L_norm) - (2 + a_values) * np.log(Ea_norm) + (1 + a_values) * np.log(1 + a_values)
            - loggamma(1 + a_values) + a_values * np.log(E_norm) - (1 + a_values) * (E_norm / Ea_norm)
    )

            # Remove bad values (NaNs, zeros)
            spectrum[np.isnan(spectrum)] = 0
            spectrum[:, E.to_value(u.erg) == 0] = 0
            spectrum = np.squeeze(spectrum)  # Remove single-dimensional entries

            initial_spectra[flavor] = spectrum  # Now completely dimensionless


        return initial_spectra

