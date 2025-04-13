# Python file to calculate event rates over time instead of energy
# Requires input .dat file with format: [time bin start, time bin end, nue_bar fluence, nue fluence, nux fluence] 
# Fluence is in [neutrinos/cm²]

import pandas as pd
import numpy as np
import os

# Paths
fluence_file = "/home/mrr62/SN_stuff/globes-3.2.18/snowglobes/fluxes/ts_combined_PDE_fluence.dat" #path to fluence file
cross_section_dir = "/home/mrr62/SN_stuff/globes-3.2.18/snowglobes/xscns" #path to SNOwGLoBES cross-section directory
output_dir = "/home/mrr62/SN_stuff/globes-3.2.18/snowglobes/out/ts_water"

# Detector config for wc100kt30prct - can be edited for any supported detector config, but must also edit allowed interaction channels and cross-sections, as well as the calculation for target count 
detector_name = "wc100kt30prct"
detector_mass_kg = 100_000 * 1e3  # 100 kton to kg
normalization = 0.1111  # from detector_configurations.dat

# Load fluence data
df = pd.read_csv(fluence_file, header=None, names=["start_time", "end_time", "ea_n", "e", "mu_tau"])

# Flavor column map (for fluence file)
fluence_map = {
    "nuebar": "ea_n",
    "nue": "e",
    "numu": "mu_tau",
    "numubar": "mu_tau",
    "nutau": "mu_tau",
    "nutaubar": "mu_tau"
}

# Track how many files we write
count = 0

# Loop over cross section files
for fname in sorted(os.listdir(cross_section_dir)):
    if not fname.endswith(".dat") or not fname.startswith("xs_"):
        continue

    parts = fname.replace(".dat", "").split("_")
    if len(parts) < 2:
        continue

    # Special case: IBD
    if parts[1] == "ibd":
        flavor = "nuebar"         # IBD is always nuebar
        interaction = "ibd"
        target = "ibd"
        fluence_col = fluence_map.get(flavor)
        if "ibd2" in fname or "ibd_he" in fname:
            continue  # Skip these

    # Special case: Neutral Current (NC)
    elif parts[1] == "nc" and len(parts) >= 4:
        flavor = parts[2]         # xs_nc_nutau_O16 → nutau
        target = parts[3]         # O16
        interaction = "_".join(parts[1:])  # nc_nutau_O16
        fluence_col = fluence_map.get(flavor)

    # Regular case
    else:
        flavor = parts[1]
        target = parts[-1]
        interaction = "_".join(parts[2:])
        fluence_col = fluence_map.get(flavor)

    # Skip if fluence not available (e.g., sterile flavors)
    if fluence_col is None:
        continue

    # Filter by allowed targets for your detector
    allowed_targets = ["e", "O16", "ibd"]
    if target not in allowed_targets:
        continue


    # Read cross section file — for multi-column files
    xs_file = os.path.join(cross_section_dir, fname)
    xs_data = pd.read_csv(xs_file, delim_whitespace=True, comment="#", header=None)
    xs_data.columns = ["logE", "nue", "numu", "nutau", "nuebar", "numubar", "nutaubar"]

    # Convert log10(E/GeV) → E [MeV]
    xs_data["E_MeV"] = 10 ** xs_data["logE"] * 1000

    # Use correct flavor column for sigma
    # Use correct flavor column for sigma
    xs_data["sigma_MeV"] = xs_data[flavor] *1e-38 / 1000  # Convert [10^-38 cm²/GeV] → [cm²/MeV]
    xs_data = xs_data[xs_data["sigma_MeV"] > 0]  # Filter physical


    if xs_data.empty:
        print(f"⚠️ Warning: No valid cross section data in {fname}. Skipping...")
        continue

    # Compute event rates per time bin
    # Compute event counts per time bin
    # Constants for target count
    NA = 6.022e23  # Avogadro's number [mol^-1]
    molar_mass_water = 18.01528e-3  # kg/mol

    # Determine number of targets per molecule
    if target == "e":
        targets_per_molecule = 10  # 10 electrons per water molecule (8 from O, 1 each from H)
    elif target == "ibd":
        targets_per_molecule = 2  # 2 free protons per water molecule (from 2 H atoms)
    else:
        targets_per_molecule = 1  # 1 O16 nucleus per water molecule

    # Total number of target particles for this interaction
    num_targets = (detector_mass_kg / molar_mass_water) * targets_per_molecule * NA

    # Compute event counts per time bin
    results = []
    for _, row in df.iterrows():
        fluence = row[fluence_col]  # [neutrinos / cm²]
        if fluence_col == "mu_tau":
            fluence /= 4.0  # Evenly distribute among 4 flavors

        # Read approximate energy spectrum shape (assume normalized shape) since we only know total fluence and
        # not the spectrum/shape of the fluence per time bin
        # Therefore use a synthetic spectrum shape — assume constant dN/dE shape
        # --> weight σ(E) with an average neutrino energy spectrum (flat for now)
        # and properly weight: integral(Φ(E) * σ(E) dE)

        # For realistic result compute:
        # fluence_spectrum(E) = fluence_total * shape(E)  [assume flat shape]

        # Normalize a flat spectrum over the energy range
        shape = np.ones_like(xs_data["E_MeV"])
        shape /= np.trapz(shape, xs_data["E_MeV"])

        # Multiply σ(E) * shape(E) and integrate → average σ
        sigma_avg = np.trapz(xs_data["sigma_MeV"] * shape, xs_data["E_MeV"])  # [cm²]
        
        #OPTIONAL DEBUGGING FEATURE — checks cross section values and fluence values 
        #print(f"Cross section average for {fname}: {sigma_avg:.3e}")
        #print(f"Fluence sample for {flavor} at time bin {row['start_time']}-{row['end_time']}: {fluence:.3e}")

        # Expected number of events in that time bin
        event_count = fluence * sigma_avg * num_targets  # unitless (just a count)

        #OPTIONAL DEBUGGING FEATURE — checks integrated cross sections
        #print(f"Cross section integral for {fname}: {sigma_avg:.3e}")

        # Expected number of events in that time bin
        event_count = fluence * sigma_avg * num_targets  # unitless (just a count)

        results.append((row["start_time"], row["end_time"], event_count))


    # Write output file
    outname = f"PDE_{flavor}_{interaction}_{detector_name}_time_binned.dat"
    outfile = os.path.join(output_dir, outname)
    with open(outfile, "w") as f:
        for start, end, rate in results:
            f.write(f"{start},{end},{rate:.6e}\n")

    print(f"✅ Saved: {outname}")
    count += 1

# Summary
print(f"\n🎉 Done! Saved {count} time-binned event rate files.")
