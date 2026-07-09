import numpy as np

def calculate_projected_sigma(angle_deg=75.0, L_cm=1.0, N_prim=34, use_heavy_tail=False):
    """
    Calculates the spatial resolution (sigma) projected onto the x-axis.
    
    Parameters:
        angle_deg (float): Angle between the track and the x-axis (0 = parallel to x).
        L_cm (float): The full 3D length of the track segment.
        N_prim (int): Mean number of primary clusters along the full length.
        use_heavy_tail (bool): If True, includes rare high-energy delta rays (n up to ~400).
                               If False, calculates the 'Core' resolution (typical events).
    """
    
    # --- 1. Calculate Geometric Projection ---
    # Convert angle to radians
    theta_rad = np.radians(angle_deg)
    
    # Calculate projected length on x-axis
    Lx_cm = L_cm * np.cos(theta_rad)
    
    # --- 2. Calculate Gas Statistics (N_eff) ---
    # Fischle et al. probabilities for Argon (Discrete part)
    fischle_data = {
        1: 0.656, 2: 0.150, 3: 0.064, 4: 0.035, 5: 0.022,
        6: 0.015, 7: 0.011, 8: 0.008, 9: 0.006
    }
    
    mean_n = 0.0
    mean_sq_n = 0.0
    prob_sum = 0.0
    
    # Discrete sums
    for n, p in fischle_data.items():
        mean_n += n * p
        mean_sq_n += (n**2) * p
        prob_sum += p
        
    # Heavy Tail Calculation (Optional)
    if use_heavy_tail:
        # Extend from n=10 up to n_max (approx 10 keV cutoff -> 384 electrons)
        n_max = 384
        remaining_prob = 1.0 - prob_sum
        
        # Normalization C for 1/n^2 distribution
        sum_inv_sq = sum(1.0/(n**2) for n in range(10, n_max + 1))
        C = remaining_prob / sum_inv_sq
        
        for n in range(10, n_max + 1):
            p = C / (n**2)
            mean_n += n * p
            mean_sq_n += (n**2) * p
            
    else:
        # Normalize the discrete part to 100% if ignoring tail (approximating the Core)
        # This effectively assumes rare delta-rays are rejected by the tracking algorithm
        normalization = 1.0 / prob_sum
        mean_n *= normalization
        mean_sq_n *= normalization

    # --- 3. Calculate Resolution ---
    
    # Calculate N_eff (Effective Electrons)
    # The statistical power is reduced by the fluctuation factor <n>^2 / <n^2>
    ratio = (mean_n**2) / mean_sq_n
    N_eff = N_prim * ratio
    
    # Final Sigma Calculation
    # Note: We use Lx_cm (projected length) but N_eff (total statistics from full track)
    sigma_cm = Lx_cm / np.sqrt(12 * N_eff)
    sigma_microns = sigma_cm * 10000

    # --- 4. Print Results ---
    mode = "Full RMS (Heavy Tail)" if use_heavy_tail else "Gaussian Core (Truncated)"
    
    print(f"--- Configuration: {angle_deg} degrees, Mode: {mode} ---")
    print(f"Full Length (L):      {L_cm:.4f} cm")
    print(f"Projected Length (Lx):{Lx_cm:.4f} cm (reduced by cos(75))")
    print(f"Stats Moments:        <n>={mean_n:.2f}, <n^2>={mean_sq_n:.2f}")
    print(f"N_eff (Effective e-): {N_eff:.2f} (out of {N_prim} primaries)")
    print(f"-------------------------------------------")
    print(f"CALCULATED SIGMA_X:   {sigma_microns:.1f} microns")
    print(f"-------------------------------------------\n")

# --- Run Scenarios ---

# Scenario 1: What matches your simulation (The "Core" fit)
calculate_projected_sigma(angle_deg=75.0, use_heavy_tail=False)

# Scenario 2: The statistical worst-case (including rare delta rays)
calculate_projected_sigma(angle_deg=75.0, use_heavy_tail=True)