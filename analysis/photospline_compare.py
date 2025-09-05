import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import LogNorm
import numpy as np
import os
import photospline
from itertools import cycle
from scipy.interpolate import PchipInterpolator,LinearNDInterpolator

lines = ["--","-.",":"]
linecycler = cycle(lines)

def compare_splines(spline1_path, spline2_path, spline1_label, spline2_label,
                    energies, x_values, y_values,
                    output_dir):
    """
    Compare two spline files by taking ratios at different energy and x points
    """

    # Load splines
    spline1 = photospline.SplineTable(spline1_path)
    spline2 = photospline.SplineTable(spline2_path)

    # Create output directory
    os.makedirs(output_dir, exist_ok=True)


    # Create colormap
    cmap = mpl.colormaps["viridis"]

    # Plot 1: Ratio vs y for different energies and x values
    for energy in energies:
        print(f"Processing energy: {energy} GeV")

        # fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 10), sharex=True)
        # fig.subplots_adjust(hspace=0)

        # plot_anything = False

        # for i, x in enumerate(x_values):


        #     # Evaluate both splines
        #     spline1_vals = spline1.evaluate_simple([
        #         np.log10(energy), np.log10(x), np.log10(y_values)
        #     ])
        #     spline2_vals = spline2.evaluate_simple([
        #         np.log10(energy), np.log10(x), np.log10(y_values)
        #     ]) / 1e4  # Convert from m^2 to cm^2

        #     # Calculate ratio (avoid division by zero)
        #     ratio = np.where(spline2_vals > 0, spline1_vals / spline2_vals, np.nan)

        #     if np.all(np.isnan(ratio)):
        #         continue

        #     plot_anything = True
        #     color = cmap(i/len(x_values))

        #     # Plot individual splines
        #     ax1.plot(y_values, spline1_vals, color=color, ls='-', alpha=0.7)
        #     ax1.plot(y_values, spline2_vals, color=color, ls='--', alpha=0.7)

        #     # Plot ratio
        #     ax2.plot(y_values, ratio, color=color, label=f"x={x:.2e}")

        # if not plot_anything:
        #     plt.close(fig)
        #     continue

        # # Formatting
        # ax1.plot([], [], 'k-', label=spline1_label)
        # ax1.plot([], [], 'k--', label=spline2_label)
        # ax1.legend()
        # ax1.set_yscale('log')
        # ax1.set_ylabel(r'$d^2\sigma/dxdy$ [cm$^2$]')
        # ax1.set_title(f'Spline Comparison at $E_\\nu$ = {energy:.2e} GeV')
        # ax1.set_ylim(1e-60,1e-30)

        # ax2.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        # ax2.set_xlabel('Bjorken y')
        # ax2.set_ylabel(f'Ratio ({spline1_label}/{spline2_label})')
        # ax2.axhline(y=1, color='black', linestyle=':', alpha=0.5)
        # ax2.set_xscale('log')
        # ax2.set_ylim(0,2)

        # plt.tight_layout()
        # plt.savefig(f"{output_dir}/ratio_E_{energy:.2e}.pdf", dpi=150, bbox_inches='tight')
        # plt.close(fig)

        # Plot 2: Ratio heatmap as a funciton of x and y

        fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(8, 15), sharex=True)

        # Create meshgrid for evaluation
        x_grid = np.logspace(np.log10(x_values.min()), np.log10(x_values.max()), 100)
        y_grid = np.logspace(np.log10(y_values.min()), np.log10(y_values.max()), 100)

        X_mesh, Y_mesh = np.meshgrid(x_grid, y_grid)
        E_mesh = np.full_like(X_mesh, energy)

        Q2 = 2 * 0.938 * energy * X_mesh * Y_mesh  # GeV^2
        Q2_mask = Q2 > 1.0  # Below 1 GeV^2 might be unreliable
        print(Q2_mask.shape)
        print(sum(Q2_mask))

        # Evaluate splines on grid
        spline1_grid = spline1.evaluate_simple([
            np.log10(E_mesh.flatten()),
            np.log10(X_mesh.flatten()),
            np.log10(Y_mesh.flatten())
        ]).reshape(E_mesh.shape)

        spline2_grid = spline2.evaluate_simple([
            np.log10(E_mesh.flatten()),
            np.log10(X_mesh.flatten()),
            np.log10(Y_mesh.flatten())
        ]).reshape(E_mesh.shape) - 4  # Convert from m^2 to cm^2
        #spline2_grid = np.where(np.isinf(spline2_grid), 0, spline2_grid)
        #spline2_grid = np.where(np.isnan(spline2_grid), 0, spline2_grid)
        #spline2_grid = np.where(spline2_grid < 0, 0, spline2_grid)
        spline2_grid = np.where(Q2_mask, spline2_grid, -100)  # Mask unreliable regions
        print("spline1_grid",spline1_grid)
        print("spline2_grid",spline2_grid)

        # Calculate ratio
        ratio_grid = np.where(Q2_mask,10**(spline1_grid - spline2_grid),-100)
        print("ratio",ratio_grid)

        # Plot heatmap
        ax1.pcolormesh(X_mesh, Y_mesh, spline1_grid, cmap='viridis',
                       shading='gouraud',vmin=-50,vmax=-35)
        ax1.set_yscale('log')
        ax1.set_xscale('log')
        ax1.set_ylabel('Bjorken y')
        ax1.set_title(f'{spline1_label} Cross Section at E = {energy} GeV')
        cbar1 = plt.colorbar(ax1.collections[0], ax=ax1)
        cbar1.set_label(r'$d^2\sigma/dxdy$ [cm$^2$]')

        ax2.pcolormesh(X_mesh, Y_mesh, spline2_grid, cmap='viridis',
                       shading='gouraud',vmin=-50,vmax=-35)
        ax2.set_yscale('log')
        ax2.set_xscale('log')
        ax2.set_ylabel('Bjorken y')
        ax2.set_title(f'{spline2_label} Cross Section at E = {energy} GeV')
        cbar2 = plt.colorbar(ax2.collections[0], ax=ax2)
        cbar2.set_label(r'$d^2\sigma/dxdy$ [cm$^2$]')

        im = ax3.pcolormesh(X_mesh, Y_mesh, ratio_grid, cmap='RdBu_r',
                           vmin=0., vmax=2.0, shading='gouraud')

        cbar = plt.colorbar(im, ax=ax3)
        cbar.set_label(f'Ratio ({spline1_label}/{spline2_label})')

        ax3.set_xscale('log')
        ax3.set_yscale('log')
        ax3.set_xlabel('Bjorken x')
        ax3.set_ylabel('Bjorken y')
        ax3.set_title(f'Cross Section Ratio at E = {energy} GeV')

        plt.tight_layout()
        plt.savefig(f"{output_dir}/ratio_heatmap_E_{energy}.pdf", dpi=150, bbox_inches='tight')
        plt.close(fig)

# Main execution
if __name__ == "__main__":
    output_dir = "figures"

    # Define paths
    spline1_path = "/n/holylfs05/LABS/arguelles_delgado_lab/Everyone/nkamp/spack/var/spack/environments/lienv/.spack-env/._view/6xqcch3tgqvuwnqnpljuykkv3uhv5u4p/lib/python3.10/site-packages/siren/resources/CrossSections/CSMSDISSplines/CSMSDISSplines-v1.0/dsdxdy_nu_CC_iso.fits"
    spline2_path = "/n/holylfs05/LABS/arguelles_delgado_lab/Everyone/pweigel/cross_sections/20241112/wcg24b_dsdxdy_CC_muon_neutrino_isoscalar.fits"

    energies = np.logspace(1,5,5)
    x_values = np.logspace(-4,0,6)
    y_values = np.logspace(-4,0,100)

    # Check if files exist
    if not (os.path.exists(spline1_path) and os.path.exists(spline2_path)):
        print("Spline files not found. Please check the paths.")

    # Labels
    spline1_label = "Philip"
    spline2_label = "SIREN"

    compare_splines(spline1_path, spline2_path, spline1_label, spline2_label,
                    energies, x_values, y_values, output_dir)
