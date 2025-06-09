# Author: Amit Sud  
# Date: 1st May 2025  
# Description: This script generates KDE density maps of spatial transcriptomics data, restricted to selected polygon regions.  
# It visualizes the spatial distribution of specific cell types (e.g., NSUN5_MUT, NSUN5_specific_TCR) and exports high-resolution plots.  
# Input:  
#   - AnnData object with spatial coordinates in `.obs[['x', 'y']]`  
#   - `sample_id` indicating the sample to subset  
#   - `categories` list indicating which cell annotations to include in KDE  
#   - `polygon_coords`: boundary of the region of interest  
# Output:  
#   - KDE plot clipped to polygon with scale bar and density colorbar  
#   - Optional PDF output (1200 dpi)  

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
from shapely.geometry import Point, Polygon
from pathlib import Path

def generate_kde_with_polygon_clipping(adata, sample_id, categories, polygon_coords, grid_size=250, percentile=90):
    """Generates a KDE density plot with clipping outside a selected polygon."""
    
    # Extract spatial coordinates for sample
    all_data = adata[adata.obs['sample_id'] == sample_id]
    all_xy_coordinates = all_data.obs[['x', 'y']].values  

    # Subset to categories of interest
    kde_data = all_data[all_data.obs['NSUN5_MUT_TCR'].isin(categories)]
    kde_coordinates = kde_data.obs[['x', 'y']].values

    if kde_coordinates.shape[0] == 0:
        print("No data available for the selected categories.")
        return None, None, None, None, None

    # Kernel density estimate
    kde = gaussian_kde(kde_coordinates.T, bw_method='scott')

    # Grid boundaries
    x_min, x_max = all_xy_coordinates[:, 0].min(), all_xy_coordinates[:, 0].max()
    y_min, y_max = all_xy_coordinates[:, 1].min(), all_xy_coordinates[:, 1].max()
    
    x_grid, y_grid = np.meshgrid(
        np.linspace(x_min, x_max, grid_size),
        np.linspace(y_min, y_max, grid_size)
    )
    grid_coords = np.vstack([x_grid.ravel(), y_grid.ravel()])
    density = kde(grid_coords).reshape(x_grid.shape)

    # Polygon mask
    polygon = Polygon(polygon_coords)
    mask = np.zeros_like(density, dtype=bool)
    for i in range(grid_size):
        for j in range(grid_size):
            point = Point(x_grid[i, j], y_grid[i, j])
            if polygon.contains(point):
                mask[i, j] = True

    density[~mask] = np.nan  
    threshold = np.nanpercentile(density, percentile)

    return density, x_grid, y_grid, threshold, (x_min, x_max, y_min, y_max)


def plot_final_kde(adata, sample_id, categories, polygon_coords, save_pdf=False, pdf_filename="kde_plot.pdf"):
    """Plots KDE with polygon clipping, adds colorbar, scale bar, and exports to high-quality PDF if requested."""

    density, x_grid, y_grid, threshold, extent = generate_kde_with_polygon_clipping(
        adata, sample_id, categories, polygon_coords)

    if density is not None:
        fig, ax = plt.subplots(figsize=(8, 8))

        # Main KDE plot
        kde_plot = ax.imshow(density, origin='lower', extent=extent, cmap='coolwarm', alpha=0.8)

        # Colorbar
        cbar = fig.colorbar(kde_plot, ax=ax, shrink=0.75, aspect=20, pad=0.02)
        cbar.set_label('KDE Density', fontsize=12)
        cbar.ax.tick_params(labelsize=10)

        # Contour threshold
        ax.contour(x_grid, y_grid, density, levels=[threshold], colors='white', linewidths=1.5)

        # Hide axes
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_frame_on(False)
        ax.set_aspect('equal')

        # Polygon outline
        poly_x, poly_y = zip(*polygon_coords)
        ax.plot(poly_x + (poly_x[0],), poly_y + (poly_y[0],), 'white', linewidth=2)

        # Scale bar (500 µm)
        bottom_right_x = extent[0] + (extent[1] - extent[0]) * 0.90
        bottom_right_y = extent[2] + (extent[3] - extent[2]) * 0.08
        size_microns = 500
        size_pixels = size_microns  # 1 pixel = 1 micron

        ax.plot([bottom_right_x, bottom_right_x - size_pixels], [bottom_right_y, bottom_right_y], 'k-', linewidth=2)
        ax.text(bottom_right_x - size_pixels / 2, bottom_right_y - (extent[3] - extent[2]) * 0.015,
                f'{size_microns} µm', ha='center', va='top', fontsize=12, color='black')

        # Save if needed
        if save_pdf:
            output_path = Path("path/to/output") / pdf_filename
            plt.savefig(output_path, dpi=1200, bbox_inches='tight', format='pdf')
            print(f"✅ High-quality KDE plot saved as: {output_path}")

        plt.show()

# Example usage (with placeholders)
plot_final_kde(
    adata=adata_merged_filtered_seurat,
    sample_id='your_sample_id_here',
    categories=['NSUN5_MUT', 'NSUN5_specific_TCR'],
    polygon_coords=polygon_coords,  # Define this before running
    save_pdf=True,
    pdf_filename="kde_plot_NSUN5_MUT_TCR.pdf"
)
