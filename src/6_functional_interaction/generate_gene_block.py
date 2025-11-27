import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib as mpl
import numpy as np
from matplotlib.lines import Line2D

mpl.rcParams['font.family'] = 'NimbusSanL'

mpl.rcParams['font.family'] = 'NimbusSanL'
def plot_gene_block(gene_name: str, measurements: np.ndarray, vmin=-0.8, vmax=0.8,
                    colormap: str = 'BrBG_r', figsize=(2,1), incluster=True, output_filename: str = None):
    """
    Generates a 1x2 gene visualization block with a central gene name overlay.

    The color of each block (Left and Right) corresponds to the value in the 
    'measurements' array, mapped via the specified colormap.

    Args:
        gene_name (str): The name of the gene (e.g., 'MAP4K1').
        measurements (np.ndarray): A 1D NumPy array of two numerical values, 
                                   representing the measurement for each block.
                                   Structure: [Left_Block, Right_Block]
        colormap (str): The matplotlib colormap to use (e.g., 'Reds', 'viridis', etc.).
    """
    # 1. Setup the figure and axis
    # Adjusted figsize for a 1x2 horizontal look
    fig, ax = plt.subplots(figsize=(figsize[0], figsize[1]))
    
    # Reshape the 1D array (2 elements) into a 1x2 2D array for imshow.
    data_2d = measurements.reshape(1, 2)

    # 2. Plot the 1x2 data using imshow (creates the colored blocks)
    # extent=[xmin, xmax, ymin, ymax]. x runs from 0 to 2, y runs from 0 to 1.
    img = ax.imshow(data_2d, cmap=colormap, vmin=vmin, vmax=vmax, 
                    interpolation='nearest', extent=[0, figsize[0], 0, figsize[1]]) 

    # 3. Apply plot aesthetics and add black line borders/separators
    ax.set_xticks([]) # Remove x ticks
    ax.set_yticks([]) # Remove y ticks

    if incluster:
        line_style = '-'
        line_width = 3.5
        font_weight = 'bold'
        custom_dashes = [2,0]
    else:
        line_style = '--'
        line_width = 3.5
        font_weight = 'normal'
        custom_dashes = [2,2]
    
    # Set all spines (outer borders) to be visible and black
    for spine in ax.spines.values():
        spine.set_linewidth(line_width)
        spine.set_linestyle(line_style)
        spine.set_edgecolor('black')
        spine.set_visible(True)
        #spine.set_linestyle((0, custom_dashes))
    
    # Add a black vertical line to separate the two blocks (at x=1.0)
    
    line = Line2D([figsize[0]/2, figsize[0]/2], [0, figsize[1]], color='black', linewidth=line_width/2, zorder=1)
    line.set_linestyle((0, custom_dashes)) # Apply custom dashes
    ax.add_line(line)

    # Ensure the aspect ratio and limits are correct
    ax.set_xlim(0, figsize[0])
    ax.set_ylim(0, figsize[1]) # Y-limit is now 1 for the 1x2 shape

    # 4. Draw the central white box (Rectangle patch)
    # The center of the 1x2 grid is now at (1.0, 0.5)
    box_width = 1.7
    box_height = 0.5
    center_x = figsize[0]/2
    center_y = figsize[1]/2
    
    # Create the white rectangle
    rect = patches.Rectangle(
        (center_x - box_width/2, center_y - box_height/2),  # (x, y) starting point (bottom-left)
        box_width, 
        box_height,
        linewidth=line_width/2,
        linestyle=line_style,
        edgecolor='black', # Border color
        facecolor='white', # Fill color
        zorder=2 # Ensure the box is above the colored blocks
    )
    rect.set_linestyle((0, custom_dashes)) 
    ax.add_patch(rect)

    if len(gene_name)<8:
        fontsize=24
    else:
        fontsize=20
    # 5. Add the gene name text
    ax.text(
        center_x, center_y*0.95,
        gene_name,
        color='black',
        fontsize=fontsize,
        fontweight=font_weight,
        ha='center', # Horizontal alignment: Center
        va='center', # Vertical alignment: Center
        zorder=3 # Ensure text is above the white box
    )

    fig.subplots_adjust(left=0, right=1, top=1, bottom=0)
    
    # 6. Show the plot (or save it)
    if output_filename:
        # Save the figure with tight bounding box and high resolution (300 DPI)
        plt.savefig(output_filename,
                    bbox_inches='tight',
                    pad_inches=0,
                    dpi=300)
        plt.close(fig) # Close the figure to free memory
        print(f"Figure successfully saved to {output_filename}")
    else:
        plt.show()