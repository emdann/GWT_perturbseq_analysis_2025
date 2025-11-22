import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Circle
import numpy as np


def plot_gene_perturbation(
    perturbed_genes,
    positive_genes=None,
    negative_genes=None,
    highlight_genes=None,
    context_label="Context",
    context_color="#d9534f",
    positive_label="Positively affected",
    negative_label="Negatively affected",
    figsize=(12, 6),
    save_path=None
):
    """
    Create a visualization of gene perturbation effects.
    
    Parameters:
    -----------
    perturbed_genes : list
        List of perturbed/upstream genes
    positive_genes : list, optional
        List of positively affected downstream genes
    negative_genes : list, optional
        List of negatively affected downstream genes
    highlight_genes : list, optional
        List of perturbed genes to highlight (will be shown in bold/red)
    context_label : str
        Label for the context annotation
    context_color : str
        Color for the context circle/border (hex or named color)
    positive_label : str
        Label for the positive effects box
    negative_label : str
        Label for the negative effects box
    figsize : tuple
        Figure size (width, height)
    save_path : str, optional
        If provided, save the figure to this path
    
    Returns:
    --------
    fig, ax : matplotlib figure and axis objects
    """
    
    fig, ax = plt.subplots(figsize=figsize)
    
    # Box parameters
    box_width = 1.8
    line_height = 0.35  # Height per gene line
    
    # Calculate box heights based on number of genes
    perturbed_box_height = max(0.8, len(perturbed_genes) * line_height + 0.5)
    
    # Calculate downstream box heights
    has_positive = positive_genes and len(positive_genes) > 0
    has_negative = negative_genes and len(negative_genes) > 0
    
    if has_positive:
        positive_box_height = max(0.8, len(positive_genes) * line_height + 0.5)
    if has_negative:
        negative_box_height = max(0.8, len(negative_genes) * line_height + 0.5)
    
    # Calculate total height needed
    if has_positive and has_negative:
        total_downstream_height = positive_box_height + negative_box_height + 0.8
    elif has_positive:
        total_downstream_height = positive_box_height
    elif has_negative:
        total_downstream_height = negative_box_height
    else:
        total_downstream_height = 0
    
    # Determine the maximum height needed
    max_content_height = max(perturbed_box_height, total_downstream_height)
    
    # Set canvas limits with padding
    canvas_height = max_content_height + 3.5  # Extra space for labels and context
    canvas_width = 10
    
    ax.set_xlim(0, canvas_width)
    ax.set_ylim(0, canvas_height)
    ax.axis('off')
    
    # Center everything vertically
    center_y = canvas_height / 2
    
    # Calculate positions
    perturbed_x = 1.5
    perturbed_y = center_y - perturbed_box_height / 2
    
    downstream_x = 5.5
    
    # Position downstream boxes centered vertically
    if has_positive and has_negative:
        # Stack them vertically, centered
        spacing = 0.8
        total_height = positive_box_height + spacing + negative_box_height
        start_y = center_y + total_height / 2
        
        positive_y = start_y - positive_box_height
        negative_y = positive_y - spacing - negative_box_height
        
    elif has_positive:
        positive_y = center_y - positive_box_height / 2
        
    elif has_negative:
        negative_y = center_y - negative_box_height / 2
    
    # Draw perturbed genes box
    perturbed_box = FancyBboxPatch(
        (perturbed_x, perturbed_y), box_width, perturbed_box_height,
        boxstyle="round,pad=0.1", 
        edgecolor='black', 
        facecolor='white', 
        linewidth=2
    )
    ax.add_patch(perturbed_box)
    
    # Add perturbed genes text
    y_offset = perturbed_y + perturbed_box_height - 0.35
    for gene in perturbed_genes:
        if highlight_genes and gene in highlight_genes:
            ax.text(perturbed_x + box_width/2, y_offset, gene, 
                   ha='center', va='top', fontsize=14, 
                   weight='bold', color='red')
        else:
            ax.text(perturbed_x + box_width/2, y_offset, gene, 
                   ha='center', va='top', fontsize=14)
        y_offset -= line_height
    
    # Draw positive effects box and arrow
    if has_positive:
        positive_box = FancyBboxPatch(
            (downstream_x, positive_y), box_width, positive_box_height,
            boxstyle="round,pad=0.1", 
            edgecolor='black', 
            facecolor='white', 
            linewidth=2
        )
        ax.add_patch(positive_box)
        
        # Add positive label above box
        ax.text(downstream_x + box_width/2, positive_y + positive_box_height + 0.15, 
               positive_label, ha='center', va='bottom', 
               fontsize=13, style='italic', color='black')
        
        # Add positive genes text
        y_offset = positive_y + positive_box_height - 0.35
        for gene in positive_genes:
            ax.text(downstream_x + box_width/2, y_offset, gene, 
                   ha='center', va='top', fontsize=14)
            y_offset -= line_height
        
        # Draw activation arrow (-->)
        arrow_y = positive_y + positive_box_height/2
        arrow = FancyArrowPatch(
            (perturbed_x + box_width + 0.1, center_y), 
            (downstream_x - 0.1, arrow_y),
            arrowstyle='->', 
            mutation_scale=25, 
            linewidth=2.5, 
            color='black'
        )
        ax.add_patch(arrow)
    
    # Draw negative effects box and inhibition line
    if has_negative:
        negative_box = FancyBboxPatch(
            (downstream_x, negative_y), box_width, negative_box_height,
            boxstyle="round,pad=0.1", 
            edgecolor='black', 
            facecolor='white', 
            linewidth=2
        )
        ax.add_patch(negative_box)
        
        # Add negative label above box
        ax.text(downstream_x + box_width/2, negative_y + negative_box_height + 0.15, 
               negative_label, ha='center', va='bottom', 
               fontsize=13, style='italic', color='black')
        
        # Add negative genes text
        y_offset = negative_y + negative_box_height - 0.35
        for gene in negative_genes:
            ax.text(downstream_x + box_width/2, y_offset, gene, 
                   ha='center', va='top', fontsize=14)
            y_offset -= line_height
        
        # Draw inhibition line (---|)
        arrow_y = negative_y + negative_box_height/2
        
        # Line part
        ax.plot([perturbed_x + box_width + 0.1, downstream_x - 0.3], 
               [center_y, arrow_y], 
               color='black', linewidth=2.5)
        
        # Perpendicular bar at the end
        ax.plot([downstream_x - 0.3, downstream_x - 0.3], 
               [arrow_y - 0.15, arrow_y + 0.15], 
               color='black', linewidth=3)
    
    # Draw context circle/border around both boxes
    if has_positive and has_negative:
        # Calculate bounding box for both downstream boxes
        top_y = positive_y + positive_box_height
        bottom_y = negative_y
        context_height = top_y - bottom_y + 0.6
        context_y = bottom_y - 0.3
    elif has_positive:
        context_height = positive_box_height + 0.6
        context_y = positive_y - 0.3
    elif has_negative:
        context_height = negative_box_height + 0.6
        context_y = negative_y - 0.3
    else:
        # No downstream genes
        context_height = 0
        context_y = 0
    
    # Extend context to include perturbed genes box as well
    context_full_height = max(context_y + context_height, perturbed_y + perturbed_box_height + 0.3) - min(context_y, perturbed_y - 0.3) + 0.6
    context_full_y = min(context_y, perturbed_y - 0.3) - 0.3
    context_full_width = (downstream_x + box_width + 0.3) - (perturbed_x - 0.3)
    context_full_x = perturbed_x - 0.3
    
    # Draw thick rounded rectangle for context
    context_border = FancyBboxPatch(
        (context_full_x, context_full_y), 
        context_full_width, 
        context_full_height,
        boxstyle="round,pad=0.3", 
        edgecolor=context_color, 
        facecolor='none', 
        linewidth=4,
        linestyle='-'
    )
    ax.add_patch(context_border)
    
    # Add context label below the context box
    ax.text(context_full_x + context_full_width/2, context_full_y - 0.5, 
           f"Context: {context_label}", 
           ha='center', va='top', 
           fontsize=14, color=context_color, weight='bold')
    
    # Add "Perturbed genes" and "Downstream genes" labels at same height
    label_y = min(negative_y if has_negative else (positive_y if has_positive else perturbed_y), 
                  perturbed_y, context_full_y) - 0.7
    
    # ax.text(perturbed_x + box_width/2, label_y, 
    #        "Perturbed\ngenes", 
    #        ha='center', va='top', 
    #        fontsize=13, style='italic', color='black')
    
    # if has_positive or has_negative:
    #     downstream_center_x = downstream_x + box_width/2
    #     ax.text(downstream_center_x, label_y, 
    #            "Downstream\ngenes", 
    #            ha='center', va='top', 
    #            fontsize=13, style='italic', color='black')
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Figure saved to {save_path}")
    
    return fig, ax 