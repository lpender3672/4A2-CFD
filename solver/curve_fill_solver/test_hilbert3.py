
import numpy as np
from collections import defaultdict

def find_hilbert_neighbors(cells):
    """
    Find all neighbors for each cell in an adaptive Hilbert grid.
    
    Args:
        cells: list of (py, px, side) tuples in Hilbert order
    
    Returns:
        neighbors: list of tuples, where neighbors[i] contains indices of 
                   all cells neighboring cells[i]
    """
    n_cells = len(cells)
    neighbors = [[] for _ in range(n_cells)]
    
    # Build spatial index: map (discretized_position, size) -> cell indices
    # We discretize to a fine grid to handle different cell sizes
    min_side = min(cell[2] for cell in cells)
    grid_resolution = min_side / 2  # Use half of smallest cell for discretization
    
    spatial_index = defaultdict(list)
    
    # Index all cells by their coverage area
    for idx, (py, px, side) in enumerate(cells):
        # Get all grid points this cell covers
        x_start = int(np.floor(px / grid_resolution))
        x_end = int(np.ceil((px + side) / grid_resolution))
        y_start = int(np.floor(py / grid_resolution))
        y_end = int(np.ceil((py + side) / grid_resolution))
        
        for gx in range(x_start, x_end):
            for gy in range(y_start, y_end):
                spatial_index[(gx, gy)].append(idx)
    
    # For each cell, find neighbors by checking adjacent grid regions
    for idx, (py, px, side) in enumerate(cells):
        x_start = int(np.floor(px / grid_resolution))
        x_end = int(np.ceil((px + side) / grid_resolution))
        y_start = int(np.floor(py / grid_resolution))
        y_end = int(np.ceil((py + side) / grid_resolution))
        
        candidate_neighbors = set()
        
        # Check grid cells around the perimeter (one layer outside)
        for gx in range(x_start - 1, x_end + 1):
            for gy in range(y_start - 1, y_end + 1):
                # Skip interior cells (only check perimeter)
                if x_start <= gx < x_end and y_start <= gy < y_end:
                    continue
                candidate_neighbors.update(spatial_index.get((gx, gy), []))
        
        # Verify actual adjacency for each candidate
        for neighbor_idx in candidate_neighbors:
            if neighbor_idx == idx:
                continue
                
            py2, px2, side2 = cells[neighbor_idx]
            
            # Check if cells are actually adjacent (share an edge or corner)
            if cells_are_adjacent(px, py, side, px2, py2, side2):
                neighbors[idx].append(neighbor_idx)
    
    # Convert to tuples for immutability
    return [tuple(sorted(set(n))) for n in neighbors]


def cells_are_adjacent(px1, py1, side1, px2, py2, side2, epsilon=1e-9):
    """
    Check if two rectangular cells are adjacent (share edge or corner).
    
    Cells are defined by bottom-left corner (px, py) and side length.
    """
    # Cell 1 bounds
    x1_min, x1_max = px1, px1 + side1
    y1_min, y1_max = py1, py1 + side1
    
    # Cell 2 bounds
    x2_min, x2_max = px2, px2 + side2
    y2_min, y2_max = py2, py2 + side2
    
    # Check if cells overlap or touch
    x_overlap = max(0, min(x1_max, x2_max) - max(x1_min, x2_min))
    y_overlap = max(0, min(y1_max, y2_max) - max(y1_min, y2_min))
    
    # Adjacent if they touch along edge (one dimension has zero gap, other has overlap)
    # or touch at corner (both dimensions have zero gap)
    x_touching = abs(x1_max - x2_min) < epsilon or abs(x2_max - x1_min) < epsilon
    y_touching = abs(y1_max - y2_min) < epsilon or abs(y2_max - y1_min) < epsilon
    
    x_overlapping = x_overlap > epsilon
    y_overlapping = y_overlap > epsilon
    
    # Adjacent if:
    # 1. Share an edge: touching in one direction, overlapping in other
    # 2. Share a corner: touching in both directions
    return (x_touching and y_overlapping) or (y_touching and x_overlapping) or (x_touching and y_touching)


def find_hilbert_neighbors_edge_only(cells, eps=1e-9):
    """
    Find edge-adjacent neighbors between variable-sized Hilbert cells.
    Works across mixed levels of detail.
    """
    n_cells = len(cells)
    neighbors = [[] for _ in range(n_cells)]

    # Base grid spacing: fine enough to resolve smallest cell edges
    min_side = min(c[2] for c in cells)
    grid_res = min_side / 2.0

    # Spatial hash: key = (gx, gy)
    spatial = defaultdict(list)

    for i, (py, px, s) in enumerate(cells):
        x0, x1 = px - s/2, px + s/2
        y0, y1 = py - s/2, py + s/2
        gx0, gx1 = int(np.floor(x0 / grid_res)), int(np.floor(x1 / grid_res))
        gy0, gy1 = int(np.floor(y0 / grid_res)), int(np.floor(y1 / grid_res))
        for gx in range(gx0, gx1 + 1):
            for gy in range(gy0, gy1 + 1):
                spatial[(gx, gy)].append(i)

    for i, (py, px, s) in enumerate(cells):
        x0, x1 = px - s/2, px + s/2
        y0, y1 = py - s/2, py + s/2
        gx0, gx1 = int(np.floor(x0 / grid_res)), int(np.floor(x1 / grid_res))
        gy0, gy1 = int(np.floor(y0 / grid_res)), int(np.floor(y1 / grid_res))

        # collect candidates in adjacent buckets
        cand = set()
        for gx in range(gx0 - 1, gx1 + 2):
            for gy in range(gy0 - 1, gy1 + 2):
                cand.update(spatial.get((gx, gy), []))

        for j in cand:
            if j == i:
                continue
            py2, px2, s2 = cells[j]
            if cells_share_edge(px, py, s, px2, py2, s2, eps):
                neighbors[i].append(j)

    # deduplicate and sort for consistency
    return neighbors #[tuple(sorted(set(n))) for n in neighbors]


def cells_share_edge(px1, py1, s1, px2, py2, s2, eps=1e-9):
    # Get bounds
    x1a, x1b = px1 - s1/2, px1 + s1/2
    y1a, y1b = py1 - s1/2, py1 + s1/2
    x2a, x2b = px2 - s2/2, px2 + s2/2
    y2a, y2b = py2 - s2/2, py2 + s2/2

    # Horizontal edge touch
    horiz_touch = (abs(y1a - y2b) < eps or abs(y1b - y2a) < eps)
    horiz_overlap = (x1b > x2a + eps) and (x2b > x1a + eps)

    # Vertical edge touch
    vert_touch = (abs(x1a - x2b) < eps or abs(x1b - x2a) < eps)
    vert_overlap = (y1b > y2a + eps) and (y2b > y1a + eps)

    return (horiz_touch and horiz_overlap) or (vert_touch and vert_overlap)


import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection


class InteractiveNeighborViewer:
    def __init__(self, cells, neighbors, poly, n, m, rotated=False):
        """
        Interactive viewer for Hilbert grid neighbors.
        
        Args:
            cells: list of (py, px, side) tuples
            neighbors: list of tuples containing neighbor indices for each cell
            n, m: domain (height, width) in same coordinate system as cells
            rotated: if True, swap x<->y for display
        """
        self.cells = cells
        self.neighbors = neighbors
        self.n = n
        self.m = m
        self.rotated = not rotated
        self.selected_cell = None
        self.poly = poly
        
        # Setup figure
        self.fig, self.ax = plt.subplots(figsize=(12, 10))
        self.fig.canvas.mpl_connect('button_press_event', self.on_click)
        
        # Draw initial state
        self.draw_grid()
        
        plt.tight_layout()
        plt.show()
    
    def draw_grid(self):
        """Draw all cells in the grid."""
        # Save current view limits
        xlim = self.ax.get_xlim()
        ylim = self.ax.get_ylim()
        
        self.ax.clear()
        
        # Determine which cells to highlight
        if self.selected_cell is not None:
            neighbor_indices = set(self.neighbors[self.selected_cell])
        else:
            neighbor_indices = set()
        
        # Prepare patches for all cells
        patches_default = []
        patches_neighbor = []
        patches_selected = []
        
        for idx, (py, px, side) in enumerate(self.cells):
            # Cell center is (px, py), half-extent is side/2
            x0 = px - 0.5 * side
            y0 = py - 0.5 * side
            w = side
            h = side
            
            if self.rotated:
                rect = Rectangle((y0, x0), h, w)
            else:
                rect = Rectangle((x0, y0), w, h)
            
            if idx == self.selected_cell:
                patches_selected.append(rect)
            elif idx in neighbor_indices:
                patches_neighbor.append(rect)
            else:
                patches_default.append(rect)
        
        # Draw default cells (gray)
        if patches_default:
            pcoll_default = PatchCollection(patches_default, 
                                           facecolor='lightgray',
                                           edgecolor='gray',
                                           alpha=0.3,
                                           linewidth=0.5)
            self.ax.add_collection(pcoll_default)
        
        # Draw neighbor cells (red)
        if patches_neighbor:
            pcoll_neighbor = PatchCollection(patches_neighbor,
                                            facecolor='red',
                                            edgecolor='darkred',
                                            alpha=0.5,
                                            linewidth=2)
            self.ax.add_collection(pcoll_neighbor)
        
        # Draw selected cell (blue)
        if patches_selected:
            pcoll_selected = PatchCollection(patches_selected,
                                            facecolor='blue',
                                            edgecolor='darkblue',
                                            alpha=0.6,
                                            linewidth=2.5)
            self.ax.add_collection(pcoll_selected)
        
        # Add cell index labels (scaled by cell size)
        for idx, (py, px, side) in enumerate(self.cells):
            text_color = 'white' if idx == self.selected_cell or idx in neighbor_indices else 'black'
            # Scale font size based on cell size (you can adjust the multiplier)
            fontsize = max(6, min(20, side * 3))  # clamp between 6 and 20
            
            if self.rotated:
                self.ax.text(py, px, str(idx),
                            ha='center', va='center',
                            fontsize=fontsize, fontweight='bold',
                            color=text_color)
            else:
                self.ax.text(px, py, str(idx),
                            ha='center', va='center',
                            fontsize=fontsize, fontweight='bold',
                            color=text_color)
        
        # Set axis properties matching your convention
        # Restore view limits if they existed, otherwise set defaults
        if xlim != (0.0, 1.0):  # matplotlib default
            self.ax.set_xlim(xlim)
            self.ax.set_ylim(ylim)
        else:
            if self.rotated:
                self.ax.set_xlim(0, self.n)
                self.ax.set_ylim(0, self.m)
            else:
                self.ax.set_xlim(0, self.m)
                self.ax.set_ylim(0, self.n)

        self.ax.plot(self.poly[:,0], self.poly[:,1], 'k-', linewidth=2)
        
        self.ax.set_aspect('equal')
        self.ax.grid(True, alpha=0.3)
        self.ax.set_xlabel('X', fontsize=10)
        self.ax.set_ylabel('Y', fontsize=10)
        
        # Add title with selection info
        if self.selected_cell is not None:
            nbr_list = list(self.neighbors[self.selected_cell])
            title = (f'Selected Cell: {self.selected_cell} | '
                    f'Neighbors: {nbr_list} ({len(nbr_list)} total)\n'
                    f'Click another cell or same cell to deselect')
        else:
            title = 'Click on a cell to highlight its neighbors\nSelected: Blue | Neighbors: Red | Others: Gray'
        
        self.ax.set_title(title, fontsize=11, pad=15)
        
        self.fig.canvas.draw()
    
    def on_click(self, event):
        """Handle mouse click events."""
        if event.inaxes != self.ax:
            return
        
        click_x, click_y = event.xdata, event.ydata
        
        # Find which cell was clicked
        clicked_cell = None
        for idx, (py, px, side) in enumerate(self.cells):
            # Cell boundaries
            x0 = px - 0.5 * side
            y0 = py - 0.5 * side
            x1 = px + 0.5 * side
            y1 = py + 0.5 * side
            
            # Check click position (accounting for rotation)
            if self.rotated:
                if y0 <= click_x <= y1 and x0 <= click_y <= x1:
                    clicked_cell = idx
                    break
            else:
                if x0 <= click_x <= x1 and y0 <= click_y <= y1:
                    clicked_cell = idx
                    break
        
        # Toggle selection
        if clicked_cell == self.selected_cell:
            # Clicking the same cell deselects it
            self.selected_cell = None
        else:
            # Select the new cell
            self.selected_cell = clicked_cell
        
        # Redraw
        self.draw_grid()


def visualize_neighbors(cells, neighbors, poly, m, n, rotated):
    """
    Launch interactive visualization of Hilbert grid neighbors.
    
    Args:
        cells: list of (py, px, side) tuples
        neighbors: list of tuples containing neighbor indices for each cell
    """
    viewer = InteractiveNeighborViewer(cells, neighbors, poly, m, n, rotated)
    return viewer


# Example usage
if __name__ == "__main__":

    from test_hilbert2 import transform_airfoil, naca4_airfoil, generate_hilbert_airfoil_lod

    m, n = 128, 64  # Grid dimensions
    alpha_deg = 5.0
    chord = 0.45 * m
    origin = (m*0.15, n*0.55)  # leading edge position
    alpha_deg = 6.0

    rotated = False
    if m > n:
        rotated = True
        n, m = m, n
        alpha_deg += 90
        origin = (origin[1], origin[0])

    poly01 = naca4_airfoil('2412', n=201, closed_te=True)
    foil = transform_airfoil(poly01, chord=chord, alpha_deg=5, origin=origin)
    if not np.allclose(foil[0], foil[-1]):
        foil = np.vstack([foil, foil[0]])

    # Generate Hilbert cells with solid fraction
    cells, fine_bits = generate_hilbert_airfoil_lod(m, n, foil, extra_global_levels=1)
    
    neighbors = find_hilbert_neighbors(cells)
    

    edge_neighbors = find_hilbert_neighbors_edge_only(cells)


    visualize_neighbors(cells, edge_neighbors, foil, m, n, rotated)

