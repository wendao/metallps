import numpy as np
from Bio.PDB import *
import numpy.linalg as la

def read_protein_coordinates(pdb_file):
    """Read main chain atom coordinates from PDB file."""
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('protein', pdb_file)
    
    main_chain_atoms = ['N', 'CA', 'C', 'O']
    coords = []
    
    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    if atom.name in main_chain_atoms:
                        coords.append(atom.coord)
    
    return np.array(coords)

def perform_pca(coords):
    """Perform PCA on coordinates."""
    # Calculate centroid
    centroid = np.mean(coords, axis=0)
    
    # Center the coordinates
    centered_coords = coords - centroid
    
    # Compute covariance matrix
    cov_matrix = np.cov(centered_coords.T)
    
    # Compute eigenvalues and eigenvectors
    eigenvalues, eigenvectors = la.eig(cov_matrix)
    
    # Sort eigenvectors by eigenvalues in descending order
    sorted_indices = np.argsort(eigenvalues)[::-1]
    principal_components = eigenvectors[:, sorted_indices]
    
    return centroid, principal_components, eigenvalues[sorted_indices]

def generate_points_along_pc(centroid, principal_component, num_points=100, step=0.4):
    """Generate points along principal component direction."""
    # Normalize the principal component
    pc_normalized = principal_component / np.linalg.norm(principal_component)
    
    # Generate points
    points = np.array([
        centroid + i * step * pc_normalized 
        for i in range(-(num_points//2), num_points//2 + 1)
    ])
    
    return points

def main(pdb_file):
    """Perform PCA analysis on protein structure."""
    # Read coordinates
    coords = read_protein_coordinates(pdb_file)
    
    # Perform PCA
    centroid, pcs, eigenvalues = perform_pca(coords)
    
    print("Centroid:", centroid)
    print("\nPrincipal Components:\n", pcs)
    print("\nEigenvalues:", eigenvalues)
    
    # Generate points along first principal component
    points = generate_points_along_pc(centroid, pcs[:, 0])
    for i, point in enumerate(points):
        print(f"HETATM {i+1:4d}  ZN   ZN X{i+1:4d}    {point[0]:8.4f}{point[1]:8.4f}{point[2]:8.4f}  1.00  0.00          Zn")

import sys
# Example usage
if __name__ == "__main__":
    main(sys.argv[1])
