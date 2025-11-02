from typing import List, Tuple, Dict, Optional, Union
import numpy as np
import numpy.typing as npt
from scipy.linalg import expm
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

class MacachorState:
    """
    Advanced Macachor Scalar Framework State Implementation
    Implements scalar-quaternion quantum states with enhanced features
    """
    
    def __init__(self, 
                 phi_0: complex, 
                 phi_vec: List[complex], 
                 orientation: str = "default") -> None:
        # Enhanced initialization with validation
        self.scalar: complex = complex(phi_0)
        self.vector: npt.NDArray[np.complex_] = np.array(phi_vec, dtype=complex)
        self.orientation: str = orientation
        self.history: List[Dict] = []  # Track state evolution
        self.time_points: List[float] = []  # Track evolution time
        self._store_state("Initialized", 0.0)
        
        # Validate dimensions
        if len(self.vector) != 3:
            raise ValueError("Vector component must have exactly 3 elements (i,j,k basis)")
    
    def _store_state(self, note: str, time: float) -> None:
        """Store state history for analysis"""
        self.history.append({
            'scalar': self.scalar,
            'vector': self.vector.copy(),
            'orientation': self.orientation,
            'note': note,
            'time': time
        })
        self.time_points.append(time)
    
    def __repr__(self) -> str:
        return (f"MacachorState(\n"
                f"  scalar={self.scalar:.3f},\n"
                f"  vector={[f'{x:.3f}' for x in self.vector]},\n"
                f"  orientation='{self.orientation}',\n"
                f"  norm={self.norm():.4f}\n"
                f")")

    def norm(self) -> float:
        """Calculate the norm of the state"""
        return np.sqrt(np.abs(self.scalar)**2 + np.sum(np.abs(self.vector)**2))
    
    def normalize(self) -> None:
        """Normalize the state to unit norm"""
        current_norm = self.norm()
        if current_norm > 0:
            self.scalar /= current_norm
            self.vector /= current_norm
            self._store_state("Normalized", self.time_points[-1] if self.time_points else 0.0)
    
    def scalar_vector_entanglement(self) -> float:
        """
        Calculate scalar-vector entanglement measure
        Returns entanglement entropy between scalar and vector components
        """
        # Construct density matrix for scalar-vector system
        psi_total = np.array([self.scalar] + self.vector.tolist())
        rho_total = np.outer(psi_total, np.conj(psi_total))
        
        # Partial trace over vector components to get scalar reduced density matrix
        rho_scalar = np.array([[np.abs(self.scalar)**2]])  # 1x1 density matrix
        
        # Calculate entanglement entropy
        if rho_scalar[0,0] > 0 and rho_scalar[0,0] < 1:
            entropy = -rho_scalar[0,0] * np.log(rho_scalar[0,0]) - (1-rho_scalar[0,0]) * np.log(1-rho_scalar[0,0])
        else:
            entropy = 0.0
            
        return entropy
    
    def macachor_energy(self, 
                       H_kinetic: npt.NDArray[np.complex_], 
                       H_scalar: npt.NDArray[np.complex_], 
                       H_vector: npt.NDArray[np.complex_], 
                       H_mix: npt.NDArray[np.complex_]) -> float:
        """
        Calculate expectation value of Macachor Hamiltonian
        """
        # Construct full state vector
        psi = np.array([self.scalar] + self.vector.tolist(), dtype=complex)
        psi_dagger = np.conj(psi)
        
        # Calculate total Hamiltonian expectation
        H_total = H_kinetic + H_scalar + H_vector + H_mix
        energy = np.dot(psi_dagger, np.dot(H_total, psi)).real
        
        return energy

    def evolve(self, 
               Hamiltonian: npt.NDArray[np.complex_], 
               dt: float, 
               method: str = 'exact') -> None:
        """
        Enhanced evolution with multiple methods
        """
        h_bar = 1.0  # Natural units
        
        # Construct state vector
        psi = np.array([self.scalar] + self.vector.tolist(), dtype=complex)
        
        if method == 'exact':
            # Exact matrix exponential
            Evolution_Operator = expm(-1j * Hamiltonian * dt / h_bar)
        elif method == 'taylor':
            # Taylor approximation (original method)
            Identity = np.identity(len(psi), dtype=complex)
            Evolution_Operator = Identity - 1j * Hamiltonian * dt / h_bar
        else:
            raise ValueError("Method must be 'exact' or 'taylor'")
        
        # Apply evolution
        psi_new = np.dot(Evolution_Operator, psi)
        
        # Update state
        self.scalar = psi_new[0]
        self.vector = psi_new[1:]
        current_time = self.time_points[-1] + dt if self.time_points else dt
        self._store_state(f"Evolved dt={dt}, method={method}", current_time)

    def measure(self, 
                operator: npt.NDArray[np.complex_], 
                collapse: bool = False) -> float:
        """
        Enhanced measurement with optional state collapse
        """
        psi = np.array([self.scalar] + self.vector.tolist(), dtype=complex)
        psi_dagger = np.conj(psi)
        
        # Calculate expectation value
        op_psi = np.dot(operator, psi)
        expectation_value = np.dot(psi_dagger, op_psi).real
        
        if collapse:
            # Simulate measurement collapse
            eigenvalues, eigenvectors = np.linalg.eigh(operator)
            probabilities = np.abs(np.dot(eigenvectors.T, psi))**2
            probabilities /= np.sum(probabilities)  # Normalize
            
            # Collapse to eigenstate based on probabilities
            chosen_idx = np.random.choice(len(eigenvalues), p=probabilities)
            collapsed_state = eigenvectors[:, chosen_idx]
            
            self.scalar = collapsed_state[0]
            self.vector = collapsed_state[1:]
            current_time = self.time_points[-1] if self.time_points else 0.0
            self._store_state(f"Collapsed to eigenvalue {eigenvalues[chosen_idx]:.3f}", current_time)
        
        return expectation_value

    def quaternion_phase(self) -> npt.NDArray[np.complex_]:
        """
        Calculate the quaternion phase factors from vector components
        """
        Q_norm = np.sqrt(np.sum(np.abs(self.vector)**2))
        if Q_norm > 0:
            phase_factors = self.vector / Q_norm
            return phase_factors
        else:
            return np.array([0+0j, 0+0j, 0+0j])

    def macachor_uncertainty(self, 
                            operator_A: npt.NDArray[np.complex_], 
                            operator_B: npt.NDArray[np.complex_]) -> Tuple[float, float]:
        """
        Calculate Macachor-enhanced uncertainty relation
        """
        # Expectation values
        exp_A = self.measure(operator_A)
        exp_B = self.measure(operator_B)
        
        # Commutator expectation
        commutator = np.dot(operator_A, operator_B) - np.dot(operator_B, operator_A)
        exp_commutator = self.measure(commutator)
        
        # Enhanced uncertainty with scalar-vector correction
        uncertainty = np.sqrt(
            (self.measure(np.dot(operator_A, operator_A)) - exp_A**2) *
            (self.measure(np.dot(operator_B, operator_B)) - exp_B**2)
        )
        
        macachor_enhancement = np.abs(exp_commutator) / 2 + 0.1 * self.scalar_vector_entanglement()
        
        return uncertainty, macachor_enhancement

    def plot_evolution(self, 
                      show_components: bool = True,
                      show_entanglement: bool = True,
                      show_quaternion_space: bool = True) -> plt.Figure:
        """
        Comprehensive visualization of state evolution history
        """
        if len(self.history) < 2:
            print("Insufficient history for plotting")
            return None
        
        # Determine subplot layout
        n_plots = sum([show_components, show_entanglement, show_quaternion_space])
        if n_plots == 0:
            return None
            
        fig, axes = plt.subplots(1, n_plots, figsize=(5*n_plots, 4))
        if n_plots == 1:
            axes = [axes]
        
        plot_idx = 0
        times = self.time_points
        
        # Plot 1: Component evolution
        if show_components:
            scalar_mags = [np.abs(state['scalar']) for state in self.history]
            vector_mags = [np.sqrt(sum(np.abs(v)**2 for v in state['vector'])) for state in self.history]
            
            axes[plot_idx].plot(times, scalar_mags, 'b-', label='|φ₀| (Scalar)', linewidth=2)
            axes[plot_idx].plot(times, vector_mags, 'r--', label='|Q| (Vector)', linewidth=2)
            axes[plot_idx].set_xlabel('Time')
            axes[plot_idx].set_ylabel('Magnitude')
            axes[plot_idx].set_title('Component Evolution')
            axes[plot_idx].legend()
            axes[plot_idx].grid(True, alpha=0.3)
            plot_idx += 1
        
        # Plot 2: Entanglement evolution
        if show_entanglement:
            # Calculate entanglement for each historical state
            entanglement_vals = []
            for state in self.history:
                temp_state = MacachorState(state['scalar'], state['vector'], "temp")
                entanglement_vals.append(temp_state.scalar_vector_entanglement())
            
            axes[plot_idx].plot(times, entanglement_vals, 'g-', linewidth=2)
            axes[plot_idx].set_xlabel('Time')
            axes[plot_idx].set_ylabel('Entanglement Entropy')
            axes[plot_idx].set_title('Scalar-Vector Entanglement')
            axes[plot_idx].grid(True, alpha=0.3)
            plot_idx += 1
        
        # Plot 3: Quaternion phase space (3D)
        if show_quaternion_space and len(self.history) > 1:
            # Extract vector components for 3D plot
            vec_x = [state['vector'][0].real for state in self.history]
            vec_y = [state['vector'][1].real for state in self.history] 
            vec_z = [state['vector'][2].real for state in self.history]
            
            # Create 3D plot
            ax_3d = fig.add_subplot(1, n_plots, plot_idx+1, projection='3d')
            scatter = ax_3d.scatter(vec_x, vec_y, vec_z, c=times, cmap='viridis', s=50)
            ax_3d.plot(vec_x, vec_y, vec_z, 'k-', alpha=0.3)
            ax_3d.set_xlabel('Q₁ (Re)')
            ax_3d.set_ylabel('Q₂ (Re)') 
            ax_3d.set_zlabel('Q₃ (Re)')
            ax_3d.set_title('Quaternion Vector Evolution')
            plt.colorbar(scatter, ax=ax_3d, label='Time')
        
        plt.tight_layout()
        return fig

    def get_state_statistics(self) -> Dict[str, float]:
        """
        Compute comprehensive statistics about state evolution
        """
        if len(self.history) < 2:
            return {}
        
        scalar_mags = [np.abs(state['scalar']) for state in self.history]
        vector_mags = [np.sqrt(sum(np.abs(v)**2 for v in state['vector'])) for state in self.history]
        
        # Calculate entanglement history
        entanglement_vals = []
        for state in self.history:
            temp_state = MacachorState(state['scalar'], state['vector'], "temp")
            entanglement_vals.append(temp_state.scalar_vector_entanglement())
        
        return {
            'scalar_mean': float(np.mean(scalar_mags)),
            'scalar_std': float(np.std(scalar_mags)),
            'vector_mean': float(np.mean(vector_mags)),
            'vector_std': float(np.std(vector_mags)),
            'entanglement_mean': float(np.mean(entanglement_vals)),
            'entanglement_max': float(np.max(entanglement_vals)),
            'correlation': float(np.corrcoef(scalar_mags, vector_mags)[0,1]),
            'total_evolution_time': self.time_points[-1] - self.time_points[0]
        }

# Enhanced Hamiltonian creation with type hints
def create_macachor_hamiltonian(mass: float = 1.0, 
                               scalar_potential: float = 0.1, 
                               vector_coupling: float = 0.05, 
                               mixing: float = 0.02) -> npt.NDArray[np.complex_]:
    """Create a sample Macachor Hamiltonian with type hints"""
    H_kinetic = np.array([
        [1.0, 0.0, 0.0, 0.0],
        [0.0, 0.5, 0.0, 0.0],
        [0.0, 0.0, 0.5, 0.0],
        [0.0, 0.0, 0.0, 0.5]
    ], dtype=complex) / (2 * mass)

    H_scalar = scalar_potential * np.array([
        [1.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.0]
    ], dtype=complex)

    H_vector = vector_coupling * np.array([
        [0.0, 0.0, 0.0, 0.0],
        [0.0, 1.0, 0.2, 0.2],
        [0.0, 0.2, 1.0, 0.2],
        [0.0, 0.2, 0.2, 1.0]
    ], dtype=complex)

    H_mix = mixing * np.array([
        [0.0, 0.1, 0.1, 0.1],
        [0.1, 0.0, 0.0, 0.0],
        [0.1, 0.0, 0.0, 0.0],
        [0.1, 0.0, 0.0, 0.0]
    ], dtype=complex)

    return H_kinetic + H_scalar + H_vector + H_mix

# Demonstration with enhanced visualization
if __name__ == "__main__":
    print("=== Enhanced Macachor Scalar Framework with Visualization ===")
    
    # Create initial state
    print("\n1. Creating initial Macachor state:")
    state = MacachorState(
        phi_0=0.8,
        phi_vec=[0.2j, 0.3, 0.1-0.1j],
        orientation="quantum_oscillator"
    )
    print(state)
    
    # Normalize state
    print("\n2. Normalizing state:")
    state.normalize()
    print(f"Normalized norm: {state.norm():.6f}")
    
    # Create Hamiltonian
    H = create_macachor_hamiltonian()
    
    # Evolve state through multiple time steps
    print("\n3. Evolving state through multiple time steps...")
    time_steps = 20
    for i in range(time_steps):
        state.evolve(H, dt=0.1, method='exact')
    
    # Generate comprehensive visualization
    print("\n4. Generating evolution plots...")
    fig = state.plot_evolution(show_components=True, 
                              show_entanglement=True, 
                              show_quaternion_space=True)
    plt.savefig('macachor_evolution.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    # Display statistics
    print("\n5. Evolution Statistics:")
    stats = state.get_state_statistics()
    for key, value in stats.items():
        print(f"   {key}: {value:.4f}")
    
    # Final measurements
    position_op = np.array([
        [1.0, 0.1, 0.1, 0.1],
        [0.1, 0.8, 0.0, 0.0],
        [0.1, 0.0, 0.8, 0.0],
        [0.1, 0.0, 0.0, 0.8]
    ], dtype=complex)
    
    momentum_op = 1j * np.array([
        [0.0, -0.5, -0.5, -0.5],
        [0.5, 0.0, -0.2, 0.2],
        [0.5, 0.2, 0.0, -0.2],
        [0.5, -0.2, 0.2, 0.0]
    ], dtype=complex)
    
    uncertainty, enhancement = state.macachor_uncertainty(position_op, momentum_op)
    print(f"\n6. Final Uncertainty Analysis:")
    print(f"   Standard uncertainty: {uncertainty:.4f}")
    print(f"   Macachor enhancement: {enhancement:.4f}")
    print(f"   Enhanced bound: ΔxΔp ≥ {enhancement:.4f}")
    
    print(f"\n7. Total evolution history: {len(state.history)} states")
    print("=== Enhanced Demonstration Complete ===")