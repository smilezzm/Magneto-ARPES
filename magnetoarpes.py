# -*- coding: utf-8 -*-
"""
MagnetoARPES Simulation Class

A comprehensive class for simulating electron trajectories in magnetic fields
for Angle-Resolved Photoemission Spectroscopy (ARPES) experiments.

Author: Original Jounghoon's codes, modified by Mengyao and Claude Sonnet 4
Date: September 2025
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.special import ellipk, ellipe
from scipy.interpolate import RegularGridInterpolator, LinearNDInterpolator
import time
from typing import Tuple, Optional, List, Union
import warnings


class MagnetoARPES:
    """
    A class for simulating electron trajectories in magnetic fields for ARPES experiments.
    
    This class handles:
    - Magnetic field calculation for solenoid coils with rotation/translation
    - Boris pusher algorithm for charged particle trajectories
    - Forward mapping: initial velocities → final velocities
    - Reverse mapping: final velocities → initial velocities (energy-sliced)
    """
    
    def __init__(self, 
                 photon_energy: float = 115.0,  # eV
                 work_function: float = 4.29,   # eV
                 binding_energy: float = 0.0,   # eV
                 flight_length: float = 0.038,  # m
                 coil_inner_d: float = 6.3e-3,  # m
                 coil_outer_d: float = 12e-3,   # m
                 coil_height: float = 5.6e-3,   # m
                 wire_diameter: float = 0.165e-3,  # m
                 current: float = 0.2,          # A
                 coil_rotation: Tuple[float, float, float] = (0.0, 0.0, 0.0),  # radians
                 coil_translation: Tuple[float, float, float] = (0.0, 0.0, 0.0),  # m
                 grid_resolution: int = 21,  # grid points for B-field interpolation
                 max_sim_steps: int = 400):
        """
        Initialize the MagnetoARPES simulator.
        
        Parameters:
        -----------
        photon_energy : float
            Photon energy in eV
        work_function : float
            Work function in eV
        binding_energy : float
            Binding energy in eV
        flight_length : float
            Distance from sample to detector in m
        coil_inner_d : float
            Inner diameter of solenoid coil in m
        coil_outer_d : float
            Outer diameter of solenoid coil in m
        coil_height : float
            Height of solenoid coil in m
        wire_diameter : float
            Wire diameter in m
        current : float
            Current through coil in A
        coil_rotation : tuple
            Rotation angles (theta_x, theta_y, theta_z) in radians
        coil_translation : tuple
            Translation vector (x, y, z) in m
        grid_resolution : int
            Number of grid points for B-field interpolation
        max_sim_steps : int
            Maximum simulation steps per trajectory
        """
        
        # Physical constants
        self.e = 1.602e-19    # electron charge
        self.m = 9.109e-31    # electron mass
        self.hbar = 1.055e-34 # reduced Planck constant
        self.mu0 = 4*np.pi*1e-7  # permeability of free space
        
        # Experimental parameters
        self.hv = photon_energy
        self.phi = work_function
        self.E_B = binding_energy
        self.L_flight = flight_length
        
        # Calculate electron kinetic energy and velocity
        self.Ekin = (self.hv - self.phi - self.E_B) * self.e
        self.v0 = np.sqrt(2 * self.Ekin / self.m)
        self.kF = self.m * self.v0 / self.hbar
        
        # Coil parameters
        self.coil_inner_d = coil_inner_d
        self.coil_outer_d = coil_outer_d
        self.coil_height = coil_height
        self.wire_d = wire_diameter
        self.I = current
        
        # Coil positioning
        self.theta_x, self.theta_y, self.theta_z = coil_rotation
        self.trans = np.array(coil_translation)
        
        # Simulation parameters
        self.grid_resolution = grid_resolution
        self.max_sim_steps = max_sim_steps
        
        # Calculate coil geometry
        self._setup_coil_geometry()
        
        # Setup magnetic field interpolation
        self._setup_magnetic_field()
        
        # Setup time stepping
        self._setup_time_stepping()
        
        # Initialize reverse mapping storage
        self.inv_map_cell = None
        self.v0_range = None
        
        print(f"MagnetoARPES initialized:")
        print(f"  Kinetic energy: {self.Ekin/self.e:.2f} eV")
        print(f"  Electron speed: {self.v0:.2e} m/s")
        print(f"  Coil turns (total number of loops): {self.total_turns}")
        print(f"  B-field at center: {self.B_mag_center*1e3:.3f} mT")
    
    def _setup_coil_geometry(self):
        """Setup coil geometry parameters"""
        self.R_inner = self.coil_inner_d / 2
        self.R_outer = self.coil_outer_d / 2
        self.z_start = -self.coil_height
        self.z_end = 0.0
        
        self.n_r = int(np.floor((self.R_outer - self.R_inner) / self.wire_d))
        self.n_z = int(np.floor(self.coil_height / self.wire_d))
        self.total_turns = self.n_r * self.n_z
        
        self.rr = np.linspace(self.R_inner + self.wire_d/2, 
                              self.R_outer - self.wire_d/2, self.n_r)
        self.zz = np.linspace(self.z_start + self.wire_d/2, 
                              self.z_end - self.wire_d/2, self.n_z)
    
    def _setup_magnetic_field(self):
        """Setup magnetic field interpolation grid"""
        # Create interpolation grid
        grid_extent = 0.05  # 5 cm
        x_vals = np.linspace(-grid_extent, grid_extent, self.grid_resolution)
        y_vals = np.linspace(-grid_extent, grid_extent, self.grid_resolution)
        z_vals = np.linspace(0.0, self.L_flight, self.grid_resolution)
        
        X, Y, Z = np.meshgrid(x_vals, y_vals, z_vals, indexing='ij')  # X, Y, Z are all (grid_resolution, grid_resolution, grid_resolution)
        X = X.flatten() # (grid_resolution**3,)
        Y = Y.flatten()
        Z = Z.flatten()
        print("Calculating 3D B-field grid (elliptic integral)...")
        B_grid = self._bfield_solenoid_rot_trans(X, Y, Z)  # (grid_resolution, grid_resolution, grid_resolution, 3)
        
        # Check for and handle singularities
        self._fix_singularities(B_grid)
        
        # Create interpolator
        self.B_interp = RegularGridInterpolator((x_vals, y_vals, z_vals), B_grid,
                                                bounds_error=False, fill_value=0.0)
        
        # Calculate reference magnetic fields
        self.B_center = self.B_interp([[0.0, 0.0, 0.0]])[0]
        self.B_mag_center = np.linalg.norm(self.B_center)
        
        self.B_center_analyzer = self.B_interp([[0.0, 0.0, self.L_flight]])[0]
        self.B_mag_center_analyzer = np.linalg.norm(self.B_center_analyzer)
    
    def _setup_time_stepping(self):
        """Setup time stepping parameters for Boris pusher"""
        omega_c = self.e * self.B_mag_center / self.m
        T_c = 2 * np.pi / omega_c
        steps_per_T = 100
        self.dt = T_c / steps_per_T
        
        n_steps_estimate = int(np.ceil(self.L_flight / self.v0 / self.dt))
        print(f"Time step: {self.dt:.3e} s")
        print(f"Estimated steps for normal emission: {n_steps_estimate}")
    
    def _fix_singularities(self, B_grid):
        """Fix singularities in magnetic field grid"""
        bad_mask = np.isnan(B_grid) | np.isinf(B_grid)
        
        if np.any(bad_mask):
            total_bad = np.sum(bad_mask.any(axis=-1))
            total_points = B_grid.shape[0] * B_grid.shape[1] * B_grid.shape[2]
            print(f"  Found {total_bad} grid points with singularities ({total_bad/total_points*100:.1f}%)")
            B_grid[bad_mask] = 0.0
            print("  Replaced singular points with zeros")
    
    def _bfield_loop_offaxis(self, x, y, z, a, I):
        """Calculate magnetic field from a single current loop"""
        x = np.atleast_1d(x).astype(float)
        y = np.atleast_1d(y).astype(float)
        z = np.atleast_1d(z).astype(float)

        rho = np.sqrt(x**2 + y**2)
        B_rho = np.zeros_like(rho)
        B_z = np.zeros_like(rho)

        mask_rho_pos = (rho > 0)
        mask_rho_zero = ~mask_rho_pos

        if np.any(mask_rho_pos):
            rp = rho[mask_rho_pos]
            zp = z[mask_rho_pos]
            num = 4 * a * rp
            denom = np.maximum((a + rp)**2 + zp**2, 1e-12)
            k2 = np.clip(num / denom, 0.0, 1.0-1e-12)
            K = ellipk(k2)
            E = ellipe(k2)
            common = np.sqrt(denom)
            
            term1 = -K
            term2 = (a*a + rp*rp + zp*zp) / np.maximum((a - rp)**2 + zp**2, 1e-12) * E
            B_rho_p = (self.mu0 * I)/(2*np.pi) * (zp / (rp*common)) * (term1 + term2)
            
            term3 = K
            term4 = (a*a - rp*rp - zp*zp)/np.maximum((a-rp)**2 + zp**2, 1e-12) * E
            B_z_p = (self.mu0 * I)/(2*np.pi*common) * (term3 + term4)
            
            # Check for finite values
            finite_mask = np.isfinite(B_rho_p) & np.isfinite(B_z_p)
            B_rho_p = np.where(finite_mask, B_rho_p, 0.0)
            B_z_p = np.where(finite_mask, B_z_p, 0.0)
            
            B_rho[mask_rho_pos] = B_rho_p
            B_z[mask_rho_pos] = B_z_p

        if np.any(mask_rho_zero):
            zz_axis = z[mask_rho_zero]
            B_z_axis = (self.mu0*I*a**2)/(2*(a**2+zz_axis**2)**1.5)
            B_z[mask_rho_zero] = B_z_axis
            B_rho[mask_rho_zero] = 0.0

        Bx = np.zeros_like(B_rho)
        By = np.zeros_like(B_rho)
        mask = (rho > 1e-12)
        Bx[mask] = B_rho[mask] * (x[mask]/rho[mask])
        By[mask] = B_rho[mask] * (y[mask]/rho[mask])
        
        return Bx, By, B_z

    def _bfield_solenoid_ellip(self, x, y, z):
        """Calculate magnetic field from solenoid using elliptic integrals"""
        x = np.atleast_1d(x).astype(float)
        y = np.atleast_1d(y).astype(float)
        z = np.atleast_1d(z).astype(float)
        
        Bx = np.zeros_like(x)
        By = np.zeros_like(y)
        Bz = np.zeros_like(z)
        
        for R_loop in self.rr:
            for z_loop in self.zz:
                bx, by, bz = self._bfield_loop_offaxis(x, y, z - z_loop, R_loop, self.I)
                Bx += bx
                By += by
                Bz += bz
        
        return np.stack([Bx, By, Bz], axis=-1)

    def _bfield_solenoid_rot_trans(self, x, y, z):
        """Calculate magnetic field with coil rotation and translation
            Input: positions in original coordinates
            x (N^3,)
            y (N^3,)
            z (N^3,)
            Return: vectors in original coordinates (grid-like)
            B-field (N,N,N,3)
        """
        if len(x.shape) > 1:
            raise ValueError("Input arrays must be 1D.")
        xshape = x.shape[0]
        grid_size = round(xshape ** (1/3))
        is_perfect_cube = grid_size ** 3 == xshape
        if not is_perfect_cube:
            raise ValueError("Input arrays must represent a perfect cube, so that the output B can be grid-like (N,N,N,3)")
        
        x_trans = x - self.trans[0]
        y_trans = y - self.trans[1]
        z_trans = z - self.trans[2]
        points = np.vstack([x_trans, y_trans, z_trans])

        # Rotation matrices
        Rx = np.array([[1,0,0],
                       [0,np.cos(self.theta_x),-np.sin(self.theta_x)],
                       [0,np.sin(self.theta_x),np.cos(self.theta_x)]])
        Ry = np.array([[np.cos(self.theta_y),0,np.sin(self.theta_y)],
                       [0,1,0],
                       [-np.sin(self.theta_y),0,np.cos(self.theta_y)]])
        Rz = np.array([[np.cos(self.theta_z),-np.sin(self.theta_z),0],
                       [np.sin(self.theta_z),np.cos(self.theta_z),0],
                       [0,0,1]])

        points_new = np.dot(Rx.T, np.dot(Ry.T, np.dot(Rz.T, points)))

        x_new = points_new[0,:]
        y_new = points_new[1,:]
        z_new = points_new[2,:]

        B_new = self._bfield_solenoid_ellip(x_new, y_new, z_new)
        B_original = np.einsum('ij,jk->ik', Rz @ Ry @ Rx, B_new.T).T  # (N^3, 3)

        # Reshape back to original grid shape + vector dimension
        B_original = B_original.reshape(grid_size, grid_size, grid_size, 3)
        
        return B_original

    def boris_push(self, r, v, q=None, m=None, dt=None):
        """
        Boris pusher algorithm for charged particle trajectory
        
        Parameters:
        -----------
        r : array_like
            Position vector [x, y, z]
        v : array_like
            Velocity vector [vx, vy, vz]
        q : float, optional
            Charge (default: -e for electron)
        m : float, optional
            Mass (default: electron mass)
        dt : float, optional
            Time step (default: self.dt)
            
        Returns:
        --------
        r_new : ndarray
            New position
        v_new : ndarray
            New velocity
        """
        if q is None:
            q = -self.e
        if m is None:
            m = self.m
        if dt is None:
            dt = self.dt
            
        B = self.B_interp([r])[0]
        t = (q*B/m)*0.5*dt
        t_mag2 = np.dot(t,t)
        s = 2*t/(1+t_mag2)
        v_minus = v
        v_prime = v_minus + np.cross(v_minus,t)
        v_plus = v_minus + np.cross(v_prime,s)
        v_new = v_plus
        r_new = r + v_new*dt
        return r_new, v_new

    def simulate_trajectory(self, initial_position, initial_velocity, 
                          max_steps=None, track_trajectory=False):
        """
        Simulate single electron trajectory
        
        Parameters:
        -----------
        initial_position : array_like
            Initial position [x, y, z]
        initial_velocity : array_like
            Initial velocity [vx, vy, vz]
        max_steps : int, optional
            Maximum simulation steps
        track_trajectory : bool
            Whether to return full trajectory
            
        Returns:
        --------
        final_position : ndarray or None
            Final position if reached detector
        final_velocity : ndarray or None
            Final velocity if reached detector
        trajectory : list, optional
            Full trajectory if track_trajectory=True
        """
        if max_steps is None:
            max_steps = self.max_sim_steps
            
        r = np.array(initial_position, dtype=float)
        v = np.array(initial_velocity, dtype=float)
        
        trajectory = [r.copy()] if track_trajectory else None
        
        for n in range(max_steps):
            r, v = self.boris_push(r, v)
            
            if track_trajectory:
                trajectory.append(r.copy())
                
            if r[2] >= self.L_flight:
                if track_trajectory:
                    return r, v, trajectory
                else:
                    return r, v
        
        # Didn't reach detector
        if track_trajectory:
            return None, None, trajectory
        else:
            return None, None

    def simulate_k_grid(self, k_range=1.0, grid_points=10):
        """
        Simulate trajectories for a grid of initial k-vectors
        
        Parameters:
        -----------
        k_range : float
            k-space range in units of 1e10 m^-1
        grid_points : int
            Number of grid points per dimension
            
        Returns:
        --------
        results : dict
            Dictionary containing simulation results
        """
        # Create k-space grid
        kx_vals = np.linspace(-k_range, k_range, grid_points) * 1e10
        ky_vals = np.linspace(-k_range, k_range, grid_points) * 1e10
        kx0, ky0 = np.meshgrid(kx_vals, ky_vals)
        
        # Calculate kz from momentum conservation
        kz2 = self.kF**2 - kx0**2 - ky0**2
        kz0 = np.where(kz2 >= 0, np.sqrt(kz2), np.nan)
        
        # Convert to velocities
        vx0 = self.hbar * kx0 / self.m
        vy0 = self.hbar * ky0 / self.m
        vz0 = self.hbar * kz0 / self.m
        
        # Initial positions (all at origin)
        pos0 = np.zeros((*kx0.shape, 3))
        vel0 = np.stack([vx0, vy0, vz0], axis=-1)
        valid = ~np.isnan(vz0)
        
        # Storage for results
        final_pos = np.full_like(pos0, np.nan)
        final_vel = np.full_like(vel0, np.nan)
        
        print(f"Simulating {grid_points}x{grid_points} k-space grid...")
        start = time.time()
        
        for i in range(grid_points):
            for j in range(grid_points):
                if not valid[i,j]:
                    continue
                    
                r_final, v_final = self.simulate_trajectory(pos0[i,j], vel0[i,j])
                
                if r_final is not None:
                    final_pos[i,j] = r_final
                    final_vel[i,j] = v_final
        
        elapsed = time.time() - start
        print(f"Simulation completed in {elapsed:.1f} s")
        
        # Convert back to k-space
        final_kx = self.m * final_vel[...,0] / self.hbar / 1e10
        final_ky = self.m * final_vel[...,1] / self.hbar / 1e10
        initial_kx = kx0 / 1e10
        initial_ky = ky0 / 1e10
        
        return {
            'initial_kx': initial_kx,
            'initial_ky': initial_ky,
            'final_kx': final_kx,
            'final_ky': final_ky,
            'initial_pos': pos0,
            'initial_vel': vel0,
            'final_pos': final_pos,
            'final_vel': final_vel,
            'valid': valid,
            'simulation_time': elapsed
        }

    def get_final_velocity(self, initial_velocity):
        """
        Get final velocity for given initial velocity
        
        Parameters:
        -----------
        initial_velocity : array_like
            Initial velocity [vx, vy, vz] or array of velocities
            
        Returns:
        --------
        final_velocity : ndarray or None
            Final velocity if trajectory reaches detector
        """
        initial_velocity = np.atleast_2d(initial_velocity)
        
        if initial_velocity.shape[0] == 1:
            # Single trajectory
            r_final, v_final = self.simulate_trajectory([0, 0, 0], initial_velocity[0])
            return v_final
        else:
            # Multiple trajectories
            results = []
            for i in range(initial_velocity.shape[0]):
                r_final, v_final = self.simulate_trajectory([0, 0, 0], initial_velocity[i])
                results.append(v_final)
            return np.array(results)

    def build_reverse_mapping(self, v0_range=[6.23e6], vx_element_num=30, vy_element_num=30, 
                            delta_vx=3e6, delta_vy=3e6):
        """
        Build reverse mapping from final to initial velocities for different energy slices
        
        Parameters:
        -----------
        v0_range : array_like, optional
            Array of total velocity magnitudes for energy slices
        vx_element_num : int
            Number of vx grid points per slice
        vy_element_num : int
            Number of vy grid points per slice
        delta_vx : float
            vx range for each slice
        delta_vy : float
            vy range for each slice

        After running this function, self.rev_map_cell will be rebuilt to be a list of functions.
        Each function will map final velocities to initial velocities for a specific energy slice.
        """
        if v0_range is None:
            v0_range = np.linspace(4e6, 8e6, 8)
        
        self.v0_range = v0_range
        self.inv_map_cell = []
        
        print("Building energy-sliced reverse mapping...")
        
        for z_slice, v_total in enumerate(v0_range):
            print(f"  Processing energy slice {z_slice+1}/{len(v0_range)}, |v| = {v_total:.1e} m/s")
            
            # Create velocity grid for this energy slice
            vx_vals_slice = np.linspace(-delta_vx/2, delta_vx/2, vx_element_num)
            vy_vals_slice = np.linspace(-delta_vy/2, delta_vy/2, vy_element_num)
            vx_grid, vy_grid = np.meshgrid(vx_vals_slice, vy_vals_slice)
            
            # Calculate vz from energy constraint
            vz_squared = v_total**2 - vx_grid**2 - vy_grid**2
            vz_grid = np.where(vz_squared >= 0, np.sqrt(vz_squared), np.nan)
            
            # Validity mask
            valid_mask = ~np.isnan(vz_grid)
            rho_squared = vx_grid**2 + vy_grid**2
            valid_mask &= (rho_squared <= (0.8 * vz_grid)**2)
            
            if np.sum(valid_mask) < 4:
                print(f"    Insufficient valid points for slice {z_slice+1}, skipping...")
                self.inv_map_cell.append(None)
                continue
            
            # Extract valid initial velocities
            initial_vx_flat = vx_grid[valid_mask]
            initial_vy_flat = vy_grid[valid_mask]
            initial_vz_flat = vz_grid[valid_mask]
            
            # Simulate trajectories
            final_vx_list = []
            final_vy_list = []
            valid_final_mask = []
            
            for idx in range(len(initial_vx_flat)):
                initial_vel = [initial_vx_flat[idx], initial_vy_flat[idx], initial_vz_flat[idx]]
                r_final, v_final = self.simulate_trajectory([0, 0, 0], initial_vel)
                
                if v_final is not None:
                    final_vx_list.append(v_final[0])
                    final_vy_list.append(v_final[1])
                    valid_final_mask.append(True)
                else:
                    valid_final_mask.append(False)
            
            valid_final_mask = np.array(valid_final_mask)
            
            if np.sum(valid_final_mask) < 4:
                print(f"    Insufficient final trajectories for slice {z_slice+1}, skipping...")
                self.inv_map_cell.append(None)
                continue
            
            # Create interpolation
            final_vx_array = np.array(final_vx_list)[valid_final_mask]
            final_vy_array = np.array(final_vy_list)[valid_final_mask]
            initial_vx_success = initial_vx_flat[valid_final_mask]
            initial_vy_success = initial_vy_flat[valid_final_mask]
            
            try:
                final_points = np.column_stack([final_vx_array, final_vy_array])
                
                interp_vx = LinearNDInterpolator(final_points, initial_vx_success, fill_value=np.nan)
                interp_vy = LinearNDInterpolator(final_points, initial_vy_success, fill_value=np.nan)
                
                def create_inverse_map(interp_vx, interp_vy):
                    def inverse_map(vx_query, vy_query):
                        if np.isscalar(vx_query):
                            points = np.array([[vx_query, vy_query]])
                        else:
                            points = np.column_stack([vx_query, vy_query])
                        
                        vx0_result = interp_vx(points)
                        vy0_result = interp_vy(points)
                        
                        if np.isscalar(vx_query):
                            return np.array([vx0_result[0], vy0_result[0]])
                        else:
                            return np.column_stack([vx0_result, vy0_result])
                    
                    return inverse_map
                
                self.inv_map_cell.append(create_inverse_map(interp_vx, interp_vy))
                print(f"    Successfully created interpolant with {len(final_vx_array)} points")
                
            except Exception as e:
                print(f"    Failed to create interpolant for slice {z_slice+1}: {e}")
                self.inv_map_cell.append(None)
        
        print("Reverse mapping construction complete!")

    def get_initial_velocity(self, vx_final, vy_final, v_total):
        """
        Get initial velocity from final velocity and total speed using reverse mapping
        
        Parameters:
        -----------
        vx_final : float or array
            Final x-velocity component
        vy_final : float or array
            Final y-velocity component
        v_total : float or array
            Total velocity magnitude

        Returns: (vx0, vy0)
        --------
        vx0 : float or array
            Initial x-velocity component
        vy0 : float or array
            Initial y-velocity component
        """
        if self.inv_map_cell is None:
            raise ValueError("Reverse mapping not built. Call build_reverse_mapping() first.")
        
        if np.isscalar(vx_final):
            return self._get_initial_velocity_single(vx_final, vy_final, v_total)
        else:
            return self._get_initial_velocity_batch(vx_final, vy_final, v_total)

    def _get_initial_velocity_single(self, vx_final, vy_final, v_total):
        """Get initial velocity for single point"""
        # Find closest energy slice
        best_slice = np.argmin(np.abs(self.v0_range - v_total))
        
        if self.inv_map_cell[best_slice] is None:
            return np.nan, np.nan
        
        try:
            result = self.inv_map_cell[best_slice](vx_final, vy_final)
            return result[0], result[1]
        except:
            return np.nan, np.nan

    def _get_initial_velocity_batch(self, vx_final_array, vy_final_array, v_total_array):
        """Get initial velocity for multiple points"""
        vx0_array = np.full_like(vx_final_array, np.nan)
        vy0_array = np.full_like(vy_final_array, np.nan)
        
        for i in range(len(vx_final_array)):
            vx0, vy0 = self._get_initial_velocity_single(
                vx_final_array[i], vy_final_array[i], v_total_array[i])
            vx0_array[i] = vx0
            vy0_array[i] = vy0
        
        return vx0_array, vy0_array

    def plot_k_mapping(self, results, show_arrows=True, show_rotation_labels=True):
        """
        Plot k-space mapping results
        
        Parameters:
        -----------
        results : dict
            Results from simulate_k_grid()
        show_arrows : bool
            Whether to show displacement arrows
        show_rotation_labels : bool
            Whether to show rotation angle and distance change labels
        """
        plt.figure(figsize=(8, 8))
        
        initial_kx = results['initial_kx']
        initial_ky = results['initial_ky']
        final_kx = results['final_kx']
        final_ky = results['final_ky']
        
        mask = np.isfinite(final_kx) & np.isfinite(final_ky)
        
        if show_arrows:
            for i in range(initial_kx.shape[0]):
                for j in range(initial_kx.shape[1]):
                    if mask[i,j]:
                        b0, a0 = initial_kx[i,j], initial_ky[i,j]
                        bf, af = final_kx[i,j], final_ky[i,j]
                        dx, dy = bf - b0, af - a0
                        
                        plt.arrow(b0, a0, dx, dy,
                                head_width=0.02, head_length=0.02,
                                fc='blue', ec='blue', alpha=0.5)
                        
                        if show_rotation_labels:
                            u = np.array([b0, a0])
                            v = np.array([bf, af])
                            if np.linalg.norm(u) > 1e-6 and np.linalg.norm(v) > 1e-6:
                                cosang = np.clip(np.dot(u, v) / (np.linalg.norm(u)*np.linalg.norm(v)), -1, 1)
                                rot_angle_deg = np.degrees(np.arccos(cosang))
                            else:
                                rot_angle_deg = 0.0
                            
                            r0 = np.linalg.norm(u)
                            rf = np.linalg.norm(v)
                            dr_percent = (rf - r0) / r0 * 100 if r0 > 1e-6 else 0.0
                            
                            plt.text(b0 + dx*0.5, a0 + dy*0.5,
                                   f'{rot_angle_deg:.1f}°\n{dr_percent:+.1f}%',
                                   fontsize=6, color='red', alpha=0.7,
                                   ha='center', va='center')
        
        plt.scatter(initial_kx, initial_ky, s=20, c='0.6', marker='o', label='Initial k (B=0)')
        
        plt.xlabel(r'$k_x\ (\mathrm{\AA^{-1}})$')
        plt.ylabel(r'$k_y\ (\mathrm{\AA^{-1}})$')
        plt.title(f'ARPES k-space mapping under solenoid B-field\n'
                 f'Coil: ID={self.coil_inner_d*1e3:.1f}mm, OD={self.coil_outer_d*1e3:.1f}mm, '
                 f'h={self.coil_height*1e3:.1f}mm, I={self.I:.2f}A\n'
                 f'hν={self.hv:.1f}eV, B(z=0)={self.B_mag_center*1e3:.3f} mT')
        plt.grid(True)
        plt.legend()
        plt.gca().set_aspect('equal')
        plt.tight_layout()
        plt.show()

    def __repr__(self):
        return (f"MagnetoARPES(hν={self.hv:.1f}eV, φ={self.phi:.1f}eV, "
                f"I={self.I:.2f}A, B_center={self.B_mag_center*1e3:.1f}mT)")

