use bevy::math::Vec3;
use itertools::Itertools;
use rayon::prelude::*;
use std::collections::HashMap;
use std::f32::consts::PI;

pub struct Particle {
    boundary: bool,
    old_position: Vec3,
    pub position: Vec3,
    velocity: Vec3,
    neighbors: Vec<usize>,
    density: f32,
    lambda: f32,
    phi: f32,
    normal: Vec3,
    force: Vec3,
}

impl Particle {
    fn new(boundary: bool, position: Vec3) -> Self {
        Self {
            boundary,
            old_position: position,
            position,
            velocity: Vec3::ZERO,
            neighbors: Vec::new(),
            density: 0.0,
            lambda: 0.0, // always 0.0 for boundary particles
            phi: 1.0,    // always 1.0 for fluid particles
            normal: Vec3::ZERO,
            force: Vec3::ZERO,
        }
    }
}

#[derive(Clone, Copy)]
pub struct FluidParams {
    pub kernel_radius: f32,
    pub target_density: f32,
    pub viscosity: f32,
    pub surface_tension: f32,
    pub adhesion: f32,
}

impl FluidParams {
    pub fn kernel(&self, dir: Vec3) -> f32 {
        let h = self.kernel_radius;
        let r = dir.length();
        let c = 315.0 / 64.0 / PI / h.powi(3);
        if r >= h {
            return 0.0;
        }
        c * (1.0 - r.powi(2) / h.powi(2)).powi(3)
    }

    pub fn kernel_gradient(&self, dir: Vec3) -> Vec3 {
        let h = self.kernel_radius;
        let r = dir.length();
        let c = -45.0 / PI / h.powi(6);
        if r >= h {
            return Vec3::ZERO;
        }
        c * (h - r).powi(2) * dir.normalize()
    }

    fn cohesion_kernel(&self, dir: Vec3) -> f32 {
        let h = self.kernel_radius;
        let r = dir.length();
        let c = 32.0 / PI / h.powi(9);
        c * if r > 0.0 && 2.0 * r <= h {
            2.0 * (h - r).powi(3) * r.powi(3) - h.powi(6) / 64.0
        } else if 2.0 * r > h && r <= h {
            (h - r).powi(3) * r.powi(3)
        } else {
            0.0
        }
    }

    fn adhesion_kernel(&self, dir: Vec3) -> f32 {
        let h = self.kernel_radius;
        let r = dir.length();
        let c = 0.007 / h.powf(3.25);
        c * if 2.0 * r > h && r <= h {
            (-4.0 * r.powi(2) / h + 6.0 * r - 2.0 * h).powf(0.25)
        } else {
            0.0
        }
    }
}

pub struct Simulator {
    pub particles: Vec<Particle>,
    boundaries: Vec<Particle>,
    step_dt: f32,
    gravity: Vec3,
    fluid_params: FluidParams,
}

impl Simulator {
    pub fn new(step_dt: f32, gravity: Vec3, fluid_params: FluidParams) -> Self {
        Self {
            particles: Vec::new(),
            boundaries: Vec::new(),
            step_dt,
            gravity,
            fluid_params,
        }
    }

    pub fn add_particle(&mut self, position: Vec3) {
        self.particles.push(Particle::new(false, position));
    }

    pub fn add_boundary(&mut self, position: Vec3) {
        self.boundaries.push(Particle::new(true, position));
    }

    pub fn step(&mut self) {
        const SOLVER_ITERATIONS: usize = 5;
        let current_time = std::time::Instant::now();

        for particle in &mut self.particles {
            particle.force = self.gravity;
        }
        self.update_neighbors();
        self.update_densities();
        self.update_normal();
        self.apply_surface_tension();
        self.apply_adhesion();
        for particle in &mut self.particles {
            particle.velocity += particle.force * self.step_dt;
        }
        for particle in &mut self.particles {
            particle.position = particle.old_position + particle.velocity * self.step_dt;
        }

        for _ in 0..SOLVER_ITERATIONS {
            self.update_neighbors();
            self.update_densities();
            self.update_lambdas();
            self.update_deltas();
        }

        for particle in &mut self.particles {
            particle.velocity = (particle.position - particle.old_position) / self.step_dt;
        }
        self.apply_viscosity();
        for particle in &mut self.particles {
            particle.old_position = particle.position;
        }

        let elapsed = current_time.elapsed();
        println!("PBD step took: {:.2?}", elapsed);
    }

    pub fn get_particle_position(&self, index: usize) -> Vec3 {
        self.particles[index].position
    }

    pub fn init(&mut self) {
        self.compute_phis();
    }

    fn update_normal(&mut self) {
        let (positions, boundaries, densities): (Vec<_>, Vec<_>, Vec<_>) = self
            .particle_chain()
            .map(|p| (p.position, p.boundary, p.density))
            .multiunzip();
        let fluid_params = self.fluid_params;

        self.particles.par_iter_mut().for_each(|particle| {
            let mut normal = Vec3::ZERO;
            for &j in &particle.neighbors {
                if boundaries[j] {
                    continue;
                }
                let dir = particle.position - positions[j];
                normal += self.fluid_params.kernel_gradient(dir) / densities[j];
            }
            normal *= fluid_params.kernel_radius;
            particle.normal = normal;
        });
    }

    fn apply_surface_tension(&mut self) {
        let (positions, boundaries, normals, densities): (Vec<_>, Vec<_>, Vec<_>, Vec<_>) = self
            .particle_chain()
            .map(|p| (p.position, p.boundary, p.normal, p.density))
            .multiunzip();
        let fluid_params = self.fluid_params;

        self.particles.par_iter_mut().for_each(|particle| {
            for &j in &particle.neighbors {
                if boundaries[j] {
                    continue;
                }
                let dir = particle.position - positions[j];
                let cohesion_force = -fluid_params.surface_tension
                    * fluid_params.cohesion_kernel(dir)
                    * dir.normalize();
                let curvature_force =
                    -fluid_params.surface_tension * (particle.normal - normals[j]);
                let k = 2.0 * fluid_params.target_density / (particle.density + densities[j]);
                particle.force += k * (cohesion_force + curvature_force);
            }
        });
    }

    fn apply_adhesion(&mut self) {
        let (positions, boundaries, phis): (Vec<_>, Vec<_>, Vec<_>) = self
            .particle_chain()
            .map(|p| (p.position, p.boundary, p.phi))
            .multiunzip();
        let fluid_params = self.fluid_params;

        self.particles.par_iter_mut().for_each(|particle| {
            for &j in &particle.neighbors {
                if !boundaries[j] {
                    continue;
                }
                let dir = particle.position - positions[j];
                particle.force += -fluid_params.adhesion
                    * phis[j]
                    * fluid_params.adhesion_kernel(dir)
                    * dir.normalize();
            }
        });
    }

    fn apply_viscosity(&mut self) {
        let (positions, velocities, boundaries, densities): (Vec<_>, Vec<_>, Vec<_>, Vec<_>) = self
            .particle_chain()
            .map(|p| (p.position, p.velocity, p.boundary, p.density))
            .multiunzip();
        let fluid_params = self.fluid_params;

        self.particles.par_iter_mut().for_each(|particle| {
            for &neighbor_index in &particle.neighbors {
                if boundaries[neighbor_index] {
                    continue;
                }
                let kernel_value =
                    fluid_params.kernel(particle.position - positions[neighbor_index]);
                let relative_velocity = particle.velocity - velocities[neighbor_index];
                particle.velocity -= fluid_params.viscosity / densities[neighbor_index]
                    * kernel_value
                    * relative_velocity;
            }
        });
    }

    fn compute_phis(&mut self) {
        println!("Precomputing boundary phi values...");
        let radius = self.fluid_params.kernel_radius;
        let cell_size = radius;
        let boundary_positions: Vec<Vec3> = self.boundaries.iter().map(|p| p.position).collect();
        let mut grid: HashMap<(i32, i32, i32), Vec<usize>> = HashMap::new();
        for (i, &position) in boundary_positions.iter().enumerate() {
            let cell = (
                (position.x / cell_size).floor() as i32,
                (position.y / cell_size).floor() as i32,
                (position.z / cell_size).floor() as i32,
            );
            grid.entry(cell).or_default().push(i);
        }
        let fluid_params = self.fluid_params;
        self.boundaries.par_iter_mut().for_each(|particle| {
            let position = particle.position;
            let cell = (
                (position.x / cell_size).floor() as i32,
                (position.y / cell_size).floor() as i32,
                (position.z / cell_size).floor() as i32,
            );
            let mut neighbors = Vec::new();
            for dx in -1..=1 {
                for dy in -1..=1 {
                    for dz in -1..=1 {
                        let neighbor_cell = (cell.0 + dx, cell.1 + dy, cell.2 + dz);
                        if let Some(cell_particles) = grid.get(&neighbor_cell) {
                            for &neighbor_index in cell_particles {
                                let dist = position.distance(boundary_positions[neighbor_index]);
                                if dist < radius {
                                    neighbors.push(neighbor_index);
                                }
                            }
                        }
                    }
                }
            }
            let mut d = 0.0;
            for &neighbor_index in &neighbors {
                let dir = position - boundary_positions[neighbor_index];
                d += fluid_params.kernel(dir);
            }
            particle.phi = fluid_params.target_density / d;
        });
        println!("Boundary phi computation completed.");
    }

    fn update_neighbors(&mut self) {
        let radius = self.fluid_params.kernel_radius;
        let cell_size = radius;
        let positions: Vec<Vec3> = self.particle_chain().map(|p| p.position).collect();
        let mut grid: HashMap<(i32, i32, i32), Vec<usize>> = HashMap::new();
        for (i, &position) in positions.iter().enumerate() {
            let cell = (
                (position.x / cell_size).floor() as i32,
                (position.y / cell_size).floor() as i32,
                (position.z / cell_size).floor() as i32,
            );
            grid.entry(cell).or_default().push(i);
        }
        self.particles
            .par_iter_mut()
            .enumerate()
            .for_each(|(i, particle)| {
                let position = particle.position;
                let cell = (
                    (position.x / cell_size).floor() as i32,
                    (position.y / cell_size).floor() as i32,
                    (position.z / cell_size).floor() as i32,
                );
                let mut neighbors = Vec::new();
                for dx in -1..=1 {
                    for dy in -1..=1 {
                        for dz in -1..=1 {
                            let neighbor_cell = (cell.0 + dx, cell.1 + dy, cell.2 + dz);
                            if let Some(cell_particles) = grid.get(&neighbor_cell) {
                                for &neighbor_index in cell_particles {
                                    if neighbor_index == i {
                                        continue;
                                    }
                                    let dist = position.distance(positions[neighbor_index]);
                                    if dist < radius {
                                        neighbors.push(neighbor_index);
                                    }
                                }
                            }
                        }
                    }
                }
                particle.neighbors = neighbors;
            });
    }

    fn update_densities(&mut self) {
        let (positions, phis): (Vec<_>, Vec<_>) = self
            .particle_chain()
            .map(|p| (p.position, p.phi))
            .multiunzip();
        let fluid_params = self.fluid_params;

        self.particles.par_iter_mut().for_each(|particle| {
            let mut density = 0.0;
            for &neighbor_index in &particle.neighbors {
                let dir = particle.position - positions[neighbor_index];
                density += phis[neighbor_index] * fluid_params.kernel(dir);
            }
            density += fluid_params.kernel(Vec3::ZERO);
            particle.density = density;
        });
    }

    fn update_lambdas(&mut self) {
        const EPSILON: f32 = 1e-6;
        let (positions, boundaries, phis): (Vec<_>, Vec<_>, Vec<_>) = self
            .particle_chain()
            .map(|p| (p.position, p.boundary, p.phi))
            .multiunzip();
        let fluid_params = self.fluid_params;

        self.particles.par_iter_mut().for_each(|particle| {
            let mut sum_gradient_squared = 0.0;
            let mut self_gradient = Vec3::ZERO;
            for &neighbor_index in &particle.neighbors {
                let dir = particle.position - positions[neighbor_index];
                let gradient = phis[neighbor_index] * fluid_params.kernel_gradient(dir)
                    / fluid_params.target_density;
                if !boundaries[neighbor_index] {
                    sum_gradient_squared += gradient.length_squared();
                }
                self_gradient += gradient;
            }
            sum_gradient_squared += self_gradient.length_squared();
            let constraint = (particle.density / fluid_params.target_density - 1.0).max(0.0);
            particle.lambda = -constraint / (sum_gradient_squared + EPSILON);
        });
    }

    fn update_deltas(&mut self) {
        let (positions, phis, lambdas): (Vec<_>, Vec<_>, Vec<_>) = self
            .particle_chain()
            .map(|p| (p.position, p.phi, p.lambda))
            .multiunzip();
        let fluid_params = self.fluid_params;

        self.particles.par_iter_mut().for_each(|particle| {
            let mut correction = Vec3::ZERO;
            for &neighbor_index in &particle.neighbors {
                let dir = particle.position - positions[neighbor_index];
                let gradient = phis[neighbor_index] * fluid_params.kernel_gradient(dir)
                    / fluid_params.target_density;
                correction += (particle.lambda + lambdas[neighbor_index]) * gradient;
            }
            particle.position += correction;
        });
    }

    fn particle_chain(&self) -> impl Iterator<Item = &Particle> {
        self.particles.iter().chain(self.boundaries.iter())
    }
}
