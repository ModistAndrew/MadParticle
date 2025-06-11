use bevy::math::Vec3;
use rayon::prelude::*;
use std::collections::HashMap;
use std::f32::consts::PI;

const GRAVITY: Vec3 = Vec3::new(0.0, -9.81, 0.0);
const SOLVER_ITERATIONS: usize = 5;
const STEP_DT: f32 = 0.04;

const KERNEL_RADIUS: f32 = 0.5;
const TARGET_DENSITY: f32 = 200.0;
const VORTICITY: f32 = 0.1;
const VISCOSITY: f32 = 0.01;
const EPSILON: f32 = 0.00001;
const DAMPING: f32 = 0.98;

pub struct Particle {
    old_position: Vec3,
    position: Vec3,
    velocity: Vec3,
    neighbors: Vec<usize>,
    lambda: f32,
    delta: Vec3,
}

impl Particle {
    pub fn new(position: Vec3) -> Self {
        Self {
            old_position: position,
            position: Vec3::ZERO,
            velocity: Vec3::ZERO,
            neighbors: Vec::new(),
            lambda: 0.0,
            delta: Vec3::ZERO,
        }
    }
}

pub struct Simulator {
    particles: Vec<Particle>,
}

impl Simulator {
    pub fn new(particles: Vec<Particle>) -> Self {
        Self { particles }
    }

    pub fn step(&mut self) {
        let current_time = std::time::Instant::now();
        for particle in &mut self.particles {
            particle.velocity += GRAVITY * STEP_DT;
            particle.position = particle.old_position + particle.velocity * STEP_DT;
        }

        for _ in 0..SOLVER_ITERATIONS {
            self.update_neighbors();
            self.update_lambdas();
            self.update_deltas();
            for particle in &mut self.particles {
                particle.position += particle.delta;
                let min_x = -2.0;
                let min_y = -1.0;
                let min_z = -2.0;
                particle.position.y = particle.position.y.max(min_y);
                particle.position.x = particle.position.x.clamp(min_x, -min_x);
                particle.position.z = particle.position.z.clamp(min_z, -min_z);
            }
        }

        for particle in &mut self.particles {
            particle.velocity = (particle.position - particle.old_position) / STEP_DT;
            // particle.velocity *= DAMPING;
            particle.old_position = particle.position;
        }
        let elapsed = current_time.elapsed();
        println!("PBD step took: {:.2?}", elapsed);
    }

    pub fn get_particle_position(&self, index: usize) -> Vec3 {
        self.particles[index].position
    }

    fn update_neighbors(&mut self) {
        let cell_size = KERNEL_RADIUS;
        let positions: Vec<Vec3> = self.particles.iter().map(|p| p.position).collect();
        let mut grid: HashMap<(i32, i32, i32), Vec<usize>> = HashMap::new();
        for (i, &position) in positions.iter().enumerate() {
            let cell = (
                (position.x / cell_size).floor() as i32,
                (position.y / cell_size).floor() as i32,
                (position.z / cell_size).floor() as i32,
            );
            grid.entry(cell).or_default().push(i);
        }
        self.particles.par_iter_mut().for_each(|particle| {
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
                                let dist = position.distance(positions[neighbor_index]);
                                if dist < KERNEL_RADIUS && dist > EPSILON {
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

    fn update_lambdas(&mut self) {
        let positions: Vec<Vec3> = self.particles.iter().map(|p| p.position).collect();
        self.particles.par_iter_mut().for_each(|particle| {
            let mut density = 0.0;
            let mut sum_gradient_sqared = 0.0;
            let mut self_gradient = Vec3::ZERO;
            for &neighbor_index in &particle.neighbors {
                let dir = particle.position - positions[neighbor_index];
                density += Self::kernel(KERNEL_RADIUS, dir);
                let gradient = Self::kernel_gradient(KERNEL_RADIUS, dir) / TARGET_DENSITY;
                sum_gradient_sqared += gradient.length_squared();
                self_gradient += gradient;
            }
            density += Self::kernel(KERNEL_RADIUS, Vec3::ZERO);
            sum_gradient_sqared += self_gradient.length_squared();
            let constraint = (density / TARGET_DENSITY - 1.0).max(0.0);
            particle.lambda = -constraint / (sum_gradient_sqared + EPSILON);
        });
    }

    fn update_deltas(&mut self) {
        let positions: Vec<Vec3> = self.particles.iter().map(|p| p.position).collect();
        let lambdas: Vec<f32> = self.particles.iter().map(|p| p.lambda).collect();
        self.particles.par_iter_mut().for_each(|particle| {
            let mut correction = Vec3::ZERO;
            for &neighbor_index in &particle.neighbors {
                let dir = particle.position - positions[neighbor_index];
                let gradient = Self::kernel_gradient(KERNEL_RADIUS, dir) / TARGET_DENSITY;
                correction += (particle.lambda + lambdas[neighbor_index]) * gradient;
            }
            particle.delta = correction;
        });
    }

    fn kernel(h: f32, dir: Vec3) -> f32 {
        let r = dir.length();
        if r >= h {
            return 0.0;
        }
        let coeff = 315.0 / 64.0 / PI / h.powi(3);
        coeff * (1.0 - r.powi(2) / h.powi(2)).powi(3)
    }

    fn kernel_gradient(h: f32, dir: Vec3) -> Vec3 {
        let r = dir.length();
        if r <= EPSILON || r >= h {
            return Vec3::ZERO;
        }
        let coeff = -45.0 / PI / h.powi(6);
        coeff * (h - r).powi(2) * dir.normalize()
    }
}
