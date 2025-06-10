use bevy::math::Vec3;
use rstar::primitives::GeomWithData;
use rstar::{PointDistance, RTree};
use std::f32::consts::PI;

const GRAVITY: Vec3 = Vec3::new(0.0, -1.0, 0.0);
const SOLVER_ITERATIONS: usize = 5;
const STEP_DT: f32 = 0.04;

const KERNEL_RADIUS: f32 = 0.5;
const TARGET_DENSITY: f32 = 100.0;
const VORTICITY: f32 = 0.1;
const VISCOSITY: f32 = 0.01;
const EPSILON: f32 = 0.0001;
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
        for particle in &mut self.particles {
            particle.velocity += GRAVITY * STEP_DT;
            particle.position = particle.old_position + particle.velocity * STEP_DT;
        }

        println!("Updating neighbors...");
        self.update_neighbors();
        println!("Finished updating neighbors.");

        for _ in 0..SOLVER_ITERATIONS {
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
    }

    pub fn get_particle_position(&self, index: usize) -> Vec3 {
        self.particles[index].position
    }

    fn update_neighbors(&mut self) {
        let tree = RTree::bulk_load(
            self.particles
                .iter()
                .enumerate()
                .map(|(i, particle)| GeomWithData::new(particle.position.to_array(), i))
                .collect(),
        );
        for particle in &mut self.particles {
            let position = particle.position.to_array();
            particle.neighbors = tree
                .locate_within_distance(position, KERNEL_RADIUS)
                .filter(|g| g.distance_2(&position) > EPSILON)
                .map(|g| g.data)
                .collect();
        }
    }

    fn update_lambdas(&mut self) {
        for i in 0..self.particles.len() {
            let particle = &self.particles[i];
            let mut density = 0.0;
            let mut sum_gradient_sqared = 0.0;
            let mut self_gradient = Vec3::ZERO;
            for &neighbor_index in &particle.neighbors {
                let neighbor = &self.particles[neighbor_index];
                let dir = particle.position - neighbor.position;
                density += Self::cubic_kernel(KERNEL_RADIUS, dir);
                let gradient = Self::cubic_kernel_gradient(KERNEL_RADIUS, dir) / TARGET_DENSITY;
                sum_gradient_sqared += gradient.length_squared();
                self_gradient += gradient;
            }
            sum_gradient_sqared += self_gradient.length_squared();
            let constraint = (density / TARGET_DENSITY - 1.0).max(0.0);
            self.particles[i].lambda = -constraint / (sum_gradient_sqared + EPSILON);
        }
    }

    fn update_deltas(&mut self) {
        for i in 0..self.particles.len() {
            let particle = &self.particles[i];
            let mut correction = Vec3::ZERO;
            for &neighbor_index in &particle.neighbors {
                let neighbor = &self.particles[neighbor_index];
                let dir = particle.position - neighbor.position;
                let gradient = Self::cubic_kernel_gradient(KERNEL_RADIUS, dir);
                correction += (particle.lambda + neighbor.lambda) * gradient;
            }
            correction /= TARGET_DENSITY;
            self.particles[i].position += correction;
        }
    }

    fn cubic_kernel(h: f32, dir: Vec3) -> f32 {
        let r = dir.length();
        if r >= h {
            return 0.0;
        }
        let coeff = 15.0 / PI / h.powi(6);
        coeff * (h - r).powi(3)
    }

    fn cubic_kernel_gradient(h: f32, dir: Vec3) -> Vec3 {
        let r = dir.length();
        if r <= EPSILON || r >= h {
            return Vec3::ZERO;
        }
        let coeff = -45.0 / PI / h.powi(6);
        coeff * (h - r).powi(2) * dir.normalize()
    }
}
