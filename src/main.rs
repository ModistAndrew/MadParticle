use bevy::prelude::*;
use std::collections::HashMap;
use std::f32::consts::PI;

// Components
#[derive(Component)]
struct Particle {
    old_position: Vec3,
    predicted_position: Vec3,
    velocity: Vec3,
    mass: f32,
    inverse_mass: f32, // Pre-computed for efficiency
}

// Optional connection between particles
#[derive(Component)]
struct Connection {
    particle1: Entity,
    particle2: Entity,
    rest_length: f32,
    stiffness: f32,
}

// Simulation parameters
#[derive(Resource)]
struct SimParams {
    gravity: Vec3,
    substeps: usize,
    collision_radius: f32,
    bounce_factor: f32,
    damping: f32,
}

// Fluid simulation parameters
#[derive(Resource)]
struct FluidParams {
    kernel_radius: f32,
    target_density: f32,
    pressure_multiplier: f32,
    spatial_hash_cell_size: f32,
}

fn main() {
    App::new()
        .add_plugins(DefaultPlugins)
        .insert_resource(SimParams {
            gravity: Vec3::new(0.0, -1.0, 0.0),
            substeps: 5,
            collision_radius: 0.12, // Slightly larger than visual radius
            bounce_factor: 0.3,
            damping: 0.98,
        })
        .insert_resource(FluidParams {
            kernel_radius: 0.5,
            target_density: 100.0,
            pressure_multiplier: 1.0,
            spatial_hash_cell_size: 0.5, // Should be at least equal to kernel_radius
        })
        .add_systems(Startup, setup)
        .add_systems(Update, pbd_step)
        .run();
}

fn setup(
    mut commands: Commands,
    mut meshes: ResMut<Assets<Mesh>>,
    mut materials: ResMut<Assets<StandardMaterial>>,
) {
    // Add camera
    commands.spawn((
        Camera3d::default(),
        Transform::from_xyz(10.0, 10.0, 10.0).looking_at(Vec3::ZERO, Vec3::Y),
    ));

    // Add light
    commands.spawn((
        PointLight {
            shadows_enabled: true,
            ..default()
        },
        Transform::from_xyz(4.0, 8.0, 4.0),
    ));

    // Particle mesh and material
    let sphere_mesh = meshes.add(Mesh::from(Sphere::new(0.1)));

    let material = materials.add(StandardMaterial {
        base_color: Color::srgb(0.8, 0.2, 0.2),
        metallic: 0.7,
        perceptual_roughness: 0.2,
        ..default()
    });

    // Spawn particles
    let mut rng = rand::rng();
    let mut particle_entities = Vec::new();

    for i in 0..5000 {
        let x = i % 10 - 5;
        let z = (i / 10) % 10 - 5;
        let y = i / 100;
        let position = Vec3::new(x as f32 * 0.2, y as f32 * 0.2, z as f32 * 0.2);

        let mass = 1.0;
        let inverse_mass = 1.0 / mass;

        let entity = commands
            .spawn((
                Mesh3d(sphere_mesh.clone()),
                MeshMaterial3d(material.clone()),
                Transform::from_translation(position),
                Particle {
                    old_position: position,
                    predicted_position: position,
                    velocity: Vec3::ZERO,
                    mass,
                    inverse_mass,
                },
            ))
            .id();

        particle_entities.push(entity);
    }
}

fn pbd_step(
    mut particles: Query<(Entity, &mut Particle, &mut Transform)>,
    connections: Query<&Connection>,
    sim_params: Res<SimParams>,
    fluid_params: Res<FluidParams>,
    time: Res<Time>,
) {
    let substep_dt = 0.04;

    // 1. Apply external forces and predict positions
    for (_, mut particle, transform) in particles.iter_mut() {
        particle.old_position = transform.translation;
        particle.velocity += sim_params.gravity * substep_dt;
        particle.predicted_position = transform.translation + particle.velocity * substep_dt;
    }

    // 2. Constraint resolution iterations
    for _ in 0..3 {
        // Multiple iterations for better convergence
        // 2a. Solve boundary constraints
        for (_, mut particle, _) in particles.iter_mut() {
            // Floor collision
            if particle.predicted_position.y < -1.0 {
                particle.predicted_position.y = -1.0;
            }

            let bounds = Vec3::new(2.0, 2.0, 2.0);
            particle.predicted_position.x =
                particle.predicted_position.x.clamp(-bounds.x, bounds.x);
            particle.predicted_position.z =
                particle.predicted_position.z.clamp(-bounds.z, bounds.z);
        }

        solve_fluid_density_constraint(&mut particles, &fluid_params);
    }

    // 3. Update velocities and positions
    for (_, mut particle, mut transform) in particles.iter_mut() {
        // Update velocity based on position change
        particle.velocity = (particle.predicted_position - particle.old_position) / substep_dt;

        // Apply damping
        particle.velocity *= sim_params.damping;

        // Update position
        transform.translation = particle.predicted_position;
    }
}

// Fluid density constraint solver function
fn solve_fluid_density_constraint(
    particles: &mut Query<(Entity, &mut Particle, &mut Transform)>,
    fluid_params: &FluidParams,
) {
    // Create a spatial hash for efficient neighbor finding
    let mut spatial_hash: HashMap<(i32, i32, i32), Vec<(Entity, Vec3, f32)>> = HashMap::new();

    // Collect all fluid particles into the spatial hash
    let fluid_particles_data: Vec<_> = particles
        .iter()
        .map(|(entity, particle, _)| {
            let pos = particle.predicted_position;
            let inv_mass = particle.inverse_mass;

            // Hash the position
            let cell_size = fluid_params.spatial_hash_cell_size;
            let cell_x = (pos.x / cell_size).floor() as i32;
            let cell_y = (pos.y / cell_size).floor() as i32;
            let cell_z = (pos.z / cell_size).floor() as i32;

            spatial_hash
                .entry((cell_x, cell_y, cell_z))
                .or_insert_with(Vec::new)
                .push((entity, pos, inv_mass));

            (entity, pos, inv_mass)
        })
        .collect();

    // Calculate density and density gradient for each fluid particle
    let mut lambdas: HashMap<Entity, f32> = HashMap::new();
    let mut density_gradients: HashMap<Entity, HashMap<Entity, Vec3>> = HashMap::new();

    for (entity_i, pos_i, _) in &fluid_particles_data {
        let mut density = 0.0;
        let mut self_gradient = Vec3::ZERO;
        let mut neighbor_gradients: HashMap<Entity, Vec3> = HashMap::new();

        // Get cell coordinates
        let cell_size = fluid_params.spatial_hash_cell_size;
        let cell_x = (pos_i.x / cell_size).floor() as i32;
        let cell_y = (pos_i.y / cell_size).floor() as i32;
        let cell_z = (pos_i.z / cell_size).floor() as i32;

        // Check neighboring cells and the cell containing the particle
        for dx in -1..=1 {
            for dy in -1..=1 {
                for dz in -1..=1 {
                    let neighbor_cell = (cell_x + dx, cell_y + dy, cell_z + dz);

                    if let Some(particles_in_cell) = spatial_hash.get(&neighbor_cell) {
                        for (entity_j, pos_j, _) in particles_in_cell {
                            let dir = pos_i - pos_j;
                            let r = dir.length();
                            if r < fluid_params.kernel_radius && r > 0.0001 {
                                // Calculate kernel value using cubic spline kernel
                                let kernel_value = cubic_kernel(fluid_params.kernel_radius, dir);

                                // Add to density
                                density += kernel_value;

                                // Calculate kernel gradient
                                let mut kernel_grad =
                                    cubic_kernel_gradient(fluid_params.kernel_radius, dir);

                                kernel_grad /= fluid_params.target_density;

                                self_gradient += kernel_grad;
                                if entity_i != entity_j {
                                    neighbor_gradients.insert(*entity_j, kernel_grad);
                                }
                            }
                        }
                    }
                }
            }
        }

        // Calculate constraint
        let constraint = (density / fluid_params.target_density - 1.0).max(0.0);

        // Calculate lambda (Lagrange multiplier)
        let mut sum_gradients_squared = 0.0;

        // Self gradient squared
        sum_gradients_squared += self_gradient.length_squared();

        // Neighbor gradients squared
        for (_, grad) in &neighbor_gradients {
            sum_gradients_squared += grad.length_squared();
        }

        let lambda = -constraint / (sum_gradients_squared + 0.001); // Add small epsilon for stability

        lambdas.insert(*entity_i, lambda);
        density_gradients.insert(*entity_i, neighbor_gradients);
    }

    // Apply position corrections based on density constraint
    let mut position_corrections: HashMap<Entity, Vec3> = HashMap::new();

    for (entity_i, _, _) in &fluid_particles_data {
        if let Some(lambda_i) = lambdas.get(entity_i) {
            if let Some(neighbor_gradients) = density_gradients.get(entity_i) {
                let mut correction = Vec3::ZERO;

                // Self position correction
                for (entity_j, grad_j) in neighbor_gradients {
                    if let Some(lambda_j) = lambdas.get(entity_j) {
                        correction += (lambda_i + lambda_j) * grad_j;
                    }
                }

                correction /= fluid_params.target_density;
                correction *= fluid_params.pressure_multiplier;
                position_corrections.insert(*entity_i, correction);
            }
        }
    }

    // Apply the calculated corrections
    for (entity, mut particle, _) in particles.iter_mut() {
        if let Some(correction) = position_corrections.get(&entity) {
            particle.predicted_position += correction;
        }
    }
}

// Cubic spline kernel function for SPH
fn cubic_kernel(h: f32, dir: Vec3) -> f32 {
    let r = dir.length();
    if r >= h {
        return 0.0;
    }
    let coeff = 15.0 / PI / h.powi(6);
    coeff * (h - r).powi(3)
}

// Gradient of the cubic spline kernel
fn cubic_kernel_gradient(h: f32, dir: Vec3) -> Vec3 {
    let r = dir.length();
    if r <= 0.0001 || r >= h {
        return Vec3::ZERO;
    }
    let coeff = -45.0 / PI / h.powi(6);
    coeff * (h - r).powi(2) * dir.normalize()
}
