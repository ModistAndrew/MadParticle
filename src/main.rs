use bevy::prelude::*;
use rand::Rng;

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

fn main() {
    App::new()
        .add_plugins(DefaultPlugins)
        .insert_resource(SimParams {
            gravity: Vec3::new(0.0, -10.0, 0.0),
            substeps: 5,
            collision_radius: 0.12, // Slightly larger than visual radius
            bounce_factor: 0.3,
            damping: 0.98,
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
        Transform::from_xyz(0.0, 5.0, 15.0).looking_at(Vec3::ZERO, Vec3::Y),
    ));

    // Add light
    commands.spawn((
        PointLight {
            shadows_enabled: true,
            ..default()
        },
        Transform::from_xyz(4.0, 8.0, 4.0),
    ));

    // Ground plane
    // commands.spawn((
    //     Mesh3d(meshes.add(Plane3d::new(
    //         Vec3::new(0.0, 1.0, 0.0),
    //         Vec2::new(20.0, 20.0),
    //     ))),
    //     MeshMaterial3d(materials.add(Color::srgb(0.3, 0.5, 0.3))),
    //     Transform::from_xyz(0.0, -0.9, 0.0),
    // ));

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

    for i in 0..1000 {
        let x = i % 10;
        let y = (i / 10) % 10;
        let z = i / 100;
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

    // Create some connections between nearby particles
    // for i in 0..particle_entities.len() - 1 {
    //     // Only create some connections
    //     commands.spawn(Connection {
    //         particle1: particle_entities[i],
    //         particle2: particle_entities[i + 1],
    //         rest_length: 0.5,
    //         stiffness: 0.8,
    //     });
    // }
}

fn pbd_step(
    mut particles: Query<(Entity, &mut Particle, &mut Transform)>,
    connections: Query<&Connection>,
    sim_params: Res<SimParams>,
    time: Res<Time>,
) {
    let dt = time.delta_secs();
    let substep_dt = dt / (sim_params.substeps as f32);

    if substep_dt <= 0.0 {
        return; // No time to simulate
    }

    for _ in 0..sim_params.substeps {
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
            }

            // 2b. Solve distance constraints (connections)
            for connection in connections.iter() {
                solve_distance_constraint(&mut particles, connection);
            }

            // 2c. Solve collision constraints between particles
            // solve_collision_constraints(&mut particles, sim_params.collision_radius);
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
}

fn solve_distance_constraint(
    particles: &mut Query<(Entity, &mut Particle, &mut Transform)>,
    connection: &Connection,
) {
    // Get both particles
    let (entity1, entity2) = (connection.particle1, connection.particle2);

    // Safety check - make sure we can get both particles
    let (mut p1, mut p2) = {
        let mut p1_opt = None;
        let mut p2_opt = None;

        for (entity, particle, _) in particles.iter() {
            if entity == entity1 {
                p1_opt = Some((entity, particle.predicted_position, particle.inverse_mass));
            } else if entity == entity2 {
                p2_opt = Some((entity, particle.predicted_position, particle.inverse_mass));
            }
        }

        if let (Some(p1), Some(p2)) = (p1_opt, p2_opt) {
            (p1, p2)
        } else {
            return; // One of the entities wasn't found
        }
    };

    // Calculate current distance
    let delta = p2.1 - p1.1;
    let current_distance = delta.length();

    // Skip if particles are at the same position
    if current_distance < 0.0001 {
        return;
    }

    // Calculate correction
    let direction = delta.normalize();
    let correction = (current_distance - connection.rest_length) * connection.stiffness;

    // Apply weighted correction based on inverse masses
    let total_inverse_mass = p1.2 + p2.2;
    if total_inverse_mass <= 0.0 {
        return; // Both particles are fixed
    }

    let p1_correction = (p1.2 / total_inverse_mass) * correction * direction;
    let p2_correction = -(p2.2 / total_inverse_mass) * correction * direction;

    // Apply corrections
    for (entity, mut particle, _) in particles.iter_mut() {
        if entity == p1.0 {
            particle.predicted_position += p1_correction;
        } else if entity == p2.0 {
            particle.predicted_position += p2_correction;
        }
    }
}

fn solve_collision_constraints(
    particles: &mut Query<(Entity, &mut Particle, &mut Transform)>,
    collision_radius: f32,
) {
    let particle_data: Vec<_> = particles
        .iter()
        .map(|(entity, particle, _)| (entity, particle.predicted_position, particle.inverse_mass))
        .collect();

    for i in 0..particle_data.len() {
        for j in (i + 1)..particle_data.len() {
            let (entity_i, pos_i, inv_mass_i) = particle_data[i];
            let (entity_j, pos_j, inv_mass_j) = particle_data[j];

            let delta = pos_j - pos_i;
            let distance = delta.length();
            let min_distance = collision_radius * 2.0;

            // Check if colliding
            if distance < min_distance && distance > 0.001 {
                let direction = delta.normalize();
                let correction = (min_distance - distance) * 0.5;

                // Calculate weighted corrections
                let total_inverse_mass = inv_mass_i + inv_mass_j;
                if total_inverse_mass <= 0.0 {
                    continue; // Both particles are fixed
                }

                let i_correction = -(inv_mass_i / total_inverse_mass) * correction * direction;
                let j_correction = (inv_mass_j / total_inverse_mass) * correction * direction;

                // Apply corrections
                for (entity, mut particle, _) in particles.iter_mut() {
                    if entity == entity_i {
                        particle.predicted_position += i_correction;
                    } else if entity == entity_j {
                        particle.predicted_position += j_correction;
                    }
                }
            }
        }
    }
}
