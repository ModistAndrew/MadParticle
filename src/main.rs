mod simulator;
use crate::simulator::{Particle, Simulator};
use bevy::prelude::*;
#[derive(Component)]
struct Index(usize);
#[derive(Component)]
struct SimulatorComponent(Simulator);
fn main() {
    App::new()
        .add_plugins(DefaultPlugins)
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

    let mut particles = Vec::new();

    for i in 0..5000 {
        let x = i % 10 - 5;
        let z = (i / 10) % 10 - 5;
        let y = i / 100;
        let position = Vec3::new(x as f32 * 0.2, y as f32 * 0.2, z as f32 * 0.2);

        particles.push(Particle::new(position));
        commands.spawn((
            Mesh3d(sphere_mesh.clone()),
            MeshMaterial3d(material.clone()),
            Transform::from_translation(position),
            Index(i as usize),
        ));
    }

    commands.spawn(SimulatorComponent(Simulator::new(particles)));
}

fn pbd_step(
    mut particles: Query<(&Index, &mut Transform)>,
    mut simulator: Single<&mut SimulatorComponent>,
) {
    simulator.0.step();
    for (id, mut transform) in particles.iter_mut() {
        transform.translation = simulator.0.get_particle_position(id.0);
    }
}
