mod generator;
mod simulator;

use crate::generator::Generator;
use crate::simulator::Simulator;
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
    commands.spawn((
        Camera3d::default(),
        Transform::from_xyz(10.0, 10.0, 10.0).looking_at(Vec3::ZERO, Vec3::Y),
    ));
    commands.spawn((
        PointLight {
            shadows_enabled: true,
            ..default()
        },
        Transform::from_xyz(4.0, 8.0, 4.0),
    ));

    let radius = 0.05;

    // Particle mesh and material
    let sphere_mesh = meshes.add(Mesh::from(Sphere::new(radius)));

    let particle_material = materials.add(StandardMaterial {
        base_color: Color::srgb(0.8, 0.2, 0.2),
        ..default()
    });
    let boundary_material = materials.add(StandardMaterial {
        base_color: Color::srgb(0.2, 0.8, 0.2),
        ..default()
    });

    let generator = Generator::new(radius);
    let mut simulator = Simulator::new(
        0.01,
        Vec3::new(0.0, -9.81, 0.0),
        simulator::FluidParams {
            kernel_radius: radius * 4.0,
            target_density: 1.0 / (2.0 * radius).powi(3),
            viscosity: 0.02,
            surface_tension: 1.0,
            adhesion: 10.0,
        },
    );

    let particles = generator.aabb(Vec3::new(-0.5, 0.0, -0.5), Vec3::new(0.5, 1.0, 0.5));
    particles.iter().enumerate().for_each(|(i, &p)| {
        simulator.add_particle(p);
        commands.spawn((
            Mesh3d(sphere_mesh.clone()),
            MeshMaterial3d(particle_material.clone()),
            Transform::from_translation(p),
            Index(i),
        ));
    });

    let boundaries = generator.open_box(Vec3::new(-2.0, -1.0, -2.0), Vec3::new(2.0, 3.0, 2.0));
    boundaries.iter().for_each(|&p| {
        simulator.add_boundary(p);
        // commands.spawn((
        //     Mesh3d(sphere_mesh.clone()),
        //     MeshMaterial3d(boundary_material.clone()),
        //     Transform::from_translation(p),
        // ));
    });
    simulator.init();
    commands.spawn(SimulatorComponent(simulator));
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
