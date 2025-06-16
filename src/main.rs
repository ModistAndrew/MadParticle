mod aabb;
mod generator;
mod grid;
mod ref_surfacer;
mod simulator;
mod surfacer;
mod wrapper;

use crate::aabb::Aabb;
use crate::generator::Generator;
use crate::wrapper::{CommonParams, SimulatorParams, SurfacerParams, Wrapper};
use bevy::prelude::*;
use std::io::Write;

#[derive(Component)]
struct Index(usize);
#[derive(Component)]
struct WrapperComponent(Wrapper);
#[derive(Resource)]
struct SurfaceMeshHandle(Handle<Mesh>);

fn main() {
    App::new()
        .add_plugins(DefaultPlugins)
        .add_systems(Startup, setup)
        .add_systems(Update, pbd_step)
        .run();
}

fn init_wrapper() -> Wrapper {
    let radius = 0.05;
    let generator = Generator::new(radius);
    let particles = generator.aabb(Vec3::new(-1.0, 0.0, -1.0), Vec3::new(1.0, 2.0, 1.0));
    let boundaries = generator.open_box(Vec3::new(-2.0, -2.0, -2.0), Vec3::new(2.0, 4.0, 2.0));
    let wrapper = Wrapper::new(
        CommonParams {
            radius,
            viscosity: 0.02,
            surface_tension: 0.1,
            adhesion: 0.1,
        },
        SimulatorParams {
            step_dt: 1.0 / 120.0,
            gravity: Vec3::new(0.0, -9.81, 0.0),
            min_max: Aabb::new(Vec3::new(-25.0, -25.0, -25.0), Vec3::new(25.0, 25.0, 25.0)),
        },
        SurfacerParams {
            density_threshold: 0.7,
            min_max: Aabb::new(Vec3::new(-2.5, -2.5, -2.5), Vec3::new(2.5, 2.5, 2.5)),
        },
        particles,
        boundaries,
        100,
        "output/exported_mesh".to_string(),
    );
    wrapper
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

    let wrapper = init_wrapper();
    let radius = wrapper.radius;
    let sphere_mesh = meshes.add(Mesh::from(Sphere::new(radius)));
    let particle_material = materials.add(StandardMaterial {
        base_color: Color::srgb(0.8, 0.2, 0.2),
        ..default()
    });
    let boundary_material = materials.add(StandardMaterial {
        base_color: Color::srgb(0.2, 0.8, 0.2),
        ..default()
    });
    wrapper
        .simulator
        .positions()
        .iter()
        .enumerate()
        .for_each(|(i, &p)| {
            commands.spawn((
                Mesh3d(sphere_mesh.clone()),
                MeshMaterial3d(particle_material.clone()),
                Transform::from_translation(p),
                Index(i),
            ));
        });

    // wrapper
    //     .simulator
    //     .boundary_positions()
    //     .iter()
    //     .for_each(|&p| {
    //         commands.spawn((
    //             Mesh3d(sphere_mesh.clone()),
    //             MeshMaterial3d(boundary_material.clone()),
    //             Transform::from_translation(p),
    //         ));
    //     });
    commands.spawn(WrapperComponent(wrapper));
}

fn pbd_step(
    mut particles: Query<(&Index, &mut Transform)>,
    mut wrapper: Single<&mut WrapperComponent>,
) {
    wrapper.0.step();
    let positions = wrapper.0.simulator.positions();
    for (id, mut transform) in particles.iter_mut() {
        transform.translation = positions[id.0];
    }
}
