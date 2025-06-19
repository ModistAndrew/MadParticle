mod aabb;
mod config;
mod generator;
mod grid;
mod simulator;
mod surfacer;
mod wrapper;

use crate::aabb::Aabb;
use crate::config::Config;
use crate::generator::Generator;
use crate::wrapper::{CommonParams, SimulatorParams, SurfacerParams, Wrapper};
use bevy::prelude::*;
use serde_json::from_reader;
use std::env::args;
use std::fs::File;
use std::io::BufReader;

#[derive(Component)]
struct Index(usize);
#[derive(Component)]
struct WrapperComponent(Wrapper);
#[derive(Resource)]
struct ConfigResource(Config);

fn main() {
    let args: Vec<String> = args().collect();
    if args.len() < 2 {
        panic!("Usage: {} <config.json> [output_path]", args[0]);
    }
    let config_path = &args[1];
    let config = load_config(config_path);
    let headless = args.len() > 2;
    if headless {
        let output_path = &args[2];
        let mut wrapper = init_wrapper(&config, output_path);
        loop {
            wrapper.step();
        }
    } else {
        let config_resource = ConfigResource(config);
        App::new()
            .add_plugins(DefaultPlugins)
            .insert_resource(config_resource)
            .add_systems(Startup, setup)
            .add_systems(Update, pbd_step)
            .run();
    }
}
fn load_config(path: &str) -> Config {
    let file = File::open(path).expect("Failed to open config file");
    let reader = BufReader::new(file);
    from_reader(reader).expect("Failed to parse config file")
}

fn init_wrapper(config: &Config, output_path: &str) -> Wrapper {
    let radius = config.common.radius;
    let generator = Generator::new(radius);
    let mut particles = Vec::new();
    for particle_gen in &config.particle_generators {
        particles.extend(particle_gen.generate_particles(&generator));
    }
    let mut boundaries = Vec::new();
    for boundary_gen in &config.boundary_generators {
        boundaries.extend(boundary_gen.generate_particles(&generator));
    }
    let wrapper = Wrapper::new(
        CommonParams {
            radius,
            viscosity: config.common.viscosity,
            surface_tension: config.common.surface_tension,
            adhesion: config.common.adhesion,
        },
        SimulatorParams {
            step_dt: 1.0 / config.simulator.step_per_sec,
            gravity: Vec3::from_array(config.simulator.gravity),
            min_max: Aabb::new(
                Vec3::from_array(config.simulator.min),
                Vec3::from_array(config.simulator.max),
            ),
        },
        SurfacerParams {
            density_threshold: config.surfacer.density_threshold,
            min_max: Aabb::new(
                Vec3::from_array(config.surfacer.min),
                Vec3::from_array(config.surfacer.max),
            ),
        },
        particles,
        boundaries,
        config.step_count,
        output_path.to_string(),
    );
    wrapper
}

fn setup(
    mut commands: Commands,
    mut meshes: ResMut<Assets<Mesh>>,
    mut materials: ResMut<Assets<StandardMaterial>>,
    config: Res<ConfigResource>,
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

    let wrapper = init_wrapper(&config.0, ""); // leave output path empty
    let radius = wrapper.radius;
    let particle_mesh = meshes.add(Mesh::from(Sphere::new(radius)));
    let boundary_mesh = meshes.add(Mesh::from(Sphere::new(radius / 5.0)));
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
                Mesh3d(particle_mesh.clone()),
                MeshMaterial3d(particle_material.clone()),
                Transform::from_translation(p),
                Index(i),
            ));
        });
    wrapper
        .simulator
        .boundary_positions()
        .iter()
        .for_each(|&p| {
            commands.spawn((
                Mesh3d(boundary_mesh.clone()),
                MeshMaterial3d(boundary_material.clone()),
                Transform::from_translation(p),
            ));
        });
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
