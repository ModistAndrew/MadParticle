mod generator;
mod simulator;
mod surfacing;

use crate::generator::Generator;
use crate::simulator::Simulator;
use bevy::prelude::*;
use bevy::render::mesh::{Indices, PrimitiveTopology, VertexAttributeValues};
use std::fs::File;
use std::io::Write;

#[derive(Component)]
struct Index(usize);
#[derive(Component)]
struct SimulatorComponent(Simulator);
#[derive(Resource)]
struct SurfaceMeshHandle(Handle<Mesh>);
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
    let params = simulator::FluidParams {
            kernel_radius: radius * 4.0,
            target_density: 1.0 / (2.0 * radius).powi(3),
            viscosity: 0.5,
            surface_tension: 5.0,
            adhesion: 5.0,
    };

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
        params
    );

    let particles = Generator::from_csv("assets/output.csv")
        .unwrap_or_else(|e| panic!("Failed to load particles: {}", e));
    particles.iter().enumerate().for_each(|(i, &p)| {
        simulator.add_particle(p);
        // commands.spawn((
        //     Mesh3d(sphere_mesh.clone()),
        //     MeshMaterial3d(particle_material.clone()),
        //     Transform::from_translation(p),
        //     Index(i),
        // ));
    });

    let boundaries = generator.aabb(Vec3::new(-5.0, -2.0, -5.0), Vec3::new(5.0, -2.0, 5.0));
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

    let surfacing_material = materials.add(StandardMaterial {
        base_color: Color::srgb(0.2, 0.5, 0.8),
        perceptual_roughness: 0.5,
        metallic: 0.5,
        reflectance: 1.0,
        ..default()
    });
    let mut mesh = Mesh::new(PrimitiveTopology::TriangleList, default());
    let (position, normal, triangle_index) = surfacing::surfacing(&Vec::new(), params);
    mesh.insert_attribute(Mesh::ATTRIBUTE_POSITION, position);
    mesh.insert_attribute(Mesh::ATTRIBUTE_NORMAL, normal);
    mesh.insert_indices(Indices::U32(triangle_index));

    let mesh_handle = meshes.add(mesh);
    commands.spawn((
        Mesh3d(mesh_handle.clone()),
        MeshMaterial3d(surfacing_material),
        Transform::from_translation(Vec3::ZERO),
        GlobalTransform::default(),
    ));
    commands.insert_resource(SurfaceMeshHandle(mesh_handle));

}

static mut CNT: usize = 0;

fn pbd_step(
    mut particles: Query<(&Index, &mut Transform)>,
    mut simulator: Single<&mut SimulatorComponent>,
    mut meshes: ResMut<Assets<Mesh>>,
    mesh_handle: Res<SurfaceMeshHandle>,
) {
    simulator.0.step();
    for (id, mut transform) in particles.iter_mut() {
        transform.translation = simulator.0.get_particle_position(id.0);
    }

    let radius = 0.05;
    let params = simulator::FluidParams {
            kernel_radius: radius * 4.0,
            target_density: 1.0 / (2.0 * radius).powi(3),
            viscosity: 0.5,
            surface_tension: 5.0,
            adhesion: 5.0,
    };
    let (position, normal, triangle_index) = surfacing::surfacing(&(simulator.0.particles), params);
    if let Some(mesh) = meshes.get_mut(&mesh_handle.0) {
        mesh.insert_attribute(Mesh::ATTRIBUTE_POSITION, position);
        mesh.insert_attribute(Mesh::ATTRIBUTE_NORMAL, normal);
        mesh.insert_indices(Indices::U32(triangle_index));
        let cnt = unsafe { CNT };
        let mut path = "output/exported_mesh".to_string();
        path += &cnt.to_string();
        path += ".obj";
        unsafe { CNT += 1; }
        if let Err(e) = export_mesh_to_obj(mesh, &path) {
            eprintln!("Failed to export mesh: {}", e);
        } else {
            println!("Mesh exported to exported_mesh.obj");
        }
    }
}

fn export_mesh_to_obj(
    mesh: &Mesh,
    path: &str,
) -> Result<(), Box<dyn std::error::Error>> {
    let vertices = match mesh.attribute(Mesh::ATTRIBUTE_POSITION) {
        Some(VertexAttributeValues::Float32x3(positions)) => positions,
        _ => return Err("Mesh missing positions".into()),
    };

    let normals = match mesh.attribute(Mesh::ATTRIBUTE_NORMAL) {
        Some(VertexAttributeValues::Float32x3(normals)) => normals,
        _ => return Err("Mesh missing normals".into()),
    };

    let indices = match mesh.indices() {
        Some(Indices::U32(indices)) => indices,
        _ => return Err("Mesh missing indices or wrong format".into()),
    };

    std::fs::create_dir_all(std::path::Path::new(path).parent().unwrap())?;
    let mut file = File::create(path)?;

    for v in vertices {
        writeln!(file, "v {} {} {}", v[0], v[1], v[2])?;
    }
    for n in normals {
        writeln!(file, "vn {} {} {}", n[0], n[1], n[2])?;
    }
    for face in indices.chunks(3) {
        if face.len() == 3 {
            // Use the same index for vertex and normal
            writeln!(
                file,
                "f {0}//{0} {1}//{1} {2}//{2}",
                face[0] + 1,
                face[1] + 1,
                face[2] + 1
            )?;
        }
    }
    Ok(())
}