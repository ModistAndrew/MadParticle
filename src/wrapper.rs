use crate::aabb::Aabb;
use crate::simulator::{FluidParams, Simulator};
use crate::surfacer::Surfacer;
use bevy::math::Vec3;
use std::error::Error;
use std::fs::{File, create_dir_all};
use std::io::Write;
use std::path::Path;

pub struct CommonParams {
    pub radius: f32,
    pub viscosity: f32,
    pub surface_tension: f32,
    pub adhesion: f32,
}

pub struct SimulatorParams {
    pub step_dt: f32,
    pub gravity: Vec3,
    pub min_max: Aabb,
}

pub struct SurfacerParams {
    pub density_threshold: f32,
    pub min_max: Aabb,
}

pub struct Wrapper {
    pub simulator: Simulator,
    surfacer: Surfacer,
    pub radius: f32,
    tick: usize,
    sub_step_count: usize,
    path: String,
    max_step: i32,
}

impl Wrapper {
    pub fn new(
        common_params: CommonParams,
        simulator_params: SimulatorParams,
        surfacer_params: SurfacerParams,
        particles: Vec<Vec3>,
        boundaries: Vec<Vec3>,
        sub_step_count: usize,
        path: String,
        max_step: i32,
    ) -> Self {
        let fluid_params = FluidParams {
            kernel_radius: common_params.radius * 4.0,
            target_density: 1.25 / (2.0 * common_params.radius).powi(3),
            viscosity: common_params.viscosity,
            surface_tension: common_params.surface_tension,
            adhesion: common_params.adhesion,
        };
        let mut simulator = Simulator::new(
            simulator_params.step_dt,
            simulator_params.gravity,
            simulator_params.min_max,
            fluid_params,
        );
        let surfacer = Surfacer::new(
            surfacer_params.density_threshold,
            surfacer_params.min_max,
            common_params.radius,
            fluid_params,
        );
        particles
            .iter()
            .for_each(|&particle| simulator.add_particle(particle));
        boundaries
            .iter()
            .for_each(|&boundary| simulator.add_boundary(boundary));
        simulator.init_boundary();
        Self {
            simulator,
            surfacer,
            radius: common_params.radius,
            tick: 0,
            sub_step_count,
            path,
            max_step
        }
    }

    pub fn step(&mut self) {
        if self.tick % self.sub_step_count == 0 && !self.path.is_empty() {
            self.surface();
        }
        self.simulator.step();
        self.tick += 1;
        if self.max_step > 0 && self.tick > self.max_step as usize {
            println!("Simulation completed at tick {}", self.tick);
            std::process::exit(0);
        }
    }

    fn surface(&mut self) {
        let mesh = self.surfacer.surface(self.simulator.positions());
        if let Ok(_) = self.export_mesh_to_obj(mesh) {
            println!("Exported mesh to OBJ at tick {}", self.tick);
        } else {
            eprintln!("Failed to export mesh to OBJ at tick {}", self.tick);
        }
    }

    fn export_mesh_to_obj(
        &self,
        mesh: (Vec<Vec3>, Vec<Vec3>, Vec<u32>),
    ) -> Result<(), Box<dyn Error>> {
        let (vertices, normals, indices) = mesh;
        let path = format!("{}_{}.obj", self.path, self.tick / self.sub_step_count);
        create_dir_all(Path::new(&path).parent().unwrap())?;
        let mut file = File::create(path)?;
        for v in vertices {
            writeln!(file, "v {} {} {}", v[0], v[1], v[2])?;
        }
        for n in normals {
            writeln!(file, "vn {} {} {}", n[0], n[1], n[2])?;
        }
        for face in indices.chunks(3) {
            if face.len() == 3 {
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
}
