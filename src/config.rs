use crate::generator::Generator;
use bevy::math::Vec3;
use serde::Deserialize;

#[derive(Debug, Deserialize)]
pub struct Config {
    pub common: CommonConfig,
    pub simulator: SimulatorConfig,
    pub surfacer: SurfacerConfig,
    pub particle_generators: Vec<ParticleGenerator>,
    pub boundary_generators: Vec<ParticleGenerator>,
    pub step_count: usize,
}

#[derive(Debug, Deserialize)]
pub struct CommonConfig {
    pub radius: f32,
    pub viscosity: f32,
    pub surface_tension: f32,
    pub adhesion: f32,
}

#[derive(Debug, Deserialize)]
pub struct SimulatorConfig {
    pub step_per_sec: f32,
    pub gravity: [f32; 3],
    pub min: [f32; 3],
    pub max: [f32; 3],
}

#[derive(Debug, Deserialize)]
pub struct SurfacerConfig {
    pub density_threshold: f32,
    pub min: [f32; 3],
    pub max: [f32; 3],
}

#[derive(Debug, Deserialize)]
#[serde(tag = "type", rename_all = "snake_case")]
pub enum ParticleGenerator {
    Aabb { min: [f32; 3], max: [f32; 3] },
    ClosedBox { min: [f32; 3], max: [f32; 3] },
    OpenBox { min: [f32; 3], max: [f32; 3] },
    Csv { path: String },
}

impl ParticleGenerator {
    pub(crate) fn generate_particles(&self, generator: &Generator) -> Vec<Vec3> {
        match self {
            ParticleGenerator::Aabb { min, max } => {
                generator.aabb(Vec3::from_array(*min), Vec3::from_array(*max))
            }
            ParticleGenerator::ClosedBox { min, max } => {
                generator.closed_box(Vec3::from_array(*min), Vec3::from_array(*max))
            }
            ParticleGenerator::OpenBox { min, max } => {
                generator.open_box(Vec3::from_array(*min), Vec3::from_array(*max))
            }
            ParticleGenerator::Csv { path } => {
                Generator::from_csv(path).expect("Failed to read particles from CSV")
            }
        }
    }
}
