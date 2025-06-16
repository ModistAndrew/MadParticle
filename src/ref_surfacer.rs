use crate::grid::Grid;
use crate::simulator::FluidParams;
use bevy::math::Vec3;
use isosurface::marching_cubes::MarchingCubes;
use isosurface::source::{HermiteSource, Source};

pub struct RefSurfacer {
    size: usize,
    fluid_params: FluidParams,
}

struct Sampler {
    grid: Grid,
    fluid_params: FluidParams,
}

impl Source for Sampler {
    fn sample(&self, x: f32, y: f32, z: f32) -> f32 {
        let position = Vec3::new(x, y, z);
        self.grid
            .neighbors(position)
            .iter()
            .map(|&p| self.fluid_params.kernel(position - self.grid.position(p)))
            .sum::<f32>()
    }
}

impl HermiteSource for Sampler {
    fn sample_normal(&self, x: f32, y: f32, z: f32) -> isosurface::math::Vec3 {
        let position = Vec3::new(x, y, z);
        let ret = self.grid
            .neighbors(position)
            .iter()
            .map(|&p| {
                self.fluid_params
                    .kernel_gradient(position - self.grid.position(p))
            })
            .sum::<Vec3>();
        isosurface::math::Vec3::new(ret.x, ret.y, ret.z)
    }
}

impl RefSurfacer {
    pub fn new(size: usize, fluid_params: FluidParams) -> Self {
        RefSurfacer { size, fluid_params }
    }

    pub fn surface(&self, positions: Vec<Vec3>) -> (Vec<Vec3>, Vec<Vec3>, Vec<u32>) {
        let current_time = std::time::Instant::now();

        let grid = Grid::new(self.fluid_params.kernel_radius, positions);
        let sampler = Sampler {
            grid,
            fluid_params: self.fluid_params,
        };
        let mut mc = MarchingCubes::new(self.size);
        let mut vertices = Vec::new();
        let mut indices = Vec::new();
        mc.extract_with_normals(&sampler, &mut vertices, &mut indices);
        let mut positions = Vec::new();
        let mut normals = Vec::new();
        for chunk in vertices.chunks(6) {
            if chunk.len() == 6 {
                positions.push(Vec3::new(chunk[0], chunk[1], chunk[2]));
                normals.push(Vec3::new(chunk[3], chunk[4], chunk[5]));
            }
        }

        let elapsed = current_time.elapsed();
        println!("Surfacing took: {:.2?}", elapsed);
        (positions, normals, indices)
    }
}
