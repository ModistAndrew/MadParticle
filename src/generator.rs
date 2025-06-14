use bevy::math::Vec3;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

pub struct Generator {
    radius: f32,
}

impl Generator {
    pub fn new(radius: f32) -> Self {
        Self { radius }
    }

    pub fn aabb(&self, from: Vec3, to: Vec3) -> Vec<Vec3> {
        let mut points = Vec::new();
        let step = self.radius * 2.0;
        let x_count = ((to.x - from.x) / step).floor() as usize + 1;
        let y_count = ((to.y - from.y) / step).floor() as usize + 1;
        let z_count = ((to.z - from.z) / step).floor() as usize + 1;
        for x in 0..x_count {
            for y in 0..y_count {
                for z in 0..z_count {
                    let point = Vec3::new(
                        from.x + x as f32 * step,
                        from.y + y as f32 * step,
                        from.z + z as f32 * step,
                    );
                    points.push(point);
                }
            }
        }
        points
    }

    pub fn open_box(&self, from: Vec3, to: Vec3) -> Vec<Vec3> {
        let mut points = Vec::new();
        let step = self.radius * 2.0;
        let x_count = ((to.x - from.x) / step).floor() as usize + 1;
        let y_count = ((to.y - from.y) / step).floor() as usize + 1;
        let z_count = ((to.z - from.z) / step).floor() as usize + 1;

        for x in 0..x_count {
            for y in 0..y_count {
                for z in 0..z_count {
                    if (x == 0 || x == x_count - 1) || (y == 0) || (z == 0 || z == z_count - 1) {
                        let point = Vec3::new(
                            from.x + x as f32 * step,
                            from.y + y as f32 * step,
                            from.z + z as f32 * step,
                        );
                        points.push(point);
                    }
                }
            }
        }
        points
    }

    pub fn from_csv(path: &str) -> Result<Vec<Vec3>, String> {
        let file =
            File::open(Path::new(path)).map_err(|e| format!("Failed to open file: {}", e))?;
        let reader = BufReader::new(file);
        let mut particles = Vec::new();
        for (i, line) in reader.lines().enumerate() {
            let line = line.map_err(|e| format!("Error reading line {}: {}", i + 1, e))?;
            let parts: Vec<&str> = line.split(',').collect();
            if parts.len() < 3 {
                return Err(format!(
                    "Invalid format at line {}: expected 3 values",
                    i + 1
                ));
            }
            let x = parts[0]
                .parse::<f32>()
                .map_err(|e| format!("Error parsing x at line {}: {}", i + 1, e))?;
            let y = parts[1]
                .parse::<f32>()
                .map_err(|e| format!("Error parsing y at line {}: {}", i + 1, e))?;
            let z = parts[2]
                .parse::<f32>()
                .map_err(|e| format!("Error parsing z at line {}: {}", i + 1, e))?;
            particles.push(Vec3::new(x, y, z));
        }
        println!("Loaded {} particles from {}", particles.len(), path);
        Ok(particles)
    }
}
