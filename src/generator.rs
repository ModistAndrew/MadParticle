use bevy::math::Vec3;

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
}
