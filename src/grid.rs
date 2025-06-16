use bevy::prelude::Vec3;
use std::collections::HashMap;

pub struct Grid {
    radius: f32,
    positions: Vec<Vec3>,
    cell_size: f32,
    hash: HashMap<(i32, i32, i32), Vec<usize>>,
}

impl Grid {
    pub fn new(radius: f32, positions: Vec<Vec3>) -> Self {
        let cell_size = radius;
        let mut hash: HashMap<(i32, i32, i32), Vec<usize>> = HashMap::new();
        for (i, &position) in positions.iter().enumerate() {
            let cell = (
                (position.x / cell_size).floor() as i32,
                (position.y / cell_size).floor() as i32,
                (position.z / cell_size).floor() as i32,
            );
            hash.entry(cell).or_default().push(i);
        }
        Self {
            radius,
            positions,
            cell_size,
            hash,
        }
    }

    pub fn neighbors(&self, position: Vec3) -> Vec<usize> {
        let cell = (
            (position.x / self.cell_size).floor() as i32,
            (position.y / self.cell_size).floor() as i32,
            (position.z / self.cell_size).floor() as i32,
        );
        let mut neighbors = Vec::new();
        for dx in -1..=1 {
            for dy in -1..=1 {
                for dz in -1..=1 {
                    let neighbor_cell = (cell.0 + dx, cell.1 + dy, cell.2 + dz);
                    if let Some(indices) = self.hash.get(&neighbor_cell) {
                        for &neighbor_index in indices {
                            let dist = position.distance(self.positions[neighbor_index]);
                            if dist < self.radius {
                                neighbors.push(neighbor_index);
                            }
                        }
                    }
                }
            }
        }
        neighbors
    }

    pub fn neighbors_exclude(&self, position: Vec3, self_index: usize) -> Vec<usize> {
        let cell = (
            (position.x / self.cell_size).floor() as i32,
            (position.y / self.cell_size).floor() as i32,
            (position.z / self.cell_size).floor() as i32,
        );
        let mut neighbors = Vec::new();
        for dx in -1..=1 {
            for dy in -1..=1 {
                for dz in -1..=1 {
                    let neighbor_cell = (cell.0 + dx, cell.1 + dy, cell.2 + dz);
                    if let Some(indices) = self.hash.get(&neighbor_cell) {
                        for &neighbor_index in indices {
                            if neighbor_index == self_index {
                                continue;
                            }
                            let dist = position.distance(self.positions[neighbor_index]);
                            if dist < self.radius {
                                neighbors.push(neighbor_index);
                            }
                        }
                    }
                }
            }
        }
        neighbors
    }

    pub fn position(&self, index: usize) -> Vec3 {
        self.positions[index]
    }
}
