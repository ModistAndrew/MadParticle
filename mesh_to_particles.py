import csv
import math
import numpy as np
import trimesh
import sys


def obj_to_particle_boundaries_sdf(obj_file_path, output_csv_path, particle_radius):
    mesh = trimesh.load_mesh(obj_file_path)
    if not isinstance(mesh, trimesh.Trimesh):
        raise ValueError("failed to load mesh")

    grid_size = 2 * particle_radius
    max_dis = math.sqrt(3) * grid_size
    bounds = mesh.bounds
    min_bound = bounds[0] - grid_size - max_dis
    max_bound = bounds[1] + grid_size + max_dis
    print("min bound: ", min_bound)
    print("max bound: ", max_bound)

    x = np.arange(min_bound[0], max_bound[0], grid_size)
    y = np.arange(min_bound[1], max_bound[1], grid_size)
    z = np.arange(min_bound[2], max_bound[2], grid_size)
    xx, yy, zz = np.meshgrid(x, y, z, indexing='ij')
    sample_points = np.vstack((xx.ravel(), yy.ravel(), zz.ravel())).T

    closest, distances, _ = mesh.nearest.on_surface(sample_points)
    surface_mask = distances <= particle_radius
    surface_points = sample_points[surface_mask]
    print("extracted {} surface points".format(len(surface_points)))

    with open(output_csv_path, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerows(surface_points)


def obj_to_particle_insides_sdf(obj_file_path, output_csv_path, particle_radius):
    mesh = trimesh.load_mesh(obj_file_path)
    if not isinstance(mesh, trimesh.Trimesh):
        raise ValueError("failed to load mesh")

    grid_size = 2 * particle_radius
    bounds = mesh.bounds
    min_bound = bounds[0] - grid_size
    max_bound = bounds[1] + grid_size
    print("min bound: ", min_bound)
    print("max bound: ", max_bound)

    x = np.arange(min_bound[0], max_bound[0], grid_size)
    y = np.arange(min_bound[1], max_bound[1], grid_size)
    z = np.arange(min_bound[2], max_bound[2], grid_size)
    xx, yy, zz = np.meshgrid(x, y, z, indexing='ij')
    sample_points = np.vstack((xx.ravel(), yy.ravel(), zz.ravel())).T

    inside_mask = mesh.contains(sample_points)
    inside_points = sample_points[inside_mask]
    print("extracted {} inside points".format(len(inside_points)))

    with open(output_csv_path, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerows(inside_points)


if __name__ == "__main__":
    argv = sys.argv
    argv = argv[argv.index("--") + 1:] if "--" in argv else []

    name = argv[0]
    b = argv[1] if len(argv) > 1 else ""
    r = argv[2] if len(argv) > 2 else "0.05"
    r = float(r)
    obj_file_path = f"assets/{name}.obj"
    output_csv_path = f"assets/{name}.csv"
    if b == "b":
        obj_to_particle_boundaries_sdf(
            obj_file_path,
            output_csv_path,
            particle_radius=0.05
        )
    else:
        obj_to_particle_insides_sdf(
            obj_file_path,
            output_csv_path,
            particle_radius=0.05
        )
