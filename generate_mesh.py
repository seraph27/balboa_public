import sys
import numpy as np
import math


def generate_torus_knot(p=2, q=3, R=1.0, r=0.3, segments=60, tube_segments=12):
    vertices = []
    colors = []
    faces = []

    for i in range(segments):
        t = 2 * math.pi * i / segments

        x_center = R * math.cos(p * t)
        y_center = R * math.sin(p * t)
        z_center = r * math.cos(q * t)

        tx = -R * p * math.sin(p * t)
        ty = R * p * math.cos(p * t)
        tz = -r * q * math.sin(q * t)
        length = math.sqrt(tx * tx + ty * ty + tz * tz)
        tx, ty, tz = tx / length, ty / length, tz / length

        nx = -math.cos(p * t)
        ny = -math.sin(p * t)

        for j in range(tube_segments):
            theta = 2 * math.pi * j / tube_segments

            bx = nx * math.cos(theta) + tz * math.sin(theta)
            by = ny * math.cos(theta)
            bz = -nx * math.sin(theta) + tz * math.cos(theta)

            tube_r = 0.1
            x = x_center + tube_r * bx
            y = y_center + tube_r * by
            z = z_center + tube_r * bz

            vertices.append([x, y, z])

            color_t = i / segments
            r = 0.3 + 0.5 * math.sin(color_t * math.pi * 2)
            g = 0.4 + 0.4 * math.cos(color_t * math.pi * 3)
            b = 0.5 + 0.3 * math.sin(color_t * math.pi * 4)
            colors.append([r, g, b])

    for i in range(segments):
        for j in range(tube_segments):
            v0 = i * tube_segments + j
            v1 = i * tube_segments + (j + 1) % tube_segments
            v2 = ((i + 1) % segments) * tube_segments + j
            v3 = ((i + 1) % segments) * tube_segments + (j + 1) % tube_segments

            faces.append([v0, v1, v2])
            faces.append([v2, v1, v3])

    return vertices, colors, faces


def generate_lattice(size=5, spacing=0.5, thickness=0.05):
    vertices = []
    colors = []
    faces = []

    def add_cube(cx, cy, cz, size_x, size_y, size_z, color):
        idx = len(vertices)

        corners = [
            [cx - size_x, cy - size_y, cz - size_z],
            [cx + size_x, cy - size_y, cz - size_z],
            [cx + size_x, cy + size_y, cz - size_z],
            [cx - size_x, cy + size_y, cz - size_z],
            [cx - size_x, cy - size_y, cz + size_z],
            [cx + size_x, cy - size_y, cz + size_z],
            [cx + size_x, cy + size_y, cz + size_z],
            [cx - size_x, cy + size_y, cz + size_z],
        ]

        for corner in corners:
            vertices.append(corner)
            colors.append(color)

        cube_faces = [
            [0, 1, 2],
            [0, 2, 3],
            [4, 6, 5],
            [4, 7, 6],
            [0, 5, 1],
            [0, 4, 5],
            [2, 6, 7],
            [2, 7, 3],
            [0, 3, 7],
            [0, 7, 4],
            [1, 5, 6],
            [1, 6, 2],
        ]

        for face in cube_faces:
            faces.append([idx + face[0], idx + face[1], idx + face[2]])

    half = size * spacing / 2

    for i in range(size + 1):
        for j in range(size + 1):
            x = i * spacing - half
            y = j * spacing - half

            t = (i + j) / (2 * size)
            r = 0.4 + 0.3 * t
            g = 0.3 + 0.4 * (1 - t)
            b = 0.5 + 0.2 * math.sin(t * math.pi)

            add_cube(x, y, -half, thickness, thickness, half, [r, g, b])
            add_cube(x, y, half, thickness, thickness, half, [r, g, b])
            add_cube(x, -half, y, thickness, half, thickness, [r, g, b])
            add_cube(x, half, y, thickness, half, thickness, [r, g, b])
            add_cube(-half, x, y, half, thickness, thickness, [r, g, b])
            add_cube(half, x, y, half, thickness, thickness, [r, g, b])

    return vertices, colors, faces


def generate_sphere_with_holes(radius=1.0, subdivisions=20):
    vertices = []
    colors = []
    faces = []

    for i in range(subdivisions):
        theta1 = math.pi * i / subdivisions
        theta2 = math.pi * (i + 1) / subdivisions

        for j in range(subdivisions * 2):
            phi1 = 2 * math.pi * j / (subdivisions * 2)
            phi2 = 2 * math.pi * (j + 1) / (subdivisions * 2)

            if (i + j) % 3 == 0:
                continue

            v0 = [
                radius * math.sin(theta1) * math.cos(phi1),
                radius * math.cos(theta1),
                radius * math.sin(theta1) * math.sin(phi1),
            ]
            v1 = [
                radius * math.sin(theta1) * math.cos(phi2),
                radius * math.cos(theta1),
                radius * math.sin(theta1) * math.sin(phi2),
            ]
            v2 = [
                radius * math.sin(theta2) * math.cos(phi1),
                radius * math.cos(theta2),
                radius * math.sin(theta2) * math.sin(phi1),
            ]
            v3 = [
                radius * math.sin(theta2) * math.cos(phi2),
                radius * math.cos(theta2),
                radius * math.sin(theta2) * math.sin(phi2),
            ]

            idx = len(vertices)

            for v in [v0, v1, v2, v3]:
                vertices.append(v)
                t = (v[1] + radius) / (2 * radius)
                r = 0.5 + 0.3 * math.sin(t * math.pi * 3)
                g = 0.4 + 0.3 * math.cos(t * math.pi * 2)
                b = 0.3 + 0.4 * t
                colors.append([r, g, b])

            faces.append([idx, idx + 1, idx + 2])
            faces.append([idx + 2, idx + 1, idx + 3])

    return vertices, colors, faces


def write_ply(filename, vertices, colors, faces):
    with open(filename, "w") as f:
        f.write("ply\n")
        f.write("format ascii 1.0\n")
        f.write(f"element vertex {len(vertices)}\n")
        f.write("property float x\n")
        f.write("property float y\n")
        f.write("property float z\n")
        f.write("property float red\n")
        f.write("property float green\n")
        f.write("property float blue\n")
        f.write(f"element face {len(faces)}\n")
        f.write("property list uchar uint vertex_indices\n")
        f.write("end_header\n")

        for i, v in enumerate(vertices):
            c = colors[i]
            f.write(f"{v[0]} {v[1]} {v[2]} {c[0]} {c[1]} {c[2]}\n")

        for face in faces:
            f.write(f"3 {face[0]} {face[1]} {face[2]}\n")

    print(f"Wrote {filename}: {len(vertices)} vertices, {len(faces)} faces")


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python generate_mesh.py [torus|lattice|sphere] output.ply")
        sys.exit(1)

    pattern = sys.argv[1]
    output = sys.argv[2]

    if pattern == "torus":
        vertices, colors, faces = generate_torus_knot()
    elif pattern == "lattice":
        vertices, colors, faces = generate_lattice(size=4, spacing=0.6)
    elif pattern == "sphere":
        vertices, colors, faces = generate_sphere_with_holes()
    else:
        print(f"Unknown pattern: {pattern}. Use torus, lattice, or sphere")
        sys.exit(1)

    write_ply(output, vertices, colors, faces)
