import sys
import struct
import numpy as np


def read_binary_ply(filename):
    with open(filename, "rb") as f:
        line = f.readline().decode("ascii").strip()
        if line != "ply":
            raise ValueError("Not a valid PLY file")

        format_line = f.readline().decode("ascii").strip()
        if "binary_little_endian" not in format_line:
            raise ValueError("Only binary little endian PLY supported")

        vertex_count = 0
        face_count = 0

        while True:
            line = f.readline().decode("ascii").strip()
            if line.startswith("element vertex"):
                vertex_count = int(line.split()[2])
            elif line.startswith("element face"):
                face_count = int(line.split()[2])
            elif line == "end_header":
                break

        print(f"Reading {vertex_count} vertices and {face_count} faces...")

        vertices = np.zeros((vertex_count, 3), dtype=np.float32)
        for i in range(vertex_count):
            data = f.read(12)
            vertices[i] = struct.unpack("fff", data)

        faces = []
        for i in range(face_count):
            n = struct.unpack("B", f.read(1))[0]
            if n == 3:
                indices = struct.unpack("III", f.read(12))
                faces.append(indices)
            else:
                f.read(n * 4)

        return vertices, np.array(faces, dtype=np.int32)


def compute_edge_collapse_cost(vertices, faces, vertex_idx):
    connected_faces = [f for f in faces if vertex_idx in f]
    if len(connected_faces) == 0:
        return float("inf"), None

    neighbors = set()
    for face in connected_faces:
        for v in face:
            if v != vertex_idx:
                neighbors.add(v)

    if len(neighbors) == 0:
        return float("inf"), None

    min_cost = float("inf")
    best_target = None
    pos = vertices[vertex_idx]

    for n in neighbors:
        cost = np.linalg.norm(vertices[n] - pos)
        if cost < min_cost:
            min_cost = cost
            best_target = n

    return min_cost, best_target


def decimate_mesh(vertices, faces, target_faces):
    print(f"Decimating from {len(faces)} to {target_faces} faces...")

    faces_list = [tuple(f) for f in faces]

    while len(faces_list) > target_faces:
        if len(faces_list) % 10000 == 0:
            print(f"  {len(faces_list)} faces remaining...")

        vertex_usage = {}
        for face in faces_list:
            for v in face:
                vertex_usage[v] = vertex_usage.get(v, 0) + 1

        min_cost = float("inf")
        best_vertex = None
        best_target = None

        for v in vertex_usage:
            if vertex_usage[v] < 3:
                continue
            cost, target = compute_edge_collapse_cost(vertices, faces_list, v)
            if cost < min_cost:
                min_cost = cost
                best_vertex = v
                best_target = target

        if best_vertex is None or best_target is None:
            break

        new_faces = []
        for face in faces_list:
            if best_vertex in face:
                new_face = tuple(best_target if v == best_vertex else v for v in face)
                if len(set(new_face)) == 3:
                    new_faces.append(new_face)
            else:
                new_faces.append(face)

        faces_list = new_faces

    return vertices, np.array(faces_list, dtype=np.int32)


def add_vertex_colors(vertices):
    y_min, y_max = vertices[:, 1].min(), vertices[:, 1].max()
    colors = np.zeros((len(vertices), 3), dtype=np.float32)

    for i, v in enumerate(vertices):
        t = (v[1] - y_min) / (y_max - y_min) if y_max != y_min else 0.5
        colors[i, 0] = 0.4 + 0.3 * t
        colors[i, 1] = 0.3 + 0.4 * (1 - t)
        colors[i, 2] = 0.4 + 0.2 * np.sin(t * 3.14159)

    return colors


def write_ascii_ply(filename, vertices, faces, colors):
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

        for i in range(len(vertices)):
            v = vertices[i]
            c = colors[i]
            f.write(f"{v[0]} {v[1]} {v[2]} {c[0]} {c[1]} {c[2]}\n")

        for face in faces:
            f.write(f"3 {face[0]} {face[1]} {face[2]}\n")

    print(f"Wrote {filename}")


if __name__ == "__main__":
    if len(sys.argv) < 4:
        print("Usage: python decimate_mesh.py input.ply output.ply target_faces")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]
    target_faces = int(sys.argv[3])

    vertices, faces = read_binary_ply(input_file)
    vertices, faces = decimate_mesh(vertices, faces, target_faces)
    colors = add_vertex_colors(vertices)
    write_ascii_ply(output_file, vertices, faces, colors)
    print("Done!")
