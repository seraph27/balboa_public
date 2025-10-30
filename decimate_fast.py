import sys
import pymeshlab
import numpy as np


def decimate_and_color(input_file, output_file, target_faces):
    print(f"Loading {input_file}...")
    ms = pymeshlab.MeshSet()
    ms.load_new_mesh(input_file)

    m = ms.current_mesh()
    print(f"Original: {m.vertex_number()} vertices, {m.face_number()} faces")

    print(f"Decimating to {target_faces} faces...")
    ms.meshing_decimation_quadric_edge_collapse(targetfacenum=target_faces)

    m = ms.current_mesh()
    print(f"After decimation: {m.vertex_number()} vertices, {m.face_number()} faces")

    print("Adding vertex colors based on height...")
    vertices = m.vertex_matrix()
    faces = m.face_matrix()

    y_min, y_max = vertices[:, 1].min(), vertices[:, 1].max()

    with open(output_file, "w") as f:
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

        for v in vertices:
            t = (v[1] - y_min) / (y_max - y_min) if y_max != y_min else 0.5
            r = 0.4 + 0.3 * t
            g = 0.3 + 0.4 * (1 - t)
            b = 0.4 + 0.2 * np.sin(t * 3.14159)
            f.write(f"{v[0]} {v[1]} {v[2]} {r} {g} {b}\n")

        for face in faces:
            f.write(f"3 {face[0]} {face[1]} {face[2]}\n")

    print(f"Saved to {output_file}")
    print("Done!")


if __name__ == "__main__":
    if len(sys.argv) < 4:
        print("Usage: python decimate_fast.py input.ply output.ply target_faces")
        sys.exit(1)

    decimate_and_color(sys.argv[1], sys.argv[2], int(sys.argv[3]))
