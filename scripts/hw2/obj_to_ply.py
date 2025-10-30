import sys
import pymeshlab
import numpy as np


def obj_to_ply(input_obj, output_ply, target_faces=5000):
    print(f"Loading {input_obj}...")
    ms = pymeshlab.MeshSet()
    ms.load_new_mesh(input_obj)

    m = ms.current_mesh()
    print(f"Original: {m.vertex_number()} vertices, {m.face_number()} faces")

    if m.face_number() > target_faces:
        print(f"Decimating to {target_faces} faces...")
        ms.meshing_decimation_quadric_edge_collapse(targetfacenum=target_faces)
        m = ms.current_mesh()
        print(
            f"After decimation: {m.vertex_number()} vertices, {m.face_number()} faces"
        )

    print("Adding vertex colors...")
    vertices = m.vertex_matrix()
    faces = m.face_matrix()

    y_min, y_max = vertices[:, 1].min(), vertices[:, 1].max()
    x_min, x_max = vertices[:, 0].min(), vertices[:, 0].max()
    z_min, z_max = vertices[:, 2].min(), vertices[:, 2].max()

    print(f"Writing {output_ply}...")
    with open(output_ply, "w") as f:
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
            ty = (v[1] - y_min) / (y_max - y_min) if y_max != y_min else 0.5
            tx = (v[0] - x_min) / (x_max - x_min) if x_max != x_min else 0.5
            tz = (v[2] - z_min) / (z_max - z_min) if z_max != z_min else 0.5

            r = 0.3 + 0.4 * ty + 0.2 * tx
            g = 0.3 + 0.3 * (1 - ty) + 0.2 * tz
            b = 0.4 + 0.3 * np.sin(ty * 3.14159) + 0.1 * tx

            r = max(0.0, min(1.0, r))
            g = max(0.0, min(1.0, g))
            b = max(0.0, min(1.0, b))

            f.write(f"{v[0]} {v[1]} {v[2]} {r} {g} {b}\n")

        for face in faces:
            f.write(f"3 {face[0]} {face[1]} {face[2]}\n")

    print(f"Saved {output_ply}")
    print("Done!")


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python obj_to_ply.py input.obj output.ply [target_faces]")
        sys.exit(1)

    input_obj = sys.argv[1]
    output_ply = sys.argv[2]
    target_faces = int(sys.argv[3]) if len(sys.argv) > 3 else 5000

    obj_to_ply(input_obj, output_ply, target_faces)
