import sys
import struct
import numpy as np


def add_colors_to_ply(input_file, output_file):
    with open(input_file, "rb") as f:
        line = f.readline().decode("ascii").strip()
        if line != "ply":
            print("Not a valid PLY file")
            return

        lines = [line]
        vertex_count = 0
        face_count = 0
        in_header = True
        properties = []

        while in_header:
            line = f.readline().decode("ascii").strip()
            lines.append(line)

            if line.startswith("element vertex"):
                vertex_count = int(line.split()[2])
            elif line.startswith("element face"):
                face_count = int(line.split()[2])
            elif line.startswith("property") and "vertex" in " ".join(lines[-5:]):
                properties.append(line)
            elif line == "end_header":
                in_header = False
                data_start = f.tell()

        print(f"Found {vertex_count} vertices, {face_count} faces")

        f.seek(data_start)
        vertex_data = f.read(vertex_count * 12)

        vertices = []
        for i in range(vertex_count):
            x, y, z = struct.unpack("fff", vertex_data[i * 12 : (i + 1) * 12])
            vertices.append([x, y, z])

        vertices = np.array(vertices)

        y_min, y_max = vertices[:, 1].min(), vertices[:, 1].max()
        colors = []
        for v in vertices:
            t = (v[1] - y_min) / (y_max - y_min) if y_max != y_min else 0.5
            r = 0.3 + 0.4 * t
            g = 0.4 + 0.3 * (1 - t)
            b = 0.5 + 0.3 * np.sin(t * 3.14159)
            colors.append([r, g, b])

        face_data = f.read()

    with open(output_file, "w") as f:
        f.write("ply\n")
        f.write("format ascii 1.0\n")
        f.write(f"element vertex {vertex_count}\n")
        f.write("property float x\n")
        f.write("property float y\n")
        f.write("property float z\n")
        f.write("property float red\n")
        f.write("property float green\n")
        f.write("property float blue\n")
        f.write(f"element face {face_count}\n")
        f.write("property list uchar uint vertex_indices\n")
        f.write("end_header\n")

        for i in range(vertex_count):
            v = vertices[i]
            c = colors[i]
            f.write(f"{v[0]} {v[1]} {v[2]} {c[0]} {c[1]} {c[2]}\n")

        faces_parsed = []
        offset = 0
        while offset < len(face_data):
            n_verts = face_data[offset]
            offset += 1
            indices = struct.unpack("III", face_data[offset : offset + 12])
            offset += 12
            faces_parsed.append((n_verts, indices))

        for n, indices in faces_parsed:
            f.write(f"{n} {indices[0]} {indices[1]} {indices[2]}\n")

    print(f"Wrote {output_file} with vertex colors")


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python add_colors.py input.ply output.ply")
        sys.exit(1)

    add_colors_to_ply(sys.argv[1], sys.argv[2])
