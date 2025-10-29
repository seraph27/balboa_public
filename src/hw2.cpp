#include "hw2.h"
#include "hw2_scenes.h"
#include <ranges>

#define sz(x) (int)(x).size()
template<class T> bool ckmin(T& a, const T& b) { return b < a ? a = b, 1 : 0; }
template<class T> bool ckmax(T& a, const T& b) { return b > a ? a = b, 1 : 0; }
using namespace hw2;

template<typename T>
Image3 supersample(int width, int height, const Vector3& bg_color, T f) {
    Image3 new_img(width * 4, height * 4);

    for(auto y = 0; y < new_img.height; y++) for(auto x = 0; x < new_img.width; x++) {
        Real xx = x + 0.5, yy = y + 0.5;
        new_img(x, y) = f(xx, yy);
    }

    Image3 res(width, height);
    for(auto y = 0; y < height; y++) for(auto x = 0; x < width; x++) {
        Vector3 sum{0, 0, 0};
        for(auto dy = 0; dy < 4; dy++) for(auto dx = 0; dx < 4; dx++) {
            sum += new_img(x * 4 + dx, y * 4 + dy);
        }
        res(x, y) = sum / 16.0;
    }
    return res;
}

Image3 hw_2_1(const std::vector<std::string> &params) {
    // Homework 2.1: render a single 3D triangle

    Image3 img(640 /* width */, 480 /* height */);

    Vector3 p0{0, 0, -1};
    Vector3 p1{1, 0, -1};
    Vector3 p2{0, 1, -1};
    Real s = 1; // scaling factor of the view frustrum
    Vector3 color = Vector3{1.0, 0.5, 0.5};
    Real z_near = 1e-6; // distance of the near clipping plane
    for (int i = 0; i < (int)params.size(); i++) {
        if (params[i] == "-s") {
            s = std::stof(params[++i]);
        } else if (params[i] == "-p0") {
            p0.x = std::stof(params[++i]);
            p0.y = std::stof(params[++i]);
            p0.z = std::stof(params[++i]);
        } else if (params[i] == "-p1") {
            p1.x = std::stof(params[++i]);
            p1.y = std::stof(params[++i]);
            p1.z = std::stof(params[++i]);
        } else if (params[i] == "-p2") {
            p2.x = std::stof(params[++i]);
            p2.y = std::stof(params[++i]);
            p2.z = std::stof(params[++i]);
        } else if (params[i] == "-color") {
            Real r = std::stof(params[++i]);
            Real g = std::stof(params[++i]);
            Real b = std::stof(params[++i]);
            color = Vector3{r, g, b};
        } else if (params[i] == "-znear") {
            z_near = std::stof(params[++i]);
        }
    }

    //reject image if vertex too close to the camera
    if(-p0.z < z_near || -p1.z < z_near || -p2.z < z_near) {
        for(int y = 0; y < img.height; y++) for(int x = 0; x < img.width; x++) {
            img(x, y) = Vector3{0.5, 0.5, 0.5};
        }
        return img;
    }

    Real x0_proj = p0.x / -p0.z, y0_proj = p0.y / -p0.z;
    Real x1_proj = p1.x / -p1.z, y1_proj = p1.y / -p1.z;
    Real x2_proj = p2.x / -p2.z, y2_proj = p2.y / -p2.z;

    Real aspect_ratio = (Real)img.width / img.height;
    Real x0_screen = img.width * (x0_proj + s * aspect_ratio) / (2 * s * aspect_ratio);
    Real y0_screen = img.height * (s - y0_proj) / (2 * s);
    Real x1_screen = img.width * (x1_proj + s * aspect_ratio) / (2 * s * aspect_ratio);
    Real y1_screen = img.height * (s - y1_proj) / (2 * s);
    Real x2_screen = img.width * (x2_proj + s * aspect_ratio) / (2 * s * aspect_ratio);
    Real y2_screen = img.height * (s - y2_proj) / (2 * s);

    auto bg = Vector3{0.5, 0.5, 0.5};

    auto [px0, py0] = Vector2{x0_screen, y0_screen};
    auto [px1, py1] = Vector2{x1_screen, y1_screen};
    auto [px2, py2] = Vector2{x2_screen, y2_screen};

    auto inside_triangle = [&](Real sx, Real sy) {
        auto e0 = (px1 - px0) * (sy - py0) - (py1 - py0) * (sx - px0);
        auto e1 = (px2 - px1) * (sy - py1) - (py2 - py1) * (sx - px1);
        auto e2 = (px0 - px2) * (sy - py2) - (py0 - py2) * (sx - px2);
        return (e0 >= 0 && e1 >= 0 && e2 >= 0) || (e0 <= 0 && e1 <= 0 && e2 <= 0);
    };

    return supersample(img.width, img.height, bg,
        [&](Real sx, Real sy) {
            if (inside_triangle(sx / 4.0, sy / 4.0)) return color;
            return bg;
        }
    );
}

Image3 hw_2_2(const std::vector<std::string> &params) {
    // Homework 2.2: render a triangle mesh

    Image3 img(640 /* width */, 480 /* height */);

    Real s = 1;
    Real z_near = 1e-6;
    int scene_id = 0;
    for (int i = 0; i < (int)params.size(); i++) {
        if (params[i] == "-s") {
            s = std::stof(params[++i]);
        } else if (params[i] == "-znear") {
            z_near = std::stof(params[++i]);
        } else if (params[i] == "-scene_id") {
            scene_id = std::stoi(params[++i]);
        }
    }

    TriangleMesh mesh = meshes[scene_id];
    Vector3 bg{0.5, 0.5, 0.5};
    Real aspect_ratio = (Real)img.width / img.height;

    auto project_to_screen = [&](const Vector3& p) -> Vector2 {
        Real x_proj = -p.x / p.z;
        Real y_proj = -p.y / p.z;
        Real x_screen = img.width * (x_proj + s * aspect_ratio) / (2 * s * aspect_ratio);
        Real y_screen = img.height * (s - y_proj) / (2 * s);
        return Vector2{x_screen, y_screen};
    };

    auto is_inside_triangle = [](Real sample_x, Real sample_y, Vector2 p0, Vector2 p1, Vector2 p2) -> bool {
        Real edge0 = (p1.x - p0.x) * (sample_y - p0.y) - (p1.y - p0.y) * (sample_x - p0.x);
        Real edge1 = (p2.x - p1.x) * (sample_y - p1.y) - (p2.y - p1.y) * (sample_x - p1.x);
        Real edge2 = (p0.x - p2.x) * (sample_y - p2.y) - (p0.y - p2.y) * (sample_x - p2.x);
        return (edge0 >= 0 && edge1 >= 0 && edge2 >= 0) || (edge0 <= 0 && edge1 <= 0 && edge2 <= 0);
    };

    auto calc_barycentric_coords = [](Real px, Real py, Real x0, Real y0, Real x1, Real y1, Real x2, Real y2) -> Vector3 {
        Real total_area = (x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0);
        Real area0 = (x1 - px) * (y2 - py) - (x2 - px) * (y1 - py);
        Real area1 = (px - x0) * (y2 - y0) - (x2 - x0) * (py - y0);
        Real area2 = (x1 - x0) * (py - y0) - (px - x0) * (y1 - y0);
        return Vector3{area0 / total_area, area1 / total_area, area2 / total_area};
    };

    return supersample(img.width, img.height, bg,
        [&](Real sx, Real sy) {
            sx /= 4.0;
            sy /= 4.0;

            auto z_min = -std::numeric_limits<Real>::infinity();
            auto res = bg;

            for (int i = 0; i < sz(mesh.faces); i++) {
                auto face = mesh.faces[i];
                auto [p0, p1, p2] = std::tuple{mesh.vertices[face[0]], mesh.vertices[face[1]], mesh.vertices[face[2]]};

                if (-p0.z < z_near || -p1.z < z_near || -p2.z < z_near) {
                    continue;
                }

                auto [px0, py0] = project_to_screen(p0);
                auto [px1, py1] = project_to_screen(p1);
                auto [px2, py2] = project_to_screen(p2);

                if (is_inside_triangle(sx, sy, {px0, py0}, {px1, py1}, {px2, py2})) {
                    auto bp = calc_barycentric_coords(sx, sy, px0, py0, px1, py1, px2, py2);

                    auto inv_d = bp[0] / p0.z + bp[1] / p1.z + bp[2] / p2.z;
                    auto b0 = (bp[0] / p0.z) / inv_d;
                    auto b1 = (bp[1] / p1.z) / inv_d;
                    auto b2 = (bp[2] / p2.z) / inv_d;
                    auto z = b0 * p0.z + b1 * p1.z + b2 * p2.z;
                    if (ckmax(z_min, z)) {
                        res = mesh.face_colors[i];
                    }
                }
            }

            return res;
        }
    );
}

Image3 hw_2_3(const std::vector<std::string> &params) {
    // Homework 2.3: render a triangle mesh with vertex colors

    Image3 img(640 /* width */, 480 /* height */);

    Real s = 1;
    Real z_near = 1e-6;
    int scene_id = 0;
    for (int i = 0; i < (int)params.size(); i++) {
        if (params[i] == "-s") {
            s = std::stof(params[++i]);
        } else if (params[i] == "-znear") {
            z_near = std::stof(params[++i]);
        } else if (params[i] == "-scene_id") {
            scene_id = std::stoi(params[++i]);
        }
    }

    TriangleMesh mesh = meshes[scene_id];
    Vector3 bg{0.5, 0.5, 0.5};
    Real aspect_ratio = (Real)img.width / img.height;

    auto project_to_screen = [&](const Vector3& p) -> Vector2 {
        Real x_proj = -p.x / p.z;
        Real y_proj = -p.y / p.z;
        Real x_screen = img.width * (x_proj + s * aspect_ratio) / (2 * s * aspect_ratio);
        Real y_screen = img.height * (s - y_proj) / (2 * s);
        return Vector2{x_screen, y_screen};
    };

    auto is_inside_triangle = [](Real sample_x, Real sample_y, Vector2 p0, Vector2 p1, Vector2 p2) -> bool {
        Real edge0 = (p1.x - p0.x) * (sample_y - p0.y) - (p1.y - p0.y) * (sample_x - p0.x);
        Real edge1 = (p2.x - p1.x) * (sample_y - p1.y) - (p2.y - p1.y) * (sample_x - p1.x);
        Real edge2 = (p0.x - p2.x) * (sample_y - p2.y) - (p0.y - p2.y) * (sample_x - p2.x);
        return (edge0 >= 0 && edge1 >= 0 && edge2 >= 0) || (edge0 <= 0 && edge1 <= 0 && edge2 <= 0);
    };

    auto calc_barycentric_coords = [](Real px, Real py, Real x0, Real y0, Real x1, Real y1, Real x2, Real y2) {
        auto tot = (x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0);
        auto a0 = (x1 - px) * (y2 - py) - (x2 - px) * (y1 - py);
        auto a1 = (px - x0) * (y2 - y0) - (x2 - x0) * (py - y0);
        auto a2 = (x1 - x0) * (py - y0) - (px - x0) * (y1 - y0);
        return Vector3{a0 / tot, a1 / tot, a2 / tot};
    };

    return supersample(img.width, img.height, bg,
        [&](Real sx, Real sy) {
            sx /= 4.0;
            sy /= 4.0;

            auto z_min = -std::numeric_limits<Real>::infinity();
            auto res = bg;

            for (int i = 0; i < sz(mesh.faces); i++) {
                auto face = mesh.faces[i];
                auto [p0, p1, p2] = std::tuple{mesh.vertices[face[0]], mesh.vertices[face[1]], mesh.vertices[face[2]]};

                if (-p0.z < z_near || -p1.z < z_near || -p2.z < z_near) {
                    continue;
                }

                auto [px0, py0] = project_to_screen(p0);
                auto [px1, py1] = project_to_screen(p1);
                auto [px2, py2] = project_to_screen(p2);

                if (is_inside_triangle(sx, sy, {px0, py0}, {px1, py1}, {px2, py2})) {
                    auto bp = calc_barycentric_coords(sx, sy, px0, py0, px1, py1, px2, py2);

                    auto inv_d = bp[0] / p0.z + bp[1] / p1.z + bp[2] / p2.z;
                    auto b0 = (bp[0] / p0.z) / inv_d;
                    auto b1 = (bp[1] / p1.z) / inv_d;
                    auto b2 = (bp[2] / p2.z) / inv_d;

                    auto z = b0 * p0.z + b1 * p1.z + b2 * p2.z;

                    if (ckmax(z_min, z)) {
                        auto c0 = mesh.vertex_colors[face[0]];
                        auto c1 = mesh.vertex_colors[face[1]];
                        auto c2 = mesh.vertex_colors[face[2]];
                        res = b0 * c0 + b1 * c1 + b2 * c2;
                    }
                }
            }

            return res;
        }
    );
}

Image3 hw_2_4(const std::vector<std::string> &params) {
     if (params.size() == 0) return Image3(0, 0);

    Scene scene = parse_scene(params[0]);
    std::cout << scene << std::endl;

    Image3 img(scene.camera.resolution.x, scene.camera.resolution.y);

    Real s = scene.camera.s;
    Real z_near = scene.camera.z_near;
    Vector3 bg = scene.background;
    Real aspect_ratio = (Real)img.width / img.height;
    Matrix4x4 view_matrix = inverse(scene.camera.cam_to_world);

    auto project_to_screen = [&](const Vector3& p) -> Vector2 {
        Real x_proj = -p.x / p.z;
        Real y_proj = -p.y / p.z;
        Real x_screen = img.width * (x_proj + s * aspect_ratio) / (2 * s * aspect_ratio);
        Real y_screen = img.height * (s - y_proj) / (2 * s);
        return Vector2{x_screen, y_screen};
    };

    auto is_inside_triangle = [](Real sample_x, Real sample_y, Vector2 p0, Vector2 p1, Vector2 p2) -> bool {
        Real edge0 = (p1.x - p0.x) * (sample_y - p0.y) - (p1.y - p0.y) * (sample_x - p0.x);
        Real edge1 = (p2.x - p1.x) * (sample_y - p1.y) - (p2.y - p1.y) * (sample_x - p1.x);
        Real edge2 = (p0.x - p2.x) * (sample_y - p2.y) - (p0.y - p2.y) * (sample_x - p2.x);
        return (edge0 >= 0 && edge1 >= 0 && edge2 >= 0) || (edge0 <= 0 && edge1 <= 0 && edge2 <= 0);
        };

    auto calc_barycentric_coords = [](Real px, Real py, Real x0, Real y0, Real x1, Real y1, Real x2, Real y2) {
        auto tot = (x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0);
        auto a0 = (x1 - px) * (y2 - py) - (x2 - px) * (y1 - py);
        auto a1 = (px - x0) * (y2 - y0) - (x2 - x0) * (py - y0);
        auto a2 = (x1 - x0) * (py - y0) - (px - x0) * (y1 - y0);
        return Vector3{a0 / tot, a1 / tot, a2 / tot};
    };

    return supersample(img.width, img.height, bg,
        [&](Real sx, Real sy) {
            sx /= 4.0;
            sy /= 4.0;

            auto z_min = -std::numeric_limits<Real>::infinity();
            auto res = bg;

            for (int mi = 0; mi < sz(scene.meshes); mi++) {
                const auto& mesh = scene.meshes[mi];
                auto xform = view_matrix * mesh.model_matrix;

                for (int i = 0; i < sz(mesh.faces); i++) {
                    auto face = mesh.faces[i];
                    auto [v0, v1, v2] = std::tuple{mesh.vertices[face[0]], mesh.vertices[face[1]], mesh.vertices[face[2]]};
                    auto p0_h = xform * Vector4(v0.x, v0.y, v0.z, Real(1));
                    auto p1_h = xform * Vector4(v1.x, v1.y, v1.z, Real(1));
                    auto p2_h = xform * Vector4(v2.x, v2.y, v2.z, Real(1));
                    auto p0 = Vector3{p0_h.x, p0_h.y, p0_h.z};
                    auto p1 = Vector3{p1_h.x, p1_h.y, p1_h.z};
                    auto p2 = Vector3{p2_h.x, p2_h.y, p2_h.z};

                    if (-p0.z < z_near || -p1.z < z_near || -p2.z < z_near) continue;

                    auto [px0, py0] = project_to_screen(p0);
                    auto [px1, py1] = project_to_screen(p1);
                    auto [px2, py2] = project_to_screen(p2);

                    if (is_inside_triangle(sx, sy, {px0, py0}, {px1, py1}, {px2, py2})) {
                        auto bp = calc_barycentric_coords(sx, sy, px0, py0, px1, py1, px2, py2);

                        auto inv_d = bp[0] / p0.z + bp[1] / p1.z + bp[2] / p2.z;
                        auto b0 = (bp[0] / p0.z) / inv_d;
                        auto b1 = (bp[1] / p1.z) / inv_d;
                        auto b2 = (bp[2] / p2.z) / inv_d;

                        auto z = b0 * p0.z + b1 * p1.z + b2 * p2.z;

                        if (ckmax(z_min, z)) {
                            auto c0 = mesh.vertex_colors[face[0]];
                            auto c1 = mesh.vertex_colors[face[1]];
                            auto c2 = mesh.vertex_colors[face[2]];
                            res = b0 * c0 + b1 * c1 + b2 * c2;
                        }
                    }
                }
            }

            return res;
        }
    );
}
