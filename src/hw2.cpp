#include "hw2.h"
#include "hw2_scenes.h"
#include <ranges>

#define sz(x) (int)(x).size()
template<class T> bool ckmin(T& a, const T& b) { return b < a ? a = b, 1 : 0; }
template<class T> bool ckmax(T& a, const T& b) { return b > a ? a = b, 1 : 0; }
using namespace hw2;

struct ClippedTri {
    Vector3 v[3];
};

int clip_triangle(Vector3 p0, Vector3 p1, Vector3 p2, Real z_near, ClippedTri out[2]) {
    auto behind = [&](const Vector3& p) { return -p.z < z_near; };
    auto lerp_z = [&](const Vector3& a, const Vector3& b) {
        auto t = (z_near + a.z) / (a.z - b.z);
        return Vector3{a.x + t * (b.x - a.x), a.y + t * (b.y - a.y), -z_near};
    };
    
    int cnt = behind(p0) + behind(p1) + behind(p2);
    if (cnt == 0) {
        out[0].v[0] = p0; out[0].v[1] = p1; out[0].v[2] = p2;
        return 1;
    }
    if (cnt == 3) return 0;
    
    if (cnt == 1) {
        Vector3 front1, front2, back;
        if (behind(p0)) { back = p0; front1 = p1; front2 = p2; }
        else if (behind(p1)) { back = p1; front1 = p0; front2 = p2; }
        else { back = p2; front1 = p0; front2 = p1; }
        
        auto i1 = lerp_z(back, front1);
        auto i2 = lerp_z(back, front2);
        
        out[0].v[0] = i1; out[0].v[1] = front1; out[0].v[2] = front2;
        out[1].v[0] = i1; out[1].v[1] = front2; out[1].v[2] = i2;
        return 2;
    }
    
    Vector3 front, back1, back2;
    if (!behind(p0)) { front = p0; back1 = p1; back2 = p2; }
    else if (!behind(p1)) { front = p1; back1 = p0; back2 = p2; }
    else { front = p2; back1 = p0; back2 = p1; }
    
    auto i1 = lerp_z(front, back1);
    auto i2 = lerp_z(front, back2);
    
    out[0].v[0] = front; out[0].v[1] = i1; out[0].v[2] = i2;
    return 1;
}

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

    ClippedTri clipped[2];
    int num_tris = clip_triangle(p0, p1, p2, z_near, clipped);
    
    auto bg = Vector3{0.5, 0.5, 0.5};
    auto aspect_ratio = (Real)img.width / img.height;
    
    auto project_to_screen = [&](const Vector3& p) -> Vector2 {
        Real x_proj = -p.x / p.z;
        Real y_proj = -p.y / p.z;
        Real x_screen = img.width * (x_proj + s * aspect_ratio) / (2 * s * aspect_ratio);
        Real y_screen = img.height * (s - y_proj) / (2 * s);
        return Vector2{x_screen, y_screen};
    };

    auto is_inside_triangle = [](Real sx, Real sy, Vector2 p0, Vector2 p1, Vector2 p2) -> bool {
        auto e0 = (p1.x - p0.x) * (sy - p0.y) - (p1.y - p0.y) * (sx - p0.x);
        auto e1 = (p2.x - p1.x) * (sy - p1.y) - (p2.y - p1.y) * (sx - p1.x);
        auto e2 = (p0.x - p2.x) * (sy - p2.y) - (p0.y - p2.y) * (sx - p2.x);
        return (e0 >= 0 && e1 >= 0 && e2 >= 0) || (e0 <= 0 && e1 <= 0 && e2 <= 0);
    };

    return supersample(img.width, img.height, bg,
        [&](Real sx, Real sy) {
            sx /= 4.0;
            sy /= 4.0;
            
            for (int ti = 0; ti < num_tris; ti++) {
                auto cp0 = clipped[ti].v[0], cp1 = clipped[ti].v[1], cp2 = clipped[ti].v[2];
                
                auto [px0, py0] = project_to_screen(cp0);
                auto [px1, py1] = project_to_screen(cp1);
                auto [px2, py2] = project_to_screen(cp2);
                
                if (is_inside_triangle(sx, sy, {px0, py0}, {px1, py1}, {px2, py2})) {
                    return color;
                }
            }
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

                ClippedTri clipped[2];
                int num_tris = clip_triangle(p0, p1, p2, z_near, clipped);
                
                for (int ti = 0; ti < num_tris; ti++) {
                    auto cp0 = clipped[ti].v[0], cp1 = clipped[ti].v[1], cp2 = clipped[ti].v[2];
                    
                    auto [px0, py0] = project_to_screen(cp0);
                    auto [px1, py1] = project_to_screen(cp1);
                    auto [px2, py2] = project_to_screen(cp2);
                    
                    if (is_inside_triangle(sx, sy, {px0, py0}, {px1, py1}, {px2, py2})) {
                        auto bp = calc_barycentric_coords(sx, sy, px0, py0, px1, py1, px2, py2);
                        
                        auto inv_d = bp[0] / cp0.z + bp[1] / cp1.z + bp[2] / cp2.z;
                        auto b0 = (bp[0] / cp0.z) / inv_d;
                        auto b1 = (bp[1] / cp1.z) / inv_d;
                        auto b2 = (bp[2] / cp2.z) / inv_d;
                        auto z = b0 * cp0.z + b1 * cp1.z + b2 * cp2.z;
                        if (ckmax(z_min, z)) {
                            res = mesh.face_colors[i];
                        }
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

                ClippedTri clipped[2];
                int num_tris = clip_triangle(p0, p1, p2, z_near, clipped);
                
                for (int ti = 0; ti < num_tris; ti++) {
                    auto cp0 = clipped[ti].v[0], cp1 = clipped[ti].v[1], cp2 = clipped[ti].v[2];
                    
                    auto [px0, py0] = project_to_screen(cp0);
                    auto [px1, py1] = project_to_screen(cp1);
                    auto [px2, py2] = project_to_screen(cp2);
                    
                    if (is_inside_triangle(sx, sy, {px0, py0}, {px1, py1}, {px2, py2})) {
                        auto bp = calc_barycentric_coords(sx, sy, px0, py0, px1, py1, px2, py2);
                        
                        auto inv_d = bp[0] / cp0.z + bp[1] / cp1.z + bp[2] / cp2.z;
                        auto b0 = (bp[0] / cp0.z) / inv_d;
                        auto b1 = (bp[1] / cp1.z) / inv_d;
                        auto b2 = (bp[2] / cp2.z) / inv_d;
                        
                        auto z = b0 * cp0.z + b1 * cp1.z + b2 * cp2.z;
                        
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

                    ClippedTri clipped[2];
                    int num_tris = clip_triangle(p0, p1, p2, z_near, clipped);
                    
                    for (int ti = 0; ti < num_tris; ti++) {
                        auto cp0 = clipped[ti].v[0], cp1 = clipped[ti].v[1], cp2 = clipped[ti].v[2];

                        auto [px0, py0] = project_to_screen(cp0);
                        auto [px1, py1] = project_to_screen(cp1);
                        auto [px2, py2] = project_to_screen(cp2);
                    
                        if (is_inside_triangle(sx, sy, {px0, py0}, {px1, py1}, {px2, py2})) {
                            auto bp = calc_barycentric_coords(sx, sy, px0, py0, px1, py1, px2, py2);
                        
                            auto inv_d = bp[0] / cp0.z + bp[1] / cp1.z + bp[2] / cp2.z;
                            auto b0 = (bp[0] / cp0.z) / inv_d;
                            auto b1 = (bp[1] / cp1.z) / inv_d;
                            auto b2 = (bp[2] / cp2.z) / inv_d;
                        
                            auto z = b0 * cp0.z + b1 * cp1.z + b2 * cp2.z;
                        
                            if (ckmax(z_min, z)) {
                                auto c0 = mesh.vertex_colors[face[0]];
                                auto c1 = mesh.vertex_colors[face[1]];
                                auto c2 = mesh.vertex_colors[face[2]];
                                res = b0 * c0 + b1 * c1 + b2 * c2;
                            }
                        }
                    }
                }
            }

            return res;
        }
    );
}
