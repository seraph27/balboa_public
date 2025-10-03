#include "hw1.h"
#include "hw1_scenes.h"
#include <ranges>

#define sz(x) (int)(x).size()
using namespace hw1;

Image3 hw_1_1(const std::vector<std::string> &params)
{
    // Homework 1.1: render a circle at the specified
    // position, with the specified radius and color.

    Image3 img(640 /* width */, 480 /* height */);

    Vector2 center = Vector2{img.width / 2 + Real(0.5), img.height / 2 + Real(0.5)};
    Real radius = 100.0;
    Vector3 color = Vector3{1.0, 0.5, 0.5};
    for (int i = 0; i < (int)params.size(); i++)
    {
        if (params[i] == "-center")
        {
            Real x = std::stof(params[++i]);
            Real y = std::stof(params[++i]);
            center = Vector2{x, y};
        }
        else if (params[i] == "-radius")
        {
            radius = std::stof(params[++i]);
        }
        else if (params[i] == "-color")
        {
            Real r = std::stof(params[++i]);
            Real g = std::stof(params[++i]);
            Real b = std::stof(params[++i]);
            color = Vector3{r, g, b};
        }
    }

    for (int y = 0; y < img.height; y++)
    {
        for (int x = 0; x < img.width; x++)
        {
            Real pixel_x = x + 0.5;
            Real pixel_y = (img.height - 1 - y) + 0.5;
            Real dx = pixel_x - center.x;
            Real dy = pixel_y - center.y;

            if (dx * dx + dy * dy <= radius * radius)
            {
                img(x, y) = color;
            }
            else
            {
                img(x, y) = Vector3{0.5, 0.5, 0.5};
            }
        }
    }
    return img;
}

Image3 hw_1_2(const std::vector<std::string> &params)
{
    // Homework 1.2: render polylines
    if (params.size() == 0)
    {
        return Image3(0, 0);
    }

    Image3 img(640 /* width */, 480 /* height */);
    std::vector<Vector2> polyline;
    // is_closed = true indicates that the last point and
    // the first point of the polyline are connected
    bool is_closed = false;
    std::optional<Vector3> fill_color;
    std::optional<Vector3> stroke_color;
    Real stroke_width = 1;
    for (int i = 0; i < (int)params.size(); i++)
    {
        if (params[i] == "-points")
        {
            while (params.size() > i + 1 &&
                   params[i + 1].length() > 0 &&
                   params[i + 1][0] != '-')
            {
                Real x = std::stof(params[++i]);
                Real y = std::stof(params[++i]);
                polyline.push_back(Vector2{x, y});
            }
        }
        else if (params[i] == "--closed")
        {
            is_closed = true;
        }
        else if (params[i] == "-fill_color")
        {
            Real r = std::stof(params[++i]);
            Real g = std::stof(params[++i]);
            Real b = std::stof(params[++i]);
            fill_color = Vector3{r, g, b};
        }
        else if (params[i] == "-stroke_color")
        {
            Real r = std::stof(params[++i]);
            Real g = std::stof(params[++i]);
            Real b = std::stof(params[++i]);
            stroke_color = Vector3{r, g, b};
        }
        else if (params[i] == "-stroke_width")
        {
            stroke_width = std::stof(params[++i]);
        }
    }

    if (fill_color && !is_closed)
    {
        std::cout << "Error: can't have a non-closed shape with fill color." << std::endl;
        return Image3(0, 0);
    }

    auto inside = [&](Vector2 point) -> bool
    {
        int winding = 0;
        for (int i = 0; i < sz(polyline); i++)
        {
            auto j = (i + 1) % sz(polyline);
            if ((polyline[i].y > point.y) != (polyline[j].y > point.y) &&
                (point.x < (polyline[j].x - polyline[i].x) * (point.y - polyline[i].y) / (polyline[j].y - polyline[i].y) + polyline[i].x))
            {
                winding += polyline[j].y > polyline[i].y ? 1 : -1;
            }
        }
        return winding != 0;
    };

    auto dist_pt_to_line = [](Vector2 v, Vector2 a, Vector2 b)
    {
        Vector2 ab = b - a;
        Vector2 av = v - a;
        Real ab_len2 = dot(ab, ab);
        if (!ab_len2)
            return length(av);
        Real t = dot(av, ab) / ab_len2;
        t = std::clamp(t, 0.0, 1.0);
        Vector2 projection = a + t * ab;
        return length(v - projection);
    };

    for (int y = 0; y < img.height; y++)
    {
        for (int x = 0; x < img.width; x++)
        {
            Real pixel_x = x + 0.5;
            Real pixel_y = (img.height - 1 - y) + 0.5;

            if (is_closed && fill_color.has_value() && inside(Vector2{pixel_x, pixel_y}))
            {
                img(x, y) = *fill_color;
            }
            else
            {
                img(x, y) = Vector3{0.5, 0.5, 0.5};
            }

            if (stroke_color.has_value())
            {
                for (int i = 0; i < (is_closed ? sz(polyline) : sz(polyline) - 1); i++)
                {
                    int j = (i + 1) % sz(polyline);
                    auto dist = dist_pt_to_line(Vector2{pixel_x, pixel_y}, polyline[i], polyline[j]);
                    if (dist <= stroke_width / 2)
                    {
                        img(x, y) = *stroke_color;
                        break;
                    }
                }
            }
        }
    }
    return img;
}

Image3 hw_1_3(const std::vector<std::string> &params)
{
    if (params.size() == 0)
    {
        return Image3(0, 0);
    }

    Scene scene = parse_scene(params[0]);
    std::cout << scene << std::endl;

    Image3 img(scene.resolution.x, scene.resolution.y);
    
    for (int y = 0; y < img.height; y++)
    {
        for (int x = 0; x < img.width; x++)
        {
            img(x, y) = scene.background;
        }
    }
    
    for (auto shape : scene.shapes | std::views::reverse) {
        if(auto *circle = std::get_if<Circle>(&shape)) {
            for (int y = 0; y < img.height; y++) {
                for (int x = 0; x < img.width; x++) {
                    Real pixel_x = x + 0.5;
                    Real pixel_y = (img.height - 1 - y) + 0.5;
                    Real dx = pixel_x - circle->center.x;
                    Real dy = pixel_y - circle->center.y;
                    Real dist = sqrt(dx * dx + dy * dy);
                    
                    if (circle->fill_color.has_value() && dist <= circle->radius) {
                        img(x, y) = *circle->fill_color;
                    }
                    
                    if (circle->stroke_color.has_value()) {
                        Real outer_radius = circle->radius + circle->stroke_width / 2;
                        Real inner_radius = circle->radius - circle->stroke_width / 2;
                        if (dist >= inner_radius && dist <= outer_radius) {
                            img(x, y) = *circle->stroke_color;
                        }
                    }
                }
            }
        } else if(auto *polyline = std::get_if<Polyline>(&shape)) {
            auto inside = [&](Vector2 point) -> bool {
                int winding = 0;
                for (int i = 0; i < sz(polyline->points); i++) {
                    auto j = (i + 1) % sz(polyline->points);
                    if ((polyline->points[i].y > point.y) != (polyline->points[j].y > point.y) &&
                        (point.x < (polyline->points[j].x - polyline->points[i].x) * (point.y - polyline->points[i].y) / 
                         (polyline->points[j].y - polyline->points[i].y) + polyline->points[i].x)) {
                        winding += polyline->points[j].y > polyline->points[i].y ? 1 : -1;
                    }
                }
                return winding != 0;
            };

            auto dist_pt_to_line = [](Vector2 v, Vector2 a, Vector2 b) {
                Vector2 ab = b - a;
                Vector2 av = v - a;
                Real ab_len2 = dot(ab, ab);
                if (!ab_len2) return length(av);
                Real t = dot(av, ab) / ab_len2;
                t = std::clamp(t, 0.0, 1.0);
                Vector2 projection = a + t * ab;
                return length(v - projection);
            };

            for (int y = 0; y < img.height; y++) {
                for (int x = 0; x < img.width; x++) {
                    Real pixel_x = x + 0.5;
                    Real pixel_y = (img.height - 1 - y) + 0.5;
                    Vector2 pixel_pos = Vector2{pixel_x, pixel_y};

                    if (polyline->is_closed && polyline->fill_color.has_value() && inside(pixel_pos)) {
                        img(x, y) = *polyline->fill_color;
                    }

                    if (polyline->stroke_color.has_value()) {
                        for (int i = 0; i < (polyline->is_closed ? sz(polyline->points) : sz(polyline->points) - 1); i++) {
                            int j = (i + 1) % sz(polyline->points);
                            auto dist = dist_pt_to_line(pixel_pos, polyline->points[i], polyline->points[j]);
                            if (dist <= polyline->stroke_width / 2) {
                                img(x, y) = *polyline->stroke_color;
                                break;
                            }
                        }
                    }
                }
            }
        }
    }
    
    return img;
}

Image3 hw_1_4(const std::vector<std::string> &params)
{
    // Homework 1.4: render transformed shapes
    if (params.size() == 0)
    {
        return Image3(0, 0);
    }

    Scene scene = parse_scene(params[0]);
    std::cout << scene << std::endl;

    Image3 img(scene.resolution.x, scene.resolution.y);

    for (int y = 0; y < img.height; y++)
    {
        for (int x = 0; x < img.width; x++)
        {
            img(x, y) = Vector3{1, 1, 1};
        }
    }
    return img;
}

Image3 hw_1_5(const std::vector<std::string> &params)
{
    // Homework 1.5: antialiasing
    if (params.size() == 0)
    {
        return Image3(0, 0);
    }

    Scene scene = parse_scene(params[0]);
    std::cout << scene << std::endl;

    Image3 img(scene.resolution.x, scene.resolution.y);

    for (int y = 0; y < img.height; y++)
    {
        for (int x = 0; x < img.width; x++)
        {
            img(x, y) = Vector3{1, 1, 1};
        }
    }
    return img;
}

Image3 hw_1_6(const std::vector<std::string> &params)
{
    // Homework 1.6: alpha blending
    if (params.size() == 0)
    {
        return Image3(0, 0);
    }

    Scene scene = parse_scene(params[0]);
    std::cout << scene << std::endl;

    Image3 img(scene.resolution.x, scene.resolution.y);

    for (int y = 0; y < img.height; y++)
    {
        for (int x = 0; x < img.width; x++)
        {
            img(x, y) = Vector3{1, 1, 1};
        }
    }
    return img;
}
