// This example is heavily based on the tutorial at https://open.gl

// OpenGL Helpers to reduce the clutter
#include "Helpers.h"
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <cmath>
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#ifdef __APPLE__
#define GL_SILENCE_DEPRECATION
// GLFW is necessary to handle the OpenGL context
#include <GLFW/glfw3.h>
#else
// GLFW is necessary to handle the OpenGL context
#include <GLFW/glfw3.h>
#endif

// Timer
#include <chrono>

// OpenGL Mathematics Library
#include <glm/glm.hpp> // glm::vec3
#include <glm/vec3.hpp> // glm::vec3
#include <glm/vec4.hpp> // glm::vec4
#include <glm/mat4x4.hpp> // glm::mat4
#include <glm/gtc/matrix_transform.hpp> // glm::translate, glm::rotate, glm::scale, glm::perspective
#include <glm/gtc/type_ptr.hpp> // glm::value_ptr

// VertexBufferObject wrapper
VertexBufferObject VBO;
VertexBufferObject VBO_M;
VertexBufferObject VBO_N_f;
VertexBufferObject VBO_N_p;

// Contains the vertex positions and its mid in triangles
std::vector<glm::vec3> V;
std::vector<glm::vec3> M;

// Contains the vertex normal
std::vector<glm::vec3> N_f;
std::vector<glm::vec3> N_p;

// Contains base color for each model
std::vector<glm::vec3> base_color;

// Contains number of triangles for each model
std::vector<int> numbers;

// Contains model matrices
std::vector<glm::mat4> model;

// Save model render mode
std::vector<unsigned short> render_mode;

// Contains camera
glm::vec3 eye = glm::vec3(4.0f, 3.0f, 3.0f);
glm::mat4 camera = glm::lookAt(eye, glm::vec3(0, 0, 0), glm::vec3(0, 1, 0));

// Contains projection
bool is_perspective = true;
float window_width = 640.0f;
float window_height = 480.0f;
glm::mat4 projection = glm::perspective(glm::radians(45.0f), window_width/window_height, 0.1f, 100.0f);

// Save light source
glm::vec3 light_source = glm::vec3(2.0f);

// light projection
float near_plane = 1.0f, far_plane = 7.5f;
glm::mat4 light_projection = glm::ortho(-10.0f, 10.0f, -10.0f, 10.0f, near_plane, far_plane);

// shadow mode
int shadow_mode = 1;

static const float unit_cube[] = {
        -0.5f,-0.5f,-0.5f, // triangle 1 : begin
        -0.5f,-0.5f, 0.5f,
        -0.5f, 0.5f, 0.5f, // triangle 1 : end
        0.5f, 0.5f,-0.5f, // triangle 2 : begin
        -0.5f,-0.5f,-0.5f,
        -0.5f, 0.5f,-0.5f, // triangle 2 : end
        0.5f,-0.5f, 0.5f,
        -0.5f,-0.5f,-0.5f,
        0.5f,-0.5f,-0.5f,
        0.5f, 0.5f,-0.5f,
        0.5f,-0.5f,-0.5f,
        -0.5f,-0.5f,-0.5f,
        -0.5f,-0.5f,-0.5f,
        -0.5f, 0.5f, 0.5f,
        -0.5f, 0.5f,-0.5f,
        0.5f,-0.5f, 0.5f,
        -0.5f,-0.5f, 0.5f,
        -0.5f,-0.5f,-0.5f,
        -0.5f, 0.5f, 0.5f,
        -0.5f,-0.5f, 0.5f,
        0.5f,-0.5f, 0.5f,
        0.5f, 0.5f, 0.5f,
        0.5f,-0.5f,-0.5f,
        0.5f, 0.5f,-0.5f,
        0.5f,-0.5f,-0.5f,
        0.5f, 0.5f, 0.5f,
        0.5f,-0.5f, 0.5f,
        0.5f, 0.5f, 0.5f,
        0.5f, 0.5f,-0.5f,
        -0.5f, 0.5f,-0.5f,
        0.5f, 0.5f, 0.5f,
        -0.5f, 0.5f,-0.5f,
        -0.5f, 0.5f, 0.5f,
        0.5f, 0.5f, 0.5f,
        -0.5f, 0.5f, 0.5f,
        0.5f,-0.5f, 0.5f
};

// Helper functions
bool in_triangle(glm::vec3 v1, glm::vec3 v2, glm::vec3 v3, glm::vec3 inter) {
    float e12 = glm::distance(v1, v2), e13 = glm::distance(v1, v3), e23 = glm::distance(v2, v3);
    float ei1 = glm::distance(v1, inter), ei2 = glm::distance(v2, inter), ei3 = glm::distance(v3, inter);

    float si12 = (ei1 + ei2 + e12)/2, si13 = (ei1 + ei3 + e13)/2, si23 = (ei2 + ei3 + e23)/2, s123 = (e12 + e13 + e23)/2;

    float ai12 = sqrt(si12*(si12 - ei1)*(si12 - ei2)*(si12 - e12));
    float ai13 = sqrt(si13*(si13 - ei1)*(si13 - ei3)*(si13 - e13));
    float ai23 = sqrt(si23*(si23 - ei2)*(si23 - ei3)*(si23 - e23));
    float a123 = sqrt(s123*(s123 - e12)*(s123 - e13)*(s123 - e23));

    return abs(ai12 + ai13 + ai23 - a123) < 0.000001f;
}

struct vecHash {
    std::size_t operator()(const glm::vec3& k) const
    {
        return std::hash<float>()(k[0]) ^
               (std::hash<float>()(k[1]) << 1) ^
               (std::hash<float>()(k[2]) << 2);
    }
};

void load_model(const std::string& path) {
    std::ifstream in = std::ifstream(path, std::fstream::in);
    if(!in.is_open()) throw std::runtime_error("no such model");

    std::vector<glm::vec3> vertices;
    std::string line;
    std::string word;
    glm::vec3 center(1.0f, 1.0f, 1.0f);
    std::unordered_map<glm::vec3, std::vector<glm::vec3>, vecHash> normal_map;
    int vertex_num, triangle_num;
    float x, y, z;
    float x_max = FLT_MIN, x_min = FLT_MAX, y_max = FLT_MIN, y_min = FLT_MAX, z_max = FLT_MIN, z_min = FLT_MAX;

    // skip head
    getline(in, line);

    // get numbers
    getline(in, line);
    std::stringstream s(line);
    s >> word;
    vertex_num = std::stoi(word);
    s >> word;
    triangle_num = std::stoi(word);

    // get all vertices and compute max and min values for x y z
    for(int i=0;i<vertex_num;i++) {
        getline(in, line);
        std::stringstream ss(line);
        ss >> word;
        x = std::stof(word);
        ss >> word;
        y = std::stof(word);
        ss >> word;
        z = std::stof(word);

        vertices.emplace_back(x, y, z);
        center[0] += x;
        center[1] += y;
        center[2] += z;
        x_max = x > x_max? x: x_max;
        x_min = x < x_min? x: x_min;
        y_max = y > y_max? y: y_max;
        y_min = y < y_min? y: y_min;
        z_max = z > z_max? z: z_max;
        z_min = z < z_min? z: z_min;
    }
    center[0] = center[0]/vertex_num;
    center[1] = center[1]/vertex_num;
    center[2] = center[2]/vertex_num;

    // get the transform and scale of the original data
    float x_diff = x_max-center[0]>center[0]-x_min? x_max-center[0]: center[0]-x_min;
    float y_diff = y_max-center[1]>center[1]-y_min? y_max-center[1]: center[1]-y_min;
    float z_diff = z_max-center[2]>center[2]-z_min? z_max-center[2]: center[2]-z_min;
    float scale = x_diff > y_diff? (x_diff>z_diff? 1/x_diff: 1/z_diff): (y_diff>z_diff? 1/y_diff: 1/z_diff);
    glm::mat4 to_origin = glm::translate(glm::mat4(1.0f), - center);
    glm::mat4 scaling = glm::scale(glm::mat4(1.0f), glm::vec3(scale/2));

    // append replaced vertices and color
    for(int i=0;i<triangle_num;i++) {
        getline(in, line);
        std::stringstream ss(line);
        ss >> word;
        glm::vec3 v1, v2, v3;

        // replace a triangle
        ss >> word;
        v1 = scaling * to_origin * glm::vec4(vertices[std::stoi(word)], 1.0f);
        ss >> word;
        v2 = scaling * to_origin * glm::vec4(vertices[std::stoi(word)], 1.0f);
        ss >> word;
        v3 = scaling * to_origin * glm::vec4(vertices[std::stoi(word)], 1.0f);

        // append replaced vertices
        V.emplace_back(v1);
        V.emplace_back(v2);
        V.emplace_back(v3);

        // append flat normal
        glm::vec3 triangle_normal = glm::normalize(glm::cross(v2 - v1, v3 - v1));
        N_f.emplace_back(triangle_normal);
        N_f.emplace_back(triangle_normal);
        N_f.emplace_back(triangle_normal);

        // append mid of triangles
        glm::vec3 mid = (v1 + v2 + v3)/3.0f;
        M.emplace_back(mid);
        M.emplace_back(mid);
        M.emplace_back(mid);

        // append each triangle normal for vertices
        normal_map[v1].emplace_back(triangle_normal);
        normal_map[v2].emplace_back(triangle_normal);
        normal_map[v3].emplace_back(triangle_normal);
    }

    // get the start of the current figure in all VBOs
    int start = 0;
    for(int i: numbers) start += i;

    for(int i=0;i<triangle_num*3;i++) {
        glm::vec3 cur = V.at(start+i);
        glm::vec3 phong_normal = glm::vec3(0.0f);
        for(auto v: normal_map[cur]) {
            phong_normal += v;
        }
        phong_normal /= normal_map[cur].size();

        // append phong normal
        N_p.emplace_back(glm::normalize(phong_normal));
    }

    base_color.emplace_back(1.0f);
    numbers.emplace_back(3*triangle_num);
    model.emplace_back(1.0f);
    render_mode.emplace_back(3);

    in.close();
}

GLuint load_cube() {
    GLuint id;
    glGenTextures(1, &id);
    glBindTexture(GL_TEXTURE_CUBE_MAP, id);

    std::string suffixes[] = {"posx", "negx", "posy", "negy", "posz", "negz"};
    GLuint targets[] = {
            GL_TEXTURE_CUBE_MAP_POSITIVE_X, GL_TEXTURE_CUBE_MAP_NEGATIVE_X,
            GL_TEXTURE_CUBE_MAP_POSITIVE_Y, GL_TEXTURE_CUBE_MAP_NEGATIVE_Y,
            GL_TEXTURE_CUBE_MAP_POSITIVE_Z, GL_TEXTURE_CUBE_MAP_NEGATIVE_Z
    };

    GLint w, h, ch;

    for(int i=0;i<6;i++) {
        std::string name = "../data/" + suffixes[i] + ".png";
        unsigned char *data = stbi_load(name.c_str(), &w, &h, &ch, 0);
        glTexImage2D(targets[i], 0, GL_RGB, w, h, 0, GL_RGB, GL_UNSIGNED_BYTE, data);
        stbi_image_free(data);
    }

    glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE);

    return id;
}

void add_plain() {
    // put the base plain
    V.emplace_back(glm::vec3(3, -0.5, 3));
    V.emplace_back(glm::vec3(3, -0.5, -3));
    V.emplace_back(glm::vec3(-3, -0.5, 3));
    V.emplace_back(glm::vec3(-3, -0.5, -3));
    V.emplace_back(glm::vec3(3, -0.5, -3));
    V.emplace_back(glm::vec3(-3, -0.5, 3));

    for(int i=0;i<6;i++) {
        N_f.emplace_back(glm::vec3(0, 1, 0));
        M.emplace_back(glm::vec3(0, -0.5, 0));
        N_p.emplace_back(glm::vec3(0, 1, 0));
    }
    VBO.update(V);
    VBO_M.update(M);
    VBO_N_f.update(N_f);
    VBO_N_p.update(N_p);

    base_color.emplace_back(1.0f);
    model.emplace_back(1.0f);
    numbers.emplace_back(6);
    render_mode.emplace_back(3);
}

void add_cube(int scale) {
    std::unordered_map<glm::vec3, std::vector<glm::vec3>, vecHash> normal_map;
    int start = 0;

    for (int i = 0; i < 108; i += 9) {
        glm::vec3 v1 = glm::vec3(unit_cube[i], unit_cube[i + 1], unit_cube[i + 2]).operator*=(scale);
        glm::vec3 v2 = glm::vec3(unit_cube[i + 3], unit_cube[i + 4], unit_cube[i + 5]).operator*=(scale);
        glm::vec3 v3 = glm::vec3(unit_cube[i + 6], unit_cube[i + 7], unit_cube[i + 8]).operator*=(scale);

        V.emplace_back(v1);
        V.emplace_back(v2);
        V.emplace_back(v3);

        // append flat normal
        glm::vec3 triangle_normal = glm::normalize(glm::cross(v2 - v1, v3 - v1));
        N_f.emplace_back(triangle_normal);
        N_f.emplace_back(triangle_normal);
        N_f.emplace_back(triangle_normal);

        // append mid of triangles
        glm::vec3 mid = (v1 + v2 + v3)/3.0f;
        M.emplace_back(mid);
        M.emplace_back(mid);
        M.emplace_back(mid);

        normal_map[v1].emplace_back(triangle_normal);
        normal_map[v2].emplace_back(triangle_normal);
        normal_map[v3].emplace_back(triangle_normal);
    }

    for(int i: numbers) start += i;

    for(int i=0;i<36;i++) {
        glm::vec3 cur = V.at(start+i);
        glm::vec3 phong_normal = glm::vec3(0.0f);
        for(auto v: normal_map[cur]) {
            phong_normal += v;
        }
        phong_normal /= normal_map[cur].size();

        N_p.emplace_back(glm::normalize(phong_normal));
    }

    base_color.emplace_back(1.0f);
    model.emplace_back(1.0f);
    numbers.emplace_back(36);
    render_mode.emplace_back(3);

    VBO.update(V);
    VBO_M.update(M);
    VBO_N_f.update(N_f);
    VBO_N_p.update(N_p);
}

// Callbacks
void framebuffer_size_callback(GLFWwindow* window, int width, int height) {
    glViewport(0, 0, width, height);
    window_width = width;
    window_height = height;

    if(is_perspective) {
        projection = glm::perspective(glm::radians(45.0f), window_width/window_height, 0.1f, 100.0f);
    } else {
        projection = glm::ortho(-window_width/window_height, window_width/window_height, -1.0f, 1.0f, 0.1f, 100.0f);
    }
}

void mouse_button_callback(GLFWwindow* window, int button, int action, int mods) {
    // Get the position of the mouse in the window
    double xpos, ypos;
    glfwGetCursorPos(window, &xpos, &ypos);

    // Get the size of the window
    int width, height;
    glfwGetWindowSize(window, &width, &height);

    // Convert screen position to world coordinates
    double xworld = ((xpos/double(width))*2)-1;
    double yworld = (((height-1-ypos)/double(height))*2)-1; // NOTE: y axis is flipped in glfw

    if(button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_PRESS) {
        if(is_perspective) {
            // shoot a ray from the camera in world ordinary
            glm::vec4 t = glm::inverse(projection) * glm::vec4(xworld, yworld, -1.0f, 1.0f);
            t = glm::inverse(camera) * glm::vec4(t[0], t[1], -1.0f, 0.0f);
            glm::vec3 ray_direction = glm::normalize(glm::vec3(t));

            int start = 42;
            for (int i = 2; i < numbers.size(); i++) {
                for (int j = 0; j < numbers[i]; j += 3) {
                    // get the figure's current position in world
                    glm::vec3 v1 = model[i] * glm::vec4(V[start + j], 1.0f);
                    glm::vec3 v2 = model[i] * glm::vec4(V[start + j + 1], 1.0f);
                    glm::vec3 v3 = model[i] * glm::vec4(V[start + j + 2], 1.0f);
                    glm::vec3 triangle_normal = glm::normalize(glm::cross(v2 - v1, v3 - v1));

                    if (glm::dot(triangle_normal, ray_direction) == 0) continue;

                    float plane_distance = glm::length(v1) * glm::dot(glm::normalize(v1), triangle_normal);
                    float ray_distance = (plane_distance - glm::dot(eye, triangle_normal)) /
                                         glm::dot(triangle_normal, ray_direction);

                    if (ray_distance < 0) continue;

                    // this is the intersection of ray and plane containing the triangle
                    glm::vec3 ray_plane_inter = eye + ray_distance * ray_direction;

                    // test if the intersection is in the triangle, if true change the base color of the whole figure
                    if (in_triangle(v1, v2, v3, ray_plane_inter)) {
                        base_color[i] = glm::vec3(1.0f, 0.0f, 0.0f);
                        break;
                    }
                }
                start += numbers[i];
            }
        } else {
            // shoot a ray from the camera in world ordinary
            glm::vec4 t = glm::inverse(projection) * glm::vec4(xworld, yworld, -1.0f, 1.0f);
            glm::vec3 ray_direction = glm::normalize(-eye);

            // distance from camera to the origin
            float camera_distance = glm::length(eye);

            // distance from camera's projection on x-z plane to origin
            float camera_xz = sqrt(eye[0]*eye[0] + eye[2]*eye[2]);

            // y position of the clicked point (calculated by the y portion of un-projected yworld)
            float o_y = eye[1] + t[1]*camera_xz/camera_distance;

            // x portion of un-projected yworld
            float o_xz_yworld = t[1]*eye[1]/camera_distance;

            // x and z position of the clicked point (calculated by the x portion of un-projected yworld and xworld)
            float o_x = eye[0] + t[0]*eye[2]/camera_xz - o_xz_yworld*eye[0]/camera_xz;
            float o_z = eye[2] - t[0]*eye[0]/camera_xz - o_xz_yworld*eye[2]/camera_xz;

            // ray origin (clicked point)
            glm::vec3 origin = glm::vec3(o_x, o_y, o_z);

            int start = 42;
            for(int i=2;i<numbers.size();i++) {
                for(int j=0;j<numbers[i];j+=3) {
                    // get the figure's current position in world
                    glm::vec3 v1 = model[i] * glm::vec4(V[start+j], 1.0f);
                    glm::vec3 v2 = model[i] * glm::vec4(V[start+j+1], 1.0f);
                    glm::vec3 v3 = model[i] * glm::vec4(V[start+j+2], 1.0f);
                    glm::vec3 triangle_normal = glm::normalize(glm::cross(v2 - v1, v3 - v1));

                    if(glm::dot(triangle_normal, ray_direction) == 0) continue;

                    float plane_distance = glm::length(v1) * glm::dot(glm::normalize(v1), triangle_normal);
                    float ray_distance = (plane_distance - glm::dot(origin, triangle_normal))/glm::dot(triangle_normal, ray_direction);

                    if(ray_distance < 0) continue;

                    // this is the intersection of ray and plane containing the triangle
                    glm::vec3 ray_plane_inter = origin + ray_distance * ray_direction;

                    // test if the intersection is in the triangle, if true change the base color of the whole figure
                    if(in_triangle(v1, v2, v3, ray_plane_inter)) {
                        base_color[i] = glm::vec3(1.0f, 0.0f, 0.0f);
                        break;
                    }
                }
                start += numbers[i];
            }
        }
    } else if(button == GLFW_MOUSE_BUTTON_RIGHT && action == GLFW_PRESS) {
        // right button to remove all selection
        for(int i=2;i<base_color.size();i++) {
            base_color[i] = glm::vec3(1.0f);
        }
    }

    // Upload the change to the GPU
    if(!V.empty()) {
        VBO.update(V);
        VBO_M.update(M);
        VBO_N_f.update(N_f);
        VBO_N_p.update(N_p);
    }
}

void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods) {
    // Get the position of the mouse in the window
    double xpos, ypos;
    glfwGetCursorPos(window, &xpos, &ypos);

    // Get the size of the window
    int width, height;
    glfwGetWindowSize(window, &width, &height);

    // Convert screen position to world coordinates
    double xworld = ((xpos/double(width))*2)-1;
    double yworld = (((height-1-ypos)/double(height))*2)-1; // NOTE: y axis is flipped in glfw

    // Update the mode
    if(action == GLFW_PRESS) {
        switch (key) {
            case GLFW_KEY_1:
                add_cube(1);
                break;
            case GLFW_KEY_2:
                load_model("../data/bumpy_cube.off");
                break;
            case GLFW_KEY_3:
                load_model("../data/bunny.off");
                break;
            case GLFW_KEY_W:
                for(int i=0;i<render_mode.size();i++) {
                    if(base_color[i] == glm::vec3(1.0f, 0.0f, 0.0f)) {
                        render_mode[i] = 1;
                    }
                }
                break;
            case GLFW_KEY_F:
                for(int i=0;i<render_mode.size();i++) {
                    if(base_color[i] == glm::vec3(1.0f, 0.0f, 0.0f)) {
                        render_mode[i] = 2;
                    }
                }
                break;
            case GLFW_KEY_P:
                for(int i=0;i<render_mode.size();i++) {
                    if(base_color[i] == glm::vec3(1.0f, 0.0f, 0.0f)) {
                        render_mode[i] = 3;
                    }
                }
                break;
            case GLFW_KEY_K:
                for(int i=0;i<render_mode.size();i++) {
                    if(base_color[i] == glm::vec3(1.0f, 0.0f, 0.0f)) {
                        render_mode[i] = 4;
                    }
                }
                break;
            case GLFW_KEY_L:
                for(int i=0;i<render_mode.size();i++) {
                    if(base_color[i] == glm::vec3(1.0f, 0.0f, 0.0f)) {
                        render_mode[i] = 5;
                    }
                }
                break;
            case GLFW_KEY_Q:
                for(int i=0;i<render_mode.size();i++) {
                    if(base_color[i] == glm::vec3(1.0f, 0.0f, 0.0f)) {
                        glm::mat4 translate = glm::translate(glm::mat4(1.0f), glm::vec3(-1.0f, 0.0f, 0.0f));
                        model[i] = translate * model[i];
                    }
                }
                break;
            case GLFW_KEY_E:
                for(int i=0;i<render_mode.size();i++) {
                    if(base_color[i] == glm::vec3(1.0f, 0.0f, 0.0f)) {
                        glm::mat4 translate = glm::translate(glm::mat4(1.0f), glm::vec3(1.0f, 0.0f, 0.0f));
                        model[i] = translate * model[i];
                    }
                }
                break;
            case GLFW_KEY_A:
                for(int i=0;i<render_mode.size();i++) {
                    if(base_color[i] == glm::vec3(1.0f, 0.0f, 0.0f)) {
                        glm::mat4 translate = glm::translate(glm::mat4(1.0f), glm::vec3(0.0f, -1.0f, 0.0f));
                        model[i] = translate * model[i];
                    }
                }
                break;
            case GLFW_KEY_D:
                for(int i=0;i<render_mode.size();i++) {
                    if(base_color[i] == glm::vec3(1.0f, 0.0f, 0.0f)) {
                        glm::mat4 translate = glm::translate(glm::mat4(1.0f), glm::vec3(0.0f, 1.0f, 0.0f));
                        model[i] = translate * model[i];
                    }
                }
                break;
            case GLFW_KEY_Z:
                for(int i=0;i<render_mode.size();i++) {
                    if(base_color[i] == glm::vec3(1.0f, 0.0f, 0.0f)) {
                        glm::mat4 translate = glm::translate(glm::mat4(1.0f), glm::vec3(0.0f, 0.0f, -1.0f));
                        model[i] = translate * model[i];
                    }
                }
                break;
            case GLFW_KEY_C:
                for(int i=0;i<render_mode.size();i++) {
                    if(base_color[i] == glm::vec3(1.0f, 0.0f, 0.0f)) {
                        glm::mat4 translate = glm::translate(glm::mat4(1.0f), glm::vec3(0.0f, 0.0f, 1.0f));
                        model[i] = translate * model[i];
                    }
                }
                break;
            case GLFW_KEY_T:
                for(int i=0;i<render_mode.size();i++) {
                    if(base_color[i] == glm::vec3(1.0f, 0.0f, 0.0f)) {
                        glm::mat4 rotate = glm::rotate(glm::mat4(1.0f), glm::radians(45.0f), glm::vec3(1.0f, 0.0f, 0.0f));
                        model[i] = model[i] * rotate;
                    }
                }
                break;
            case GLFW_KEY_U:
                for(int i=0;i<render_mode.size();i++) {
                    if(base_color[i] == glm::vec3(1.0f, 0.0f, 0.0f)) {
                        glm::mat4 rotate = glm::rotate(glm::mat4(1.0f), glm::radians(-45.0f), glm::vec3(1.0f, 0.0f, 0.0f));
                        model[i] = model[i] * rotate;
                    }
                }
                break;
            case GLFW_KEY_G:
                for(int i=0;i<render_mode.size();i++) {
                    if(base_color[i] == glm::vec3(1.0f, 0.0f, 0.0f)) {
                        glm::mat4 rotate = glm::rotate(glm::mat4(1.0f), glm::radians(45.0f), glm::vec3(0.0f, 1.0f, 0.0f));
                        model[i] = model[i] * rotate;
                    }
                }
                break;
            case GLFW_KEY_J:
                for(int i=0;i<render_mode.size();i++) {
                    if(base_color[i] == glm::vec3(1.0f, 0.0f, 0.0f)) {
                        glm::mat4 rotate = glm::rotate(glm::mat4(1.0f), glm::radians(-45.0f), glm::vec3(0.0f, 1.0f, 0.0f));
                        model[i] = model[i] * rotate;
                    }
                }
                break;
            case GLFW_KEY_B:
                for(int i=0;i<render_mode.size();i++) {
                    if(base_color[i] == glm::vec3(1.0f, 0.0f, 0.0f)) {
                        glm::mat4 rotate = glm::rotate(glm::mat4(1.0f), glm::radians(45.0f), glm::vec3(0.0f, 0.0f, 1.0f));
                        model[i] = model[i] * rotate;
                    }
                }
                break;
            case GLFW_KEY_M:
                for(int i=0;i<render_mode.size();i++) {
                    if(base_color[i] == glm::vec3(1.0f, 0.0f, 0.0f)) {
                        glm::mat4 rotate = glm::rotate(glm::mat4(1.0f), glm::radians(-45.0f), glm::vec3(0.0f, 0.0f, 1.0f));
                        model[i] = model[i] * rotate;
                    }
                }
                break;
            case GLFW_KEY_I:
                for(int i=0;i<render_mode.size();i++) {
                    if(base_color[i] == glm::vec3(1.0f, 0.0f, 0.0f)) {
                        glm::mat4 scale = glm::scale(glm::mat4(1.0f), glm::vec3(1.25f));
                        model[i] = model[i] * scale;
                    }
                }
                break;
            case GLFW_KEY_O:
                for(int i=0;i<render_mode.size();i++) {
                    if(base_color[i] == glm::vec3(1.0f, 0.0f, 0.0f)) {
                        glm::mat4 scale = glm::scale(glm::mat4(1.0f), glm::vec3(0.75f));
                        model[i] = model[i] * scale;
                    }
                }
                break;
            case GLFW_KEY_4:
                eye = glm::rotate(glm::mat4(1.0f), glm::radians(45.0f), glm::vec3(1.0f, 0.0f, 0.0f)) * glm::vec4(eye, 1.0f);
                camera = glm::lookAt(eye, glm::vec3(0.0f, 0.0f, 0.0f), glm::vec3(0.0f, 1.0f, 0.0f));
                break;
            case GLFW_KEY_5:
                eye = glm::rotate(glm::mat4(1.0f), glm::radians(-45.0f), glm::vec3(1.0f, 0.0f, 0.0f)) * glm::vec4(eye, 1.0f);
                camera = glm::lookAt(eye, glm::vec3(0.0f, 0.0f, 0.0f), glm::vec3(0.0f, 1.0f, 0.0f));
                break;
            case GLFW_KEY_6:
                eye = glm::rotate(glm::mat4(1.0f), glm::radians(45.0f), glm::vec3(0.0f, 1.0f, 0.0f)) * glm::vec4(eye, 1.0f);
                camera = glm::lookAt(eye, glm::vec3(0.0f, 0.0f, 0.0f), glm::vec3(0.0f, 1.0f, 0.0f));
                break;
            case GLFW_KEY_7:
                eye = glm::rotate(glm::mat4(1.0f), glm::radians(-45.0f), glm::vec3(0.0f, 1.0f, 0.0f)) * glm::vec4(eye, 1.0f);
                camera = glm::lookAt(eye, glm::vec3(0.0f, 0.0f, 0.0f), glm::vec3(0.0f, 1.0f, 0.0f));
                break;
            case GLFW_KEY_8:
                eye = glm::rotate(glm::mat4(1.0f), glm::radians(45.0f), glm::vec3(0.0f, 0.0f, 1.0f)) * glm::vec4(eye, 1.0f);
                camera = glm::lookAt(eye, glm::vec3(0.0f, 0.0f, 0.0f), glm::vec3(0.0f, 1.0f, 0.0f));
                break;
            case GLFW_KEY_9:
                eye = glm::rotate(glm::mat4(1.0f), glm::radians(-45.0f), glm::vec3(0.0f, 0.0f, 1.0f)) * glm::vec4(eye, 1.0f);
                camera = glm::lookAt(eye, glm::vec3(0.0f, 0.0f, 0.0f), glm::vec3(0.0f, 1.0f, 0.0f));
                break;
            case GLFW_KEY_KP_ADD:
                eye *= 1.25f;
                camera = glm::lookAt(eye, glm::vec3(0.0f, 0.0f, 0.0f), glm::vec3(0.0f, 1.0f, 0.0f));
                break;
            case GLFW_KEY_EQUAL:
                eye *= 1.25f;
                camera = glm::lookAt(eye, glm::vec3(0.0f, 0.0f, 0.0f), glm::vec3(0.0f, 1.0f, 0.0f));
                break;
            case GLFW_KEY_KP_SUBTRACT:
                eye *= 0.75f;
                camera = glm::lookAt(eye, glm::vec3(0.0f, 0.0f, 0.0f), glm::vec3(0.0f, 1.0f, 0.0f));
                break;
            case GLFW_KEY_MINUS:
                eye *= 0.75f;
                camera = glm::lookAt(eye, glm::vec3(0.0f, 0.0f, 0.0f), glm::vec3(0.0f, 1.0f, 0.0f));
                break;
            case GLFW_KEY_SPACE:
                if(is_perspective) {
                    projection = glm::ortho(-window_width/window_height, window_width/window_height, -1.0f, 1.0f, 0.1f, 100.0f);
                    is_perspective = false;
                } else {
                    projection = glm::perspective(glm::radians(45.0f), window_width/window_height, 0.1f, 100.0f);
                    is_perspective = true;
                }
                break;
            case GLFW_KEY_ESCAPE:
                V.clear();
                N_f.clear();
                N_p.clear();
                M.clear();
                base_color.clear();
                model.clear();
                numbers.clear();

                add_cube(2);
                add_plain();
                break;
            case GLFW_KEY_S:
                if(shadow_mode == 1)
                    shadow_mode = 0;
                else
                    shadow_mode = 1;
                break;
            default:
                break;
        }
    }

    // Upload the change to the GPU
    if(!V.empty()) {
        VBO.update(V);
        VBO_M.update(M);
        VBO_N_f.update(N_f);
        VBO_N_p.update(N_p);
    }
}

int main() {
    GLFWwindow* window;

    // Initialize the library
    if (!glfwInit())
        return -1;

    // Activate supersampling
    glfwWindowHint(GLFW_SAMPLES, 8);

    // Ensure that we get at least a 3.2 context
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 2);

    // On apple we have to load a core profile with forward compatibility
    #ifdef __APPLE__
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
    #endif

    // Create a windowed mode window and its OpenGL context
    window = glfwCreateWindow(640, 480, "Main", nullptr, nullptr);
    if (!window) {
        glfwTerminate();
        return -1;
    }

    // Make the window's context current
    glfwMakeContextCurrent(window);

    #ifndef __APPLE__
      glewExperimental = true;
      GLenum err = glewInit();
      if(GLEW_OK != err) {
        /* Problem: glewInit failed, something is seriously wrong. */
       fprintf(stderr, "Error: %s\n", glewGetErrorString(err));
      }
      glGetError(); // pull and savely ignonre unhandled errors like GL_INVALID_ENUM
      fprintf(stdout, "Status: Using GLEW %s\n", glewGetString(GLEW_VERSION));
    #endif

    int major, minor, rev;
    major = glfwGetWindowAttrib(window, GLFW_CONTEXT_VERSION_MAJOR);
    minor = glfwGetWindowAttrib(window, GLFW_CONTEXT_VERSION_MINOR);
    rev = glfwGetWindowAttrib(window, GLFW_CONTEXT_REVISION);
    printf("OpenGL version recieved: %d.%d.%d\n", major, minor, rev);
    printf("Supported OpenGL is %s\n", (const char*)glGetString(GL_VERSION));
    printf("Supported GLSL is %s\n", (const char*)glGetString(GL_SHADING_LANGUAGE_VERSION));

    // Initialize the VAO
    // A Vertex Array Object (or VAO) is an object that describes how the vertex
    // attributes are stored in a Vertex Buffer Object (or VBO). This means that
    // the VAO is not the actual object storing the vertex data,
    // but the descriptor of the vertex data.
    VertexArrayObject VAO;
    VAO.init();
    VAO.bind();

    // Initialize the VBO with the vertices data
    // A VBO is a data container that lives in the GPU memory
    VBO.init();
    VBO_M.init();
    VBO_N_f.init();
    VBO_N_p.init();

    // Initialize the OpenGL Program
    // A program controls the OpenGL pipeline and it must contains
    // at least a vertex shader and a fragment shader to be valid
    Program program;
    const GLchar* vertex_shader =
            "#version 420 core\n"
                    "in vec3 position;"
                    "in vec3 mid;"
                    "in vec3 normal_f;"
                    "in vec3 normal_p;"
                    "out vec3 vertex_normal_camera;"
                    "out vec3 light_to_vertex_camera;"
                    "out vec3 pos_model;"
                    "out vec3 normal_model;"
                    "out vec4 vertex_light;"
                    "uniform mat4 model;"
                    "uniform mat4 camera;"
                    "uniform mat4 projection;"
                    "uniform mat4 light_vp;"
                    "uniform vec3 light_pos;"
                    "uniform int flat_phong;"
                    "void main()"
                    "{"
                    "    pos_model = (model * vec4(position, 1.0)).xyz;"
                    "    normal_model = normalize((model * vec4(normal_p, 1.0)).xyz);"
                    "    gl_Position = projection * camera * vec4(pos_model, 1.0);"
                    "    vec3 vertex_camera = (camera * model * vec4(position, 1.0)).xyz;"
                    "    vec3 mid_camera = (camera * model * vec4(mid, 1.0)).xyz;"
                    "    mat2x3 point = mat2x3(mid_camera, vertex_camera);"
                    "    vec3 light_camera = (camera * vec4(light_pos, 1.0)).xyz;"
                    "    light_to_vertex_camera = light_camera - point[flat_phong];"
                    "    mat2x3 normal = mat2x3(normal_f, normal_p);"
                    "    vertex_normal_camera = (camera * model * vec4(normal[flat_phong], 0.0)).xyz;"
                    "    vertex_light = light_vp * model * vec4(position, 1.0);"
                    "}";
    const GLchar* fragment_shader =
            "#version 420 core\n"
                    "in vec3 vertex_normal_camera;"
                    "in vec3 light_to_vertex_camera;"
                    "in vec3 pos_model;"
                    "in vec3 normal_model;"
                    "in vec4 vertex_light;"
                    "out vec4 outColor;"
                    "uniform int shadow_mode;"
                    "uniform int render_mode;"
                    "uniform vec3 base_color;"
                    "uniform vec3 eye;"
                    "layout(binding=0) uniform sampler2D shadowMap;"
                    "layout(binding=1) uniform samplerCube skybox;"
                    ""
                    "float ShadowCalculation(vec4 fragPosLightSpace, float cos)"
                    "{"
                    "    vec3 projCoords = fragPosLightSpace.xyz / fragPosLightSpace.w;"
                    "    if(projCoords.z > 1.0) return 0.0;"
                    "    projCoords = projCoords * 0.5 + 0.5;"
                    "    float closestDepth = texture(shadowMap, projCoords.xy).r;"
                    "    float currentDepth = projCoords.z;"
                    "    float bias = max(0.05 * (1.0 - cos), 0.005); "
                    "    float shadow = currentDepth - bias > closestDepth  ? 1.0 : 0.0;"
                    "    return shadow;"
                    "}"
                    ""
                    "void main()"
                    "{"
                    "    vec3 n = normalize(vertex_normal_camera);"
                    "    vec3 l = normalize(light_to_vertex_camera);"
                    "    float cosTheta = clamp(dot(n, l), 0.0, 1.0);"
                    "    float shadow = ShadowCalculation(vertex_light, cosTheta);"
                    "    if(shadow == 1.0) {"
                    "        if(shadow_mode == 1)"
                    "            outColor = vec4(0.0, 0.0, 0.0, 1.0) + vec4(0.2*base_color, 1.0);"
                    "        else "
                    "            outColor = vec4(1.0, 0.0, 0.0, 1.0) + vec4(0.2*base_color, 1.0);"
                    "    } else {"
                    "        if(render_mode == 4) {"
                    "            vec3 i = normalize(pos_model - eye);"
                    "            vec3 r = reflect(i, normal_model);"
                    "            outColor = vec4(texture(skybox, r).rgb, 1.0);"
                    "        } else if(render_mode == 5) {"
                    "            vec3 i = normalize(pos_model - eye);"
                    "            vec3 r = refract(i, normal_model, 0.658);"
                    "            outColor = vec4(texture(skybox, r).rgb, 1.0);"
                    "        } else {"
                    "            outColor = vec4((10 * cosTheta/(length(light_to_vertex_camera)*length(light_to_vertex_camera)) + 0.2) * base_color, 1.0);"
                    "        }"
                    "    }"
                    "}";

    // Compile the two shaders and upload the binary to the GPU
    // Note that we have to explicitly specify that the output "slot" called outColor
    // is the one that we want in the fragment buffer (and thus on screen)
    program.init(vertex_shader,fragment_shader,"outColor");

    Program depth;
    const GLchar* vertex_shader_depth =
            "#version 150 core\n"
            "in vec3 position;"
            "uniform mat4 model;"
            "uniform mat4 light_vp;"
            "void main()"
            "{"
            "    gl_Position = light_vp * model * vec4(position, 1.0);"
            "}";
    const GLchar* fragment_shader_depth =
            "#version 150 core\n"
            "out float fragmentdepth;"
            "void main()"
            "{"
            "     fragmentdepth = gl_FragCoord.z;"
            "}";

    depth.init(vertex_shader_depth, fragment_shader_depth, "fragmentdepth");

    Program skybox;
    const GLchar* vertex_shader_skybox =
            "#version 150 core\n"
            "in vec3 position;"
            "out vec3 tex_pos;"
            "uniform mat4 camera;"
            "uniform mat4 projection;"
            "void main()"
            "{"
            "    gl_Position = projection * camera * vec4(position, 1.0);"
            "    tex_pos = position;"
            "}";
    const GLchar* fragment_shader_skybox =
            "#version 150 core\n"
            "in vec3 tex_pos;"
            "out vec4 outColor;"
            "uniform samplerCube cube;"
            "void main()"
            "{"
            "     outColor = texture(cube, tex_pos);"
            "}";

    skybox.init(vertex_shader_skybox, fragment_shader_skybox, "outColor");

    // Register the keyboard callback
    glfwSetKeyCallback(window, key_callback);

    // Register the mouse button callback
    glfwSetMouseButtonCallback(window, mouse_button_callback);

    // Update viewport
    glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);

    // enable depth rendering
    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LESS);

    // save previous number of vertices
    unsigned short pre_size = V.size();

    // Save the current time --- it will be used to dynamically change light position
    auto t_start = std::chrono::high_resolution_clock::now();

    // depth map
    unsigned int depthMapFBO;
    glGenFramebuffers(1, &depthMapFBO);

    const unsigned int SHADOW_WIDTH = 1024, SHADOW_HEIGHT = 1024;

    unsigned int depthMap;
    glGenTextures(1, &depthMap);
    glBindTexture(GL_TEXTURE_2D, depthMap);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_DEPTH_COMPONENT,
                 SHADOW_WIDTH, SHADOW_HEIGHT, 0, GL_DEPTH_COMPONENT, GL_FLOAT, nullptr);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);

    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_BORDER);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_BORDER);
    float borderColor[] = { 1.0f, 1.0f, 1.0f, 1.0f };
    glTexParameterfv(GL_TEXTURE_2D, GL_TEXTURE_BORDER_COLOR, borderColor);

    glBindFramebuffer(GL_FRAMEBUFFER, depthMapFBO);
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_TEXTURE_2D, depthMap, 0);
    glDrawBuffer(GL_NONE);
    glReadBuffer(GL_NONE);
    glBindFramebuffer(GL_FRAMEBUFFER, 0);

    // skybox
    GLuint box = load_cube();

    add_cube(2);
    add_plain();

    // Loop until the user closes the window
    while (!glfwWindowShouldClose(window)) {
        // Set the light position depending on the time difference
        auto t_now = std::chrono::high_resolution_clock::now();
        float time = std::chrono::duration_cast<std::chrono::duration<float>>(t_now - t_start).count();
        light_source = glm::vec3(std::sin(time), 4, std::cos(time));
        glm::mat4 light_view = glm::lookAt(light_source, glm::vec3(0.0f), glm::vec3(0.0f, 1.0f, 0.0f));
        glm::mat4 light_vp = light_projection * light_view;

        // only pass VBOs to gpu when number of vertices changed
        if(!V.empty() && V.size() != pre_size) {
            depth.bindVertexAttribArray("position", VBO);

            skybox.bindVertexAttribArray("position", VBO);

            program.bindVertexAttribArray("position", VBO);
            program.bindVertexAttribArray("mid", VBO_M);
            program.bindVertexAttribArray("normal_f", VBO_N_f);
            program.bindVertexAttribArray("normal_p", VBO_N_p);
        }
        pre_size = V.size();

        // Clear the framebuffer
        glClearColor(0.5f, 0.5f, 0.5f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        int start = 36;

        // bind depth
        depth.bind();
        glCullFace(GL_FRONT);
        glUniformMatrix4fv(depth.uniform("light_vp"), 1, GL_FALSE, glm::value_ptr(light_vp));

        glViewport(0, 0, SHADOW_WIDTH, SHADOW_HEIGHT);
        glBindFramebuffer(GL_FRAMEBUFFER, depthMapFBO);
        glClear(GL_DEPTH_BUFFER_BIT);
        for(int i=1;i<numbers.size();i++) {
            glUniformMatrix4fv(depth.uniform("model"), 1, GL_FALSE, glm::value_ptr(model[i]));
            for (int j = 0; j < numbers[i]; j += 3) {
                glDrawArrays(GL_TRIANGLES, start + j, 3);
            }
            start += numbers[i];
        }
        glBindFramebuffer(GL_FRAMEBUFFER, 0);
        glViewport(0, 0, window_width, window_height);
        glCullFace(GL_BACK);

        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        // bind skybox
        skybox.bind();
        glDepthMask(GL_FALSE);
        glActiveTexture(GL_TEXTURE0);
        glBindTexture(GL_TEXTURE_CUBE_MAP, box);
        glUniformMatrix4fv(skybox.uniform("projection"), 1, GL_FALSE, glm::value_ptr(projection));
        glUniformMatrix4fv(skybox.uniform("camera"), 1, GL_FALSE, glm::value_ptr(glm::mat4(glm::mat3(camera))));
        glDrawArrays(GL_TRIANGLES, 0, 36);
        glDepthMask(GL_TRUE);

        // Bind your program
        program.bind();
        glBindTextureUnit(0, depthMap);
        glBindTextureUnit(1, box);

        // Bind the camera and projection and light source
        glUniformMatrix4fv(program.uniform("projection"), 1, GL_FALSE, glm::value_ptr(projection));
        glUniformMatrix4fv(program.uniform("camera"), 1, GL_FALSE, glm::value_ptr(camera));
        glUniformMatrix4fv(program.uniform("light_vp"), 1, GL_FALSE, glm::value_ptr(light_vp));
        glUniform3fv(program.uniform("light_pos"), 1, glm::value_ptr(light_source));
        glUniform1i(program.uniform("shadow_mode"), shadow_mode);
        glUniform3fv(program.uniform("eye"), 1, glm::value_ptr(eye));

        start = 36;
        for(int i=1;i<numbers.size();i++) {
            glUniform3fv(program.uniform("base_color"), 1, glm::value_ptr(base_color[i]));
            glUniformMatrix4fv(program.uniform("model"), 1, GL_FALSE, glm::value_ptr(model[i]));
            glUniform1i(program.uniform("render_mode"), render_mode[i]);
            if(render_mode[i] == 1) {
                glUniform1i(program.uniform("flat_phong"), 0);
                for (int j = 0; j < numbers[i]; j += 3) {
                    glDrawArrays(GL_LINE_LOOP, start + j, 3);
                }
            } else if(render_mode[i] == 2) {
                glUniform1i(program.uniform("flat_phong"), 0);
                for (int j = 0; j < numbers[i]; j += 3) {
                    glDrawArrays(GL_TRIANGLES, start + j, 3);
                    glDrawArrays(GL_LINE_LOOP, start + j, 3);
                }
            } else {
                glUniform1i(program.uniform("flat_phong"), 1);
                for (int j = 0; j < numbers[i]; j += 3) {
                    glDrawArrays(GL_TRIANGLES, start + j, 3);
                }
            }
            start += numbers[i];
        }

        // Swap front and back buffers
        glfwSwapBuffers(window);

        // Poll for and process events
        glfwPollEvents();
    }

    // Deallocate opengl memory
    program.free();
    depth.free();
    skybox.free();
    VAO.free();
    VBO.free();
    VBO_M.free();
    VBO_N_f.free();
    VBO_N_p.free();

    // Deallocate glfw internals
    glfwTerminate();
    return 0;
}
