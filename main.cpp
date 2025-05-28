#include "lib/gif/gif.h"
#include <cstdlib>
#include <ctime>
#include <vector>
#include <string>
#include <sstream>
#include <iostream>
#include <cmath>
#include <eigen/Eigen/Core>
#include <algorithm>
#include "lib/SimplexNoise/SimplexNoise.h"

const int SCALE = 8;
const int WIDTH = pow(2, SCALE);
const int HEIGHT = pow(2, SCALE);
const int NUM_FRAMES = 120;
const int NUM_GIFS = 1;
const int RADIUS = pow(2, SCALE - 2);
const int SHADING_LEVELS = -1; // -1 for smooth shading

// Global Random Variables
double CLOUD_PROBABILTIY;
double CLOUD_DENSITY;

std::vector<double> generate_normal_map(int r)
{
    std::vector<double> normal(WIDTH * HEIGHT * 3); // RGBA
    int midw = WIDTH / 2;
    int midh = HEIGHT / 2;
    for (int py = 0; py < HEIGHT; py++)
    {
        for (int px = 0; px < WIDTH; px++)
        {
            int idx = ((py)*WIDTH + (px)) * 3;

            // calculate the coords as if (0,0) is the center of the image
            double sphere_x = px - midw;
            double sphere_y = py - midh;

            // shrink all the points so that its a unit circle (radius is 1)
            double nx = sphere_x / r;
            double ny = sphere_y / r;

            // calculate the distance each point is from the center of the circle (minus the final sqrt)
            double dist2 = nx * nx + ny * ny;

            // temporarily make nz the difference between the dist and the radius (just for convenience)
            double nz = 1 - dist2;

            // if we are outside the circle/spheres radius then mark this point as outside using an impossible number (we should never have a normal vector with a value < -1)
            if (nz < 0)
            {
                normal[idx + 0] = -5;
                continue;
            }

            // calculate the actual nz (by finally adding in that sqrt)
            nz = std::sqrt(nz);

            normal[idx + 0] = nx;
            normal[idx + 1] = ny;
            normal[idx + 2] = nz;
        }
    }
    return normal;
}

Eigen::Matrix3d rotation_matrix(int rot_idx)
{

    double angle = rot_idx * 2 / double(NUM_FRAMES) * M_PI;

    Eigen::Matrix3d y_rot;
    y_rot << std::cos(angle), 0, std::sin(angle),
        0, 1, 0,
        -std::sin(angle), 0, std::cos(angle);

    return y_rot;
}

std::vector<int> get_terrain_value(const SimplexNoise *s, double x, double y, double z)
{
    SimplexNoise simplex = *s;
    double terrain_value = 0.5 * simplex.noise(x, y, z);
    terrain_value += 0.25 * (simplex.noise(x * 4, y * 4, z * 4) + 1) / 2;
    terrain_value += 0.125 * (simplex.noise(x * 8, y * 8, z * 8) + 1) / 2;
    terrain_value += 0.125 * (simplex.noise(x * 16, y * 16, z * 16) + 1) / 2;

    std::vector<int> rgb(3);
    if (terrain_value < 0.3)
    {
        // blue = ocean
        rgb[2] = 255;
    }
    else if (terrain_value < 0.35)
    {
        // yellow = sand
        rgb[0] = 255, rgb[1] = 255;
    }
    else if (terrain_value < 0.45)
    {
        // green = land
        rgb[1] = 255;
    }
    else
    {
        // white = snow
        rgb[0] = 0, rgb[1] = 100, rgb[2] = 0;
    }
    return rgb;
}

double get_cloud_value(const SimplexNoise *s, double x, double y, double z)
{
    SimplexNoise simplex = *s;
    double terrain_value = 0.5 * simplex.noise(x, y, z);
    terrain_value += 0.25 * (simplex.noise(x * 4, y * 4, z * 4) + 1) / 2;
    terrain_value += 0.125 * (simplex.noise(x * 8, y * 8, z * 8) + 1) / 2;
    terrain_value += 0.125 * (simplex.noise(x * 16, y * 16, z * 16) + 1) / 2;
    return terrain_value;
}

// Choose if Global Random Variables are random and if now what value they are
void init_rand_global_variables()
{
    CLOUD_PROBABILTIY = true ? double(rand()) / RAND_MAX : 0.4;
    CLOUD_DENSITY = true ? 2 * double(rand()) / RAND_MAX : 1.5;
}

// Generates a simple animation: a color-shifting square
void generate_gif(const std::string &filename)
{
    // Randomized scalers used to ensure randomness of generation
    int simplexScalerX = rand() % WIDTH - (WIDTH / 2);
    int simplexScalerY = rand() % WIDTH - (WIDTH / 2);
    int simplexScalerZ = rand() % WIDTH - (WIDTH / 2);

    const SimplexNoise simplex(0.01f, 0.5f, 3.0f, 0.5f);
    std::vector<unsigned char> frame(WIDTH * HEIGHT * 4); // RGBA
    // precalculate normal maps for planet and cloud to save on calculation time
    std::vector<double> planet_normal_map = generate_normal_map(int(double(RADIUS) * .95));
    std::vector<double> cloud_normal_map = generate_normal_map(RADIUS);

    // std::vector<unsigned char> texture_map = generate_texture();
    GifWriter writer;
    GifBegin(&writer, filename.c_str(), WIDTH, HEIGHT, 4);

    int midw = WIDTH / 2;
    int midh = HEIGHT / 2;

    for (int frame_idx = 0; frame_idx < NUM_FRAMES; ++frame_idx)
    {
        frame.clear();

        // loop through assuming (0,0) os the center of the circle
        for (int x = -RADIUS; x < RADIUS; x++)
        {
            for (int y = -RADIUS; y < RADIUS; y++)
            {

                // since sin & cos go back and forth between 1 & -1 we can use these to smoothly loop overtime
                double cos_loop = std::cos(double(frame_idx) * 2 / double(NUM_FRAMES) * M_PI);
                double sin_loop = std::sin(double(frame_idx) * 2 / double(NUM_FRAMES) * M_PI);

                // calculate where the sphere pixels will be on the actual image and the corresponding index for the array
                double px = x + midw;
                double py = y + midh;

                // set up index for the final image and also grab the normal vector values from our map
                int idx = ((py)*WIDTH + (px)) * 4;
                int nidx = ((py)*WIDTH + (px)) * 3;
                double planet_nx = planet_normal_map[nidx + 0];
                double planet_ny = planet_normal_map[nidx + 1];
                double planet_nz = planet_normal_map[nidx + 2];

                double cloud_nx = cloud_normal_map[nidx + 0];
                double cloud_ny = cloud_normal_map[nidx + 1];
                double cloud_nz = cloud_normal_map[nidx + 2];

                // if the normal vector is <-1 we are outside of the cloud radius and planet radius (since cloud > planet) so just skip the pixel
                if (cloud_nx < -1)
                {
                    continue;
                }

                // make our actual vector and normalize it just in case
                Eigen::Vector3d planet_normal_vec(planet_nx, planet_ny, planet_nz);
                Eigen::Vector3d cloud_normal_vec(cloud_nx, cloud_ny, cloud_nz);
                planet_normal_vec.normalize();
                cloud_normal_vec.normalize();

                // rotate our vector according by what frame we are on
                Eigen::Vector3d planet_rotated_normal_vec;
                Eigen::Vector3d cloud_rotated_normal_vec;
                planet_rotated_normal_vec = rotation_matrix(frame_idx) * planet_normal_vec;
                cloud_rotated_normal_vec = rotation_matrix(frame_idx) * cloud_normal_vec;

                // create a biased normal vector (x,y,z) = ([0,1], [0,1], [0,1])
                // Eigen::Vector3d planet_biased_normal_vec(0.5 * (planet_rotated_normal_vec[0] + 1), 0.5 * (planet_rotated_normal_vec[1] + 1), planet_rotated_normal_vec[2]);

                // light calculations
                double light_intensity = 10;
                Eigen::Vector3d light_direction(std::cos(double(frame_idx) * 4 / double(NUM_FRAMES) * M_PI) * 2, 0, std::sin(double(frame_idx) * 4 / double(NUM_FRAMES) * M_PI) * 2);
                light_direction.normalize();
                double light_power = planet_normal_vec.dot(light_direction) * light_intensity;
                double diffuse = (0 < light_power) ? (light_power > light_intensity) ? light_intensity : light_power : 0;
                diffuse *= 0.8 / light_intensity;
                if (SHADING_LEVELS != -1)
                {
                    double shading_step_size = 0.8 / (SHADING_LEVELS - 1);
                    diffuse = double(int((diffuse + shading_step_size / 2) / shading_step_size));
                    diffuse = diffuse * shading_step_size;
                }
                double ambient = 0.2f;
                double shade = diffuse + ambient;

                // calculate the probability that this is a cloud using random scaling and cyclic noise so that the clouds move in smooth way and are random
                double cloud_val = get_cloud_value(&simplex,
                                                   cloud_rotated_normal_vec[0] + simplexScalerX - cos_loop,
                                                   cloud_rotated_normal_vec[1] + simplexScalerY - sin_loop,
                                                   cloud_rotated_normal_vec[2] + simplexScalerZ - sin_loop);

                // determine if this point is a cloud, using some cyclic noise so the size of the clouds grow and shrink
                double is_cloud = cloud_val >= CLOUD_PROBABILTIY + (sin_loop / 10);

                std::vector<int> terrain_color(3);

                // If the x value of the normal vector is less than -1 (not normalized) we must be outside of the planet's radius
                if (!(planet_nx < -1))
                {
                    terrain_color = get_terrain_value(&simplex,
                                                      planet_rotated_normal_vec[0] - simplexScalerX,
                                                      planet_rotated_normal_vec[1] - simplexScalerY,
                                                      planet_rotated_normal_vec[2] - simplexScalerZ);
                }

                // calculate how thick the clouds should be at this point
                double cloud_opacity = (1 < cloud_val * CLOUD_DENSITY) ? 1 : (0 > cloud_val * CLOUD_DENSITY) ? 0
                                                                                                             : cloud_val * CLOUD_DENSITY;

                // if this pixel is determined to have a cloud on it then draw it as so
                // ( (1-p)R1 + p*R2, (1-p)*G1 + p*G2, (1-p)*B1 + p*B2 )
                if (is_cloud)
                {
                    terrain_color[0] = (1 - cloud_opacity) * terrain_color[0] + cloud_opacity * 255;
                    terrain_color[1] = (1 - cloud_opacity) * terrain_color[1] + cloud_opacity * 255;
                    terrain_color[2] = (1 - cloud_opacity) * terrain_color[2] + cloud_opacity * 255;
                }

                frame[idx + 0] = shade * terrain_color[0]; // R
                frame[idx + 1] = shade * terrain_color[1]; // G
                frame[idx + 2] = shade * terrain_color[2]; // B
                frame[idx + 3] = 255;                      // A
            }
        }

        // use one radius value to represent total radius (including clouds), if the dist from center is <0.9 then do planet terrain and cloud coloring,
        // otherwise if <= 1 then only do cloud coloring this way we only have to do one main loop and dont need to handle multiple radius

        // print_frame(frame, frame_idx);

        GifWriteFrame(&writer, frame.data(), WIDTH, HEIGHT, 8);
    }

    GifEnd(&writer);
}

int main()
{
    std::srand(std::time(0));
    init_rand_global_variables();

    for (int i = 0; i < NUM_GIFS; ++i)
    {
        std::ostringstream filename;
        filename << "output_" << i << ".gif";

        std::cout << "Generating " << filename.str() << "...\n";
        generate_gif(filename.str());
    }

    std::cout << "Done.\n";
    return 0;
}
