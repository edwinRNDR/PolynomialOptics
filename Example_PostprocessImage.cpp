/* Polynomial Optics
* (C) 2012, Matthias Hullin <hullin@cs.ubc.ca>
* (University of British Columbia)
*
* This file is in the public domain.
* Different licenses may apply for the included libraries.
*
* ============================
* Example_PostprocessImage.cpp
* ============================
* Example application filtering a HDR image by simulating an
* achromatic doublet lens (Paper figure 13).
* Includes somewhat fancier tricks such as re-arrangement of equations
* to reduce the effort of computing Lambertian cosine term, and
* introduction of wavelength as additional variable in the system.
*/

#include <TruncPoly/TruncPolySystem.hh>

#include <OpticalElements/OpticalMaterial.hh>
#include <OpticalElements/Spherical5.hh>
#include <OpticalElements/Cylindrical5.hh>
#include <OpticalElements/Propagation5.hh>
#include <OpticalElements/TwoPlane5.hh>

#include <OpticalElements/FindFocus.hh>

#include <iostream>
#include <stdlib.h>
#include <chrono>

#define cimg_display 0

#include <include/CImg.h>
#include <include/spectrum.h>

using namespace cimg_library;

void showUsage(char *s);

Transform4f get_system(float lambda, int degree) {
    // Let's simulate Edmund Optics achromat #NT32-921:
    /* Clear Aperture CA (mm) 	39.00
    Eff. Focal Length EFL (mm) 	120.00
    Back Focal Length BFL (mm) 	111.00
    Center Thickness CT 1 (mm) 	9.60
    Center Thickness CT 2 (mm) 	4.20
    Radius R1 (mm) 	65.22
    Radius R2 (mm) 	-62.03
    Radius R3 (mm) 	-1240.67
    Substrate 	N-SSK8/N-SF10
    */

    OpticalMaterial glass1("N-SSK8", true);
    OpticalMaterial glass2("N-SF10", true);

    // Also try: const float d0 = 5000; // Scene is 5m away
    const float d0 = 5000000; // Scene is 5km away
    const float R1 = 65.22;
    const float d1 = 9.60;
    const float R2 = -62.03;
    const float d2 = 4.20;
    const float R3 = -1240.67;

    return two_plane_5(d0, degree)
            >> refract_spherical_5(R1, 1.f, glass1.get_index(lambda), degree)
            >> propagate_5(d1, degree)
            >> refract_spherical_5(R2, glass1.get_index(lambda), glass2.get_index(lambda), degree)
            >> propagate_5(d2, degree)
            >> refract_spherical_5(R3, glass2.get_index(lambda), 1.f, degree);
}

Transform4f get_system_from_file(char *filename, float lambda, int degree, float distance) {
    std::ifstream infile(filename);
    std::string line;

    Transform4f system;
    while (std::getline(infile, line)) {
        std::istringstream ls(line);
        std::string op;
        ls >> op;

        if (op == "two_plane") {
            system = two_plane_5(distance, degree);
            //cout << "two_plane" << " " << d << endl;
        }
        else if (op == "cylindrical_x") {
            float radius;
            std::string glassName1;
            std::string glassName2;
            ls >> radius;
            ls >> glassName1;
            ls >> glassName2;
            float n1 = 1.0f;
            float n2 = 1.0f;
            if (glassName1[0] >= '0' && glassName1[0] <= '9') {
                n1 = atof(glassName1.c_str());

            } else {
                OpticalMaterial glass1(glassName1.c_str());
                n1 = glass1.get_index(lambda);
            }

            if (glassName2[0] >= '0' && glassName2[0] <= '9') {
                n2 = atof(glassName2.c_str());
            } else {
                OpticalMaterial glass2(glassName2.c_str());
                n2 = glass2.get_index(lambda);
            }

            system = system >> refract_cylindrical_x_5(radius, n1, n2);
        }
        else if (op == "cylindrical_y") {
            float radius;
            std::string glassName1;
            std::string glassName2;
            ls >> radius;
            ls >> glassName1;
            ls >> glassName2;
            float n1 = 1.0f;
            float n2 = 1.0f;
            if (glassName1[0] >= '0' && glassName1[0] <= '9') {
                n1 = atof(glassName1.c_str());

            } else {
                OpticalMaterial glass1(glassName1.c_str());
                n1 = glass1.get_index(lambda);
            }

            if (glassName2[0] >= '0' && glassName2[0] <= '9') {
                n2 = atof(glassName2.c_str());
            } else {
                OpticalMaterial glass2(glassName2.c_str());
                n2 = glass2.get_index(lambda);
            }

            system = system >> refract_cylindrical_y_5(radius, n1, n2);
        }
        else if (op == "reflect_spherical") {
            float radius;
            ls >> radius;
            system = system >> reflect_spherical_5(radius, degree);
        }
        else if (op == "refract_spherical") {
                float radius;
                std::string glassName1;
                std::string glassName2;
                ls >> radius;
                ls >> glassName1;
                ls >> glassName2;

                float n1 = 1.0f;
                float n2 = 1.0f;

                if (glassName1[0] >= '0' && glassName1[0] <= '9') {
                    n1 = atof(glassName1.c_str());

                } else {
                    OpticalMaterial glass1(glassName1.c_str());
                    n1 = glass1.get_index(lambda);
                }

                if (glassName2[0] >= '0' && glassName2[0] <= '9') {
                    n2 = atof(glassName2.c_str());
                } else {
                    OpticalMaterial glass2(glassName2.c_str());
                    n2 = glass2.get_index(lambda);
                }

                system = system >> refract_spherical_5(radius, n1, n2, degree);
                //system = system >> refract_cylindrical_x_5(radius, n1, n2, degree);
                //cout << "refract_spherical" << " " << radius << " " << n1 << " " << n2 << endl;

        } else if (op == "propagate") {
            float d;
            ls >> d;
            system = system >> propagate_5(d, degree);
            //cout << "propagate" << " " << d << endl;

        } else {
            cout << "invalid op: " << op << endl;
        }
    }
    return system;
}


int main(int argc, char *argv[]) {

	auto start = std::chrono::steady_clock::now();

    if (argc < 2) {
        showUsage(argv[0]);
        // printf("usage: %s <infile.pfm> <outfile.pfm> [exposure :1.0] [degree :3] [sample-mul :1000] [entrance :19.5] [defocus :0.0] [lambda-count :12] [filter-size :1]\n", argv[0]);
        exit(1);
    }

    FILE *fp = std::fopen(argv[1], "r");
    if (!fp) {
        cerr << argv[1] << " not found";
        exit(1);
    } else {
        std::fclose(fp);
    }

    int blade_count = 0;
	int anamorphic = 1;
    int degree = 3;
    float distance = 5000000;
    float sample_mul = 1000;
    float r_entrance = 19.5;
    float defocus = 0.0;
    int num_lambdas = 12;
    int filter_size = 1;
    float exposure = 1.0;

    const char *out_file = "out.exr";

    char *system_definition_file = 0;

    char tmp;

    if (argc >= 4) {
        while ((tmp = getopt(argc - 1, &argv[1], "a:b:c:d:e:f:o:p:s:x:i:z:")) != -1) {
            switch (tmp) {

				case 'a'://system definition
                    anamorphic = atoi(optarg);
                    break;

                case 'b'://system definition
                    blade_count = atoi(optarg);
                    break;

                case 'c'://curve(degree)
                    degree = atol(optarg);
                    break;

                case 'd'://defocus
                    defocus = atof(optarg);
                    break;

                case 'e'://defocus
                    r_entrance = atof(optarg);
                    break;

                case 'f'://filter
                    filter_size = atol(optarg);
                    break;

                case 'o'://out render file
                    out_file = optarg;
                    break;

                case 'p'://lambda passes
                    num_lambdas = atol(optarg);
                    break;

                case 's'://samples(multiplier)
                    sample_mul = atof(optarg);
                    break;

                case 'x'://exposure
                    exposure = atof(optarg);
                    break;

                case 'i'://system definition
                    system_definition_file = strdup(optarg);
                    break;

                case 'z': // distance
                    distance = atof(optarg);
                    break;

                default:
                    showUsage(argv[0]);
                    break;
            }
        }
    }

    const float lambda_from = 440;
    const float lambda_to = 660;

    // cout << "-- config: " << endl;
    // cout << "exposure: " << exposure << endl;
    // cout << "sample-mul: " << sample_mul << endl;
    // cout << "entrance: " << r_entrance << endl;
    // cout << "defocus:" << defocus << endl;
    // cout << "lambda-count: " << num_lambdas << endl;
    cout << "lambda-from: " << lambda_from << endl;
    cout << "lambda-to: " << lambda_to << endl;
	cout << "filter-size: " << filter_size << endl;

    // Sensor scaling
    const float sensor_width = 80;
    const int sensor_xres = 1920;
    const int sensor_yres = 1080;
    const float sensor_scaling = sensor_xres / sensor_width;
    cout << "sensor scaling: " << sensor_scaling << endl;


    CImg<float> img_in(argv[1]);
    int width = img_in.width();
    int height = img_in.height();

    cout << "opened image with dimensions: " << width << "*" << height << endl;

    if (system_definition_file) {
        cout << "system definition file: " << system_definition_file << endl;
    }


    // cout << "degree: " << degree << endl;
    float r_pupil = r_entrance;

    // Focus on 550nm
    Transform4f system = system_definition_file ? get_system_from_file(system_definition_file, 550, degree, distance)
                                                : get_system(550, degree);


    // Determine back focal length from degree-1 terms (matrix optics)
    float d3 = find_focus_X(system);
    cout << "system focus: " << d3 << endl;
    cout << "effective focus: " << (d3 + defocus) << endl;
    // Compute magnification and output equation system
    float magnification = get_magnification_X(system >> propagate_5(d3));
    cout << "magnification: " << magnification << endl;
    //cout << "System: " << system << endl<<endl;

    // Add that propagation, plus a little animated defocus to the overall system;
    Transform4f prop = propagate_5(d3 + defocus, degree);
    system = system >> prop;

    CImg<float> img_out(sensor_xres, sensor_yres, 1, 3, 0);

    // Precompute spectrum
    float rgb[3 * num_lambdas];
    for (int ll = 0; ll < num_lambdas; ++ll) {
        float lambda = lambda_from + (lambda_to - lambda_from) * (ll / (float) (num_lambdas - 1));
        if (num_lambdas == 1) lambda = 550;
        // Convert wavelength to spectral power
        spectrum_p_to_rgb(lambda, 1, rgb + 3 * ll);
    }

    // Sample optical system at two spectral locations
    Transform4d system_spectral_center =
            (system_definition_file ? get_system_from_file(system_definition_file, 500, degree, distance) : get_system(500, degree))
                    >> prop;
    Transform4d system_spectral_right =
            (system_definition_file ? get_system_from_file(system_definition_file, 600, degree, distance) : get_system(600,
                                                                                                             degree))
                    >> prop;

    // Obtain (xyworld + xyaperture + lambda) -> (ray) mapping including chromatic effects,
    // by linear interpolation of the two sample systems drawn above
    System54f system_spectral = system_spectral_center.lerp_with(system_spectral_right, 550, 600);

    // dx and dy after propagation are really only needed for Lambertian
    // term; hence: combine them to obtain sin^2 = 1 - cosine^2 term in equation 2:
    system_spectral[2] = (system_spectral[2] * system_spectral[2] + system_spectral[3] * system_spectral[3]);
    system_spectral[2] %= 2;
    System53d system_lambert_cos2 = system_spectral.drop_equation(3);

    // Support of an input image pixel in world plane
    float pixel_size = sensor_width / (float) width / magnification;



    float blade_positions[blade_count * 2];

    for (int b = 0; b < blade_count; ++b) {
        blade_positions[b * 2] = std::cos(b * 3.141596535 * 2.0 / blade_count) * r_pupil;
        blade_positions[b * 2 + 1] = std::sin(b * 3.141596535 * 2.0 / blade_count) * r_pupil;
        cout << blade_positions[b * 2] << " " << blade_positions[b * 2 + 1] << endl;
    }


    for (int ll = 0; ll < num_lambdas; ++ll) {
        float lambda = lambda_from + (lambda_to - lambda_from) * (ll / (float) (num_lambdas - 1));
        if (num_lambdas == 1) lambda = 550;
        cout << "[" << lambda << "nm]" << flush;

        // Bake lambda dependency
        System43f system_lambda = system_lambert_cos2.bake_input_variable(4, lambda);
        system_lambda %= degree;

        // Bake lambda into derivatives as well:

        for (int j = 0; j < height; j++) {
            if (!(j % 10)) cout << "." << flush;
            const float y_sensor = ((j - height / 2) / (float) width) * sensor_width;
            const float y_world = y_sensor / magnification;

            // Bake y dependency
            System33f system_y = system_lambda.bake_input_variable(1, y_world);

            for (int i = 0; i < width; i++) {
                const float x_sensor = (i / (float) width - 0.5) * sensor_width;
                const float x_world = x_sensor / magnification;

                // Sample intensity at wavelength lambda from source image
                const float rgbin[3] = {
                        exposure * img_in.linear_atXY(i, j, 0, 0, 0),
                        exposure * img_in.linear_atXY(i, j, 0, 1, 0),
                        exposure * img_in.linear_atXY(i, j, 0, 2, 0)};
                float L_in = spectrum_rgb_to_p(lambda, rgbin);

                // Quasi-importance sampling:
                // pick number of samples according to pixel intensity
                int num_samples = max(1, (int) (L_in * sample_mul));

                float sample_weight = L_in / num_samples;

                // With that, we can now start sampling the aperture:
                for (int sample = 0; sample < num_samples; ++sample) {

                    // Rejection-sample points from lens aperture:
                    float x_ap, y_ap;

                    if (blade_count == 0) {
                        do {
                            x_ap = (rand() / (float) RAND_MAX - 0.5) * 2 * r_pupil;
                            y_ap = (rand() / (float) RAND_MAX - 0.5) * 2 * r_pupil;
                        } while (x_ap * x_ap + y_ap * y_ap > r_pupil * r_pupil);
                    }
                    else {
                        bool inside;
                        do {
                            inside = true;
                            x_ap = (rand() / (float) RAND_MAX - 0.5) * 2 * r_pupil;
                            y_ap = (rand() / (float) RAND_MAX - 0.5) * 2 * r_pupil;

                            for (int b = 0; b < blade_count; ++b) {
                                float bx = blade_positions[((b + 1) % blade_count) * 2] - blade_positions[b * 2];
                                float by =
                                        blade_positions[((b + 1) % blade_count) * 2 + 1] - blade_positions[b * 2 + 1];

                                float px = x_ap - blade_positions[b * 2];
                                float py = y_ap - blade_positions[b * 2 + 1];


                                float det = (px * by) - (py * bx); //bx * py - px * by;
                                if (det > 0) {

                                    inside = false;
                                    break;
                                }
                            }


                        } while (!inside);
                    }


                    float in[5], out[4];

                    // Fill in variables and evaluate systems:
                    in[0] = x_world + pixel_size * (rand() / (float) RAND_MAX - 0.5);
                    in[1] = x_ap / anamorphic;
                    in[2] = y_ap;

                    system_y.evaluate(in, out);

                    // Scale to pixel size:
                    out[0] = out[0] * sensor_scaling + sensor_xres / 2;
                    out[1] = out[1] * sensor_scaling + sensor_yres / 2;

                    // out[2] contains one minus square of Lambertian cosine
                    float lambert = sqrt(1 - out[2]);
                    if (lambert != lambert) lambert = 0; // NaN check

                    img_out.set_linear_atXY(lambert * sample_weight * rgb[0 + 3 * ll], out[0], out[1], 0, 0, true);
                    img_out.set_linear_atXY(lambert * sample_weight * rgb[1 + 3 * ll], out[0], out[1], 0, 1, true);
                    img_out.set_linear_atXY(lambert * sample_weight * rgb[2 + 3 * ll], out[0], out[1], 0, 2, true);

                }
            }
        }
        cout << endl;
    }


    // Fix gamut problem (pure wavelengths sometimes result in negative RGB)
    for (int j = 0; j < sensor_yres; ++j) {
        for (int i = 0; i < sensor_xres; ++i) {
            float max_value = max(img_out.atXY(i, j, 0, 0),
                                  max(img_out.atXY(i, j, 0, 1),
                                      img_out.atXY(i, j, 0, 2)));
            img_out.atXY(i, j, 0, 0) = max(img_out.atXY(i, j, 0, 0), 0.02f * max_value);
            img_out.atXY(i, j, 0, 1) = max(img_out.atXY(i, j, 0, 1), 0.02f * max_value);
            img_out.atXY(i, j, 0, 2) = max(img_out.atXY(i, j, 0, 2), 0.02f * max_value);

        }
    }

    img_out.save(out_file);
	auto end = chrono::steady_clock::now();
	double elapsed_time_ns = double(chrono::duration_cast <chrono::nanoseconds> (end - start).count());
	cout << "" << out_file << " — " << elapsed_time_ns / 1e9 << "s elapsed." << endl;
}


void showUsage(char *s) {
    cout << "" << endl;
    cout << "DYLE — Polynomial Optics" << endl;
    cout << "Usage:   " << s << " sourcefile outputfile [-option] [argument]" << endl;
    cout << "         " << "-a (fake)anamorphic" << endl;
    cout << "         " << "-b iris blade count" << endl;
    cout << "         " << "-c curve/degree" << endl;
    cout << "         " << "-d defocus" << endl;
    cout << "         " << "-e entrance" << endl;
    cout << "         " << "-f filter" << endl;
    cout << "         " << "-o output file" << endl;
    cout << "         " << "-p lambda passes" << endl;
    cout << "         " << "-s samples" << endl;
    cout << "         " << "-x exposure" << endl;
    cout << "         " << "-i .lens file" << endl;
    cout << "" << endl;
}
