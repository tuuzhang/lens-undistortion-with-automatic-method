/*
   Copyright (c) 2010-2014, AMI RESEARCH GROUP <lalvarez@dis.ulpgc.es>
   License : CC Creative Commons "Attribution-NonCommercial-ShareAlike"
   see http://creativecommons.org/licenses/by-nc-sa/3.0/es/deed.en
   */


/**
 * @file lens_distortion_correction_2p_iterative_optimization.cpp
 * @brief Distortion correction using two-parameter polynomial and division lens
 * distortion models, including an iterative optimization process and ensuring
 * the invertibility of the models
 *
 * @author Luis Alvarez <lalvarez@dis.ulpgc.es> and Daniel Santana-Cedrés <dsantana@ctim.es>
 */


//Included libraries
#include "ami_image.h"
#include "filters.h"
#include "image_contours.h"
#include "line_extraction.h"
#include "ami_image_primitives.h"
#include "lens_distortion_procedures.h"
#include "ami_utilities.h"
#include <iostream>
#include <fstream>
#include <string>
using namespace std;

static const int neighborhood_radius = 2; /*radius of neighborhood to take into
										  account*/
static const int min_neighbor_points = 2; /*min number of contour points in a
										  neighborhood*/
static const double min_orientation_value = 0.95; /*  min average scalar product
											 of neigborhood point orientation */
static const int min_distance_point = 1; //minimum distance between contour points
static const int max_lines = 100; //maximun number of lines estimated
static const float angle_resolution = 0.1; //angle discretization step(in degrees)
static const float distance_resolution = 1.; //line distance discretization step
static const float distortion_parameter_resolution = 0.1;/*distortion parameter
														 discretization step*/
static const ImageAmplification amp = FIT_WIDTH; /*Image amplification factor
											  for undistorting the input image*/

//------------------------------------------------------------------------------
/** \fn double energy_minimization(lens_distortion_model &ldm,
 *                                 image_primitives &ip, int w, int h,
 *                                 bool opt_center)
 * \brief Function to compute the energy minimization according to the type of
 the lens distortion model and the center optimization
 * \param [in,out] ldm Lens distortion model
 * \param [in] ip Set of primitives
 * \param [in] w Image width
 * \param [in] h Image height
 * \param [in] opt_center Option for center optimization
 * \return Returns the obtained value for the minimized energy
 * \author Luis Alvarez and Daniel Santana-Cedrés
 */
double energy_minimization(lens_distortion_model &ldm, const image_primitives &ip,
	int w, int h, bool opt_center)
{
	vector<bool> vtf(4); vtf[0] = vtf[1] = true; vtf[2] = vtf[3] = false;
	vector<bool> vtt(4, true);
	double error = 0., tt_error = 0.;
	lens_distortion_model tf_ldm, tt_ldm;
	vector<line_points> l = ip.get_lines();
	error = model_center_estimation_2p(l, ldm, w, h, vtf);
	tf_ldm = ldm;
	if (opt_center)
	{
		tt_error = model_center_estimation_2p(l, ldm, w, h, vtt);
		tt_ldm = ldm;
	}

	if ((fabs(tf_ldm.get_distortion_center().x -
		tt_ldm.get_distortion_center().x) < 0.2*w) ||
		(fabs(tf_ldm.get_distortion_center().y -
		tt_ldm.get_distortion_center().y) < 0.2*h) ||
		check_invertibility(tt_ldm, w, h))
	{
		ldm = tt_ldm;
		error = tt_error;
	}

	return error;
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
/** \fn double iterative_optimization(
								  ami::image_contours &contours,
								  image_primitives &i_primitives,
								  float distance_point_line_max_hough,
								  int max_lines,
								  float angle_resolution,
								  float distance_resolution,
								  float distortion_parameter_resolution,
								  float angle_point_orientation_max_difference,
								  bool opt_center,
								  int width, int height
								  )
								  * \brief Function to compute the iterative optimization
								  * \param [in] contours image_contours object with the edges and its
								  orientations
								  * \param [in,out] i_primitives Set of primitives detected in the previous step
								  * \param [in] distance_point_line_max_hough Maximum distance between the point
								  and the line for which votes
								  * \param [in] max_lines Maximum number of lines to detect
								  * \param [in] angle_resolution Angle resolution for the Hough transform
								  * \param [in] distance_resolution Distance resolution for the Hough transform
								  * \param [in] distortion_parameter_resolution Resolution for the normalized
								  distortion parameter
								  * \param [in] angle_point_orientation_max_difference Maximum difference in
								  degrees between the orientation of the point and the orientation
								  of the line
								  * \param [in] opt_center Option for the center optimization
								  * \return Returns the final error obtained
								  * \author Luis Alvarez and Daniel Santana-Cedrés
								  */
double iterative_optimization(
	const ami::image_contours &contours,
	image_primitives &i_primitives,
	float distance_point_line_max_hough,
	int max_lines,
	float angle_resolution,
	float distance_resolution,
	float distortion_parameter_resolution,
	float angle_point_orientation_max_difference,
	bool opt_center
	)
{
	double final_error = 0.;
	int width = contours.get_width();
	int height = contours.get_height();
	//We initialize the previous model, the best one and the previous set of 
	//primitives
	i_primitives.get_distortion().get_d().resize(3);
	i_primitives.get_distortion().get_d()[2] = 0.;
	lens_distortion_model previous_model = i_primitives.get_distortion();
	lens_distortion_model best_model = previous_model;
	image_primitives previous_ip = i_primitives;
	//Tolerance for convergence
	double TOL = 1e-2;
	//Fails counter
	int fail_count = 0;
	//Number of points: current, best and next
	int num_points = count_points(i_primitives);
	int best_num_points = num_points;
	int next_num_points = num_points;
	//We apply the process until the number of points is not significantly greater
	//or until the process fails three times 
	while ((next_num_points >= (num_points*(1 + TOL))) || (fail_count < 3))
	{
		double error = energy_minimization(previous_model, i_primitives, width,
			height, opt_center);
		i_primitives.clear();
		//CALL TO IMPROVED HOUGH WITH THE MODEL COMPUTED BEFORE
		line_equation_distortion_extraction_improved_hough(
			contours, i_primitives, distance_point_line_max_hough,
			max_lines, angle_resolution,
			distance_resolution, 0., 0.,
			distortion_parameter_resolution,
			angle_point_orientation_max_difference,
			true, previous_model);

		int local_num_points = count_points(i_primitives);
		if (local_num_points > next_num_points)
		{
			//We update the primitives only if the result is better
			if (local_num_points > best_num_points)
			{
				previous_ip = i_primitives;
				best_num_points = local_num_points;
				best_model = previous_model;
				final_error = error;
			}
		}
		else
		{
			fail_count++;
		}
		num_points = next_num_points;
		next_num_points = local_num_points;
	}
	//We take the last and best image primitives object and model
	i_primitives = previous_ip;
	i_primitives.set_distortion(best_model);

	//We return the average error
	return (final_error / (double)best_num_points);
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//Auxiliary procedure for filling the image with the result of the edge
//detection
void fill_edges3c(const ami::image_contours &contours,
	ami::image<unsigned char> &edges3c,
	const string image_name)
{
	int size_ = contours.get_width() * contours.get_height();
	const vector<int> &index = contours.get_index();
	for (int i = 0; i < (int)index.size(); i++)
	{
		edges3c[index[i]] = 0;
		edges3c[index[i] + size_] = 0;
		edges3c[index[i] + 2 * size_] = 0;
	}
	//Writing Canny detector output after the cleaning process
	edges3c.write(image_name);
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//Procedure for writing in a file the lens distortion model and the lines and 
//associated points
void write_output_file(const input_params &i_params,
	const image_primitives &i_primitives,
	int width, int height, double final_error)
{
	i_primitives.write(i_params.primitives_file.c_str());
	// writing function parameters and basic outputs :
	ofstream fs("output.txt");//Output file
	fs << "Selected parameters:" << endl;
	fs << "\t High Canny's threshold: " << i_params.canny_high_threshold << endl;
	fs << "\t Initial normalized distortion parameter: "
		<< i_params.initial_distortion_parameter << endl;
	fs << "\t Final normalized distortion parameter: "
		<< i_params.final_distortion_parameter << endl;
	fs << "\t Maximum distance between points and line: "
		<< i_params.distance_point_line_max_hough << endl;
	fs << "\t Maximum differente between edge point and line orientations: "
		<< i_params.angle_point_orientation_max_difference << endl;
	fs << "\t Model applied: " << i_params.tmodel << endl;
	string copt = (i_params.opt_center) ? "True" : "False";
	fs << "\t Center optimization: " << copt << endl;
	fs << "-------------------------" << endl;
	fs << "Results: " << endl;
	fs << "\t Number of detected lines: " << i_primitives.get_lines().size()
		<< endl;
	int count = count_points(i_primitives);
	fs << "\t Total amount of line points: " << count << endl;
	fs << "\t Distortion center: (" << i_primitives.get_distorsion_center().x <<
		", " << i_primitives.get_distorsion_center().y << ")" << endl;
	double p1 = 0., p2 = 0.;
	compute_ps(p1, p2, i_primitives.get_distortion(), width, height);
	fs << "\t Estimated normalized distortion parameters: p1 = " << p1
		<< " p2 = " << p2 << endl;
	fs << "\t Average squared error distance in pixels between line and"
		<< " associated points = " << final_error << endl;
	fs.close();
}
//------------------------------------------------------------------------------

int set_params_lens_distortion_correction_2p_iterative_optimization(input_params &i_params)
{
	//Messages
	string error_message_head("The argument ");
	string error_message_image(" must be a png image.\n");
	string error_message_numgt(" must be greater than ");
	string error_message_numgte(" must be greater or equal to ");
	string error_message("");

	//Input image
	i_params.input_name = "example/building.png";

	//Canny image
	i_params.canny_name = "example/building_canny.png";

	//Hough image
	i_params.hough_name = "example/building_hough.png";

	//Corrected image
	i_params.result_name = "example/building_corrected_image.png";

	//High threshold Canny
	i_params.canny_high_threshold = 0.8;
	if (i_params.canny_high_threshold <= 0.7)
		error_message += error_message_head + "5 (high threshold Canny)" +
		error_message_numgt + "0.7\n";

	//Initial distortion parameter
	i_params.initial_distortion_parameter = 0.0;
	if (i_params.initial_distortion_parameter < -0.5)
		error_message += error_message_head + "6 (initial distortion parameter)" +
		error_message_numgte + "-0.5\n";

	//Final distortion parameter
	i_params.final_distortion_parameter = 3.0;
	if (i_params.final_distortion_parameter < i_params.initial_distortion_parameter)
		error_message += error_message_head + "7 (final distortion parameter)" +
		error_message_numgte + to_string(i_params.initial_distortion_parameter) + "\n";

	//Maximum distance between points and line
	i_params.distance_point_line_max_hough = 3.0;
	if (i_params.distance_point_line_max_hough < 0.0)
		error_message += error_message_head +
		"8 (maximum distance between points and line)" +
		error_message_numgte + "0.0\n";

	//Maximum difference between point angle and line angle
	i_params.angle_point_orientation_max_difference = 10.0;
	if (i_params.angle_point_orientation_max_difference < 0.0 ||
		i_params.angle_point_orientation_max_difference > 45.0)
		error_message += error_message_head + "9 (maximum difference between point" +
		" angle and line angle) must be between 0.0 and 45.0\n";

	//Type of the lens distortion model, pol or div
	i_params.tmodel = "pol";

	//Optimization of the center of the lens distortion model
	i_params.opt_center = true;

	//Name for the primitives file
	i_params.primitives_file = "example/primitives.txt";

	if (error_message.length() > 0)
	{
		cout << error_message;
		return -1;
	}

	return 0;
}

//------------------------------------------------------------------------------
int main(int argc, char *argv[])
{
	//if (argc != 13){
	//	print_function_syntax_lens_distortion_correction_2p_iterative_optimization();
	//	exit(EXIT_FAILURE);
	//}

	//input_params i_params;

	//if (check_params_lens_distortion_correction_2p_iterative_optimization(argv,
	//	i_params)
	//	!= 0){
	//	manage_failure(argv, 0);
	//	exit(EXIT_SUCCESS);
	//}

	// set the parameters
	input_params i_params;
	if (set_params_lens_distortion_correction_2p_iterative_optimization(i_params) != 0)
	{
		exit(EXIT_FAILURE);
	}
	

	//We read the input image and some variables are initialized
	ami::image<unsigned char> input(i_params.input_name); // input image
	int width = input.width(), height = input.height();//input image dimensions
	int size_ = width*height; // image size
	ami::image<unsigned char> gray(width, height, 1, 0);//gray-level image to call 
	//canny
	ami::image<unsigned char> edges(width, height, 1, 0);//image to store edge 
	//information
	bool print_msg = true;

	//Converting the input image to gray level
	for (int i = 0; i < size_; i++)
		gray[i] = 0.3*input[i] + 0.59*input[i + size_] + 0.11*input[i + size_ * 2];

	//ALGORITHM STAGE 1 : Detecting edges with Canny
	if (print_msg) cout << "Detecting edges with Canny..." << endl;
	i_params.canny_low_threshold = 0.7; //default value for canny lower threshold
	//ami::image_contours contours=canny(gray,edges,
	//                                            i_params.canny_low_threshold,
	//                                            i_params.canny_high_threshold);

	ami::image_contours contours(gray.width(), gray.height());
	canny_with_contours(
		contours, gray, edges, 
		i_params.canny_low_threshold, 
		i_params.canny_high_threshold);

	//We create writable 3 channel images for edges and gray level
	ami::image<unsigned char> edges3c(width, height, 3, 255);

	//We clean the contours
	contours.clean(
		neighborhood_radius,
		min_neighbor_points,
		min_orientation_value,
		min_distance_point
		);

	fill_edges3c(contours, edges3c, i_params.canny_name);
	if (print_msg) cout << "...edges detected" << endl;

	//ALGORITHM STAGE 2 : Detecting lines with improved_hough_quotient
	if (print_msg) cout << "Detecting lines with improved Hough and "
		<< i_params.tmodel << " model..." << endl;
	image_primitives i_primitives;//object to store output edge line structure
	lens_distortion_model ini_ldm;
	if (i_params.tmodel == string("pol"))
		ini_ldm.set_type(POLYNOMIAL);
	else
		ini_ldm.set_type(DIVISION);
	//we call 3D Hough line extraction
	line_equation_distortion_extraction_improved_hough(
		contours,
		i_primitives,
		i_params.distance_point_line_max_hough,
		max_lines,
		angle_resolution,
		distance_resolution,
		i_params.initial_distortion_parameter,
		i_params.final_distortion_parameter,
		distortion_parameter_resolution,
		i_params.angle_point_orientation_max_difference,
		true,
		ini_ldm
		);

	//ALGORITHM STAGE 3 : We apply the iterative optimization process
	double final_error = iterative_optimization(
		contours,
		i_primitives,
		i_params.distance_point_line_max_hough,
		max_lines,
		angle_resolution,
		distance_resolution,
		distortion_parameter_resolution,
		i_params.angle_point_orientation_max_difference,
		i_params.opt_center);

	//We check if the iterative optimization process finishes properly
	if (i_primitives.get_lines().size() == 0){
		manage_failure(argv, 0);
		exit(EXIT_SUCCESS);
	}

	if (print_msg) cout << "...lines detected: " << i_primitives.get_lines().size()
		<< " with " << count_points(i_primitives) << " points"
		<< endl;

	//Drawing the detected lines on the original image to illustrate the results
	drawHoughLines(i_primitives, edges3c);
	edges3c.write(i_params.hough_name);

	//ALGORITHM STAGE 4: Correcting the image distortion using the estimated model
	if (i_primitives.get_distortion().get_d().size() > 0){
		if (print_msg) cout << "Correcting the distortion..." << endl;
		ami::image<unsigned char> undistorted(width, height, 3, 0);
		undistorted = undistort_image_inverse(
			input, // input image
			i_primitives.get_distortion(), // lens distortion model
			amp // value to fix the way the corrected image is scaled to fit input
			//size image
			);
		//Writing the distortion corrected image
		undistorted.write(i_params.result_name);
		if (print_msg) cout << "...distortion corrected." << endl;
	}

	// WRITING OUTPUT TEXT DOCUMENTS
	write_output_file(i_params, i_primitives, width, height, final_error);
	exit(EXIT_SUCCESS);
}