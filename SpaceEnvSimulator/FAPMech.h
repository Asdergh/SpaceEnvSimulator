#pragma once

#include<vector>
#include<iostream>
#include<string>
#include<math.h >
#include<iostream>
#include<fstream>


#define PI  3.14159265358979323846
#define KA_mass 1250.0
#define KA_SA 1.5 * 1.8
#define Cx_coeff 2.0 / 2.5
#define Sigma_coeff (Cx_coeff * KA_SA) / (2 * KA_mass)

#define H_peregee 350.0 + 6371.0
#define H_apogee 450.0 + 6371.0
#define I_param (20.0 * PI) / 180.0

#define Semimajor_axis (H_peregee + H_apogee) / 2
#define Ext_param (H_apogee - H_peregee) / (H_apogee + H_peregee)
#define Focal_param Semimajor_axis * (1 - pow(Ext_param, 2.0))
#define Venus_grav_param 324859.0
#define Earth_grav_param 398600.4415
#define Venus_anguler_vel 6.52 * pow(10, -5)

#define Omega_param (10.0 * PI) / 180.0
#define omega_param 0.0
#define M_param (45.0 * PI) / 180.0

#define Earth_anguler_vel 7.2921159 * pow(10, -5)
#define Static_density 1.58868 * pow(10, -8)
#define Accuracy (0.001 * 3.14) / 180.0

#define Ext_param3 0.0067385254
#define ND_upper_120 1.58868 * pow(10, -8)

#define SM_Earth  6378136.0


class FAPMech
{
private:


	
	float* position_AGSK = new float[3];
	float* position_GSK = new float[3];
	float* position_LBH = new float[3];

	float transversial_velocity;
	float radial_velocity;
	float velocity_normal;
	float* velocity_vector = new float[3];
	double** indignant_acceleration_tensor = new double* [7];


	float r_normal_AGSK;
	float theta_param;
	float ext_param;
	float ext_past_param;
	float u_param;
	float time;

	float N_param;
	float B_param;
	float H_param;
	float La_param;
	float L_param;
	float D_param;

	double* density_night_param = new double[7];
	double* density_param = new double[7];


	double* a0_first_level = new double[7]{26.8629, 27.4598, 28.6395, 29.6418, 30.1671, 29.7578, 30.7854};
	double* a1_first_level = new double[7]{ -0.451674, -0.463668, -0.490987, -05.149557, -0.527837, -0.517915, -0.545695 };
	double* a2_first_level = new double[7]{ 0.002903397, 0.002974, 0.00320649, 0.0034126, 0.00353211, 0.00342699, 0.00370328 };
	double* a3_first_level = new double[7]{ -1.06953 * pow(10.0, -5.0), -1.0753 * pow(10.0, -5.0), -1.1681 * pow(10.0, -5.0), -1.25785 * pow(10.0, -5.0), -1.30227 * pow(10.0, -5.0), -1.24137 * pow(10.0, -5.0), -1.37072 * pow(10.0, -5.0) };
	double* a4_first_level = new double[7]{ 2.21598 * pow(10.0, -8.0), 2.17059 * pow(10.0, -8.0), 2.36847 * pow(10.0, -8.0), 2.5727 * pow(10.0, -8.0), 2.66455 * pow(10.0, -8.0), 2.48209 * pow(10.0, -8.0), 2.80614 * pow(10.0, -8.0) };
	double* a5_first_level = new double[7] { -2.42941 * pow(10.0, -11.0), -230249 * pow(10.0, -11.0), -2.51809 * pow(10.0, -11.0) * pow(10.0, -11.0), -285432 * pow(10.0, -11.0) * pow(10.0, -11.0), -2.58413 * pow(10.0, -11.0), -3.00184 * pow(10.0, -11.0) };
	double* a6_first_level = new double[7]{ 1.09926 * pow(10.0, -14.0), 1.00123 * pow(10.0, -14.0), 1.09536 * pow(10.0, -14.0), 1.21091 * pow(10.0, -14.0), 1.25009 * pow(10.0, -14.0), 1.09383 * pow(10.0, -14.0), 1.31142 * pow(10.0, -14.0) };

	double* a0_second_level = new double[7]{ 17.8781, -2.54909, -13.9599, -23.3079, -14.7264, -4.912, -5.40952 };
	double* a1_second_level = new double[7]{ -0.132025, 0.0140064, 0.0844951, 0.135141, 0.0713256, 0.0108326, 0.00550749 };
	double* a2_second_level = new double[7]{ 0.000227717, -0.00016946, -0.000328875, -0.000420802, -0.000228015, -8.10546e-5, -3.78851e-5 };
	double* a3_second_level = new double[7]{ -2.2543 * pow(10.0, -7), 3.27196 * pow(10.0, -7.0), 5.05918 * pow(10.0, -7.0), 5.73717 * pow(10.0, -7.0), 2.8487 * pow(10.0, -7.0), 1.15712 * pow(10.0, -7.0), 2.4808 * pow(10.0, -8.0) };
	double* a4_second_level = new double[7]{ 1.33574 * pow(10.0, -10.0), -2.8763 * pow(10.0, -10.0), -3.92299 * pow(10.0, -10.0), -4.03238 * pow(10.0, -10.0), -1.74282 * pow(10.0, -10.0), -813296 * pow(10.0, -11.0), 4.92183 * pow(10.0, -12.0) };
	double* a5_second_level = new double[7]{ -4.50458 * pow(10.0, -14.0), 1.22625 * pow(10.0, -13.0), 1.52279 * pow(10.0, -13.0), 1.42846 * pow(10.0, -13.0), 5.08071 * pow(10.0, -14.0), 3.04913 * pow(10.0, -14.0), -8.45011 * pow(10.0, -18.0) };
	double* a6_second_level = new double[7]{ 6.72086 * pow(10.0, -18.0), -2.05736 * pow(10.0, -17.0), -2.35576 * pow(10.0, -17.0), -2.01726 * pow(10.0, -17.0), -5.34955 * pow(10.0, -18.0), -4.94989 * pow(10.0, -18.0), 1.9849 * pow(10.0, -18.0) };

	float max_time = 2 * 3.14 * sqrt(pow(Semimajor_axis, 3) / Venus_grav_param);

public:

	FAPMech();
	~FAPMech();
	void calculate_anomalies();
	void calculate_AGSK_pos();
	void calculate_velocities();
	void calculate_GSK_pos();
	float calculate_LBH_params();
	void calculate_density();
	void calculate_acceleration();
	void run_simulation();


};

