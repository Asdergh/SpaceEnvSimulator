#pragma once

#include<vector>
#include<iostream>
#include<string>
#include<math.h >
#include<iostream>
#include<fstream>


#define KA_mass 1250.0
#define KA_SA 1.5 * 1.8
#define Cx_coeff 2.0 * 2.5
#define Sigma_coeff (Cx_coeff * KA_SA) / (2 * KA_mass)

#define H_peregee 350.0
#define H_apogee 450.0
#define I_param 20.0

#define Semimajor_axis (H_peregee + H_apogee) / 2
#define Ext_param (H_apogee - H_peregee) / (H_apogee + H_peregee)
#define Focal_param (H_apogee - H_peregee) / (2 * Semimajor_axis)
#define Venus_grav_param 324859.0
#define Venus_anguler_vel 6.52 * pow(10, -5)

#define Omega_param 10.0
#define omega_param 0.0
#define M_param 45.0

#define Earth_anguler_vel 7.2921159 * pow(10, -5)
#define Static_density 1.58868 * pow(10, -8)
#define Accuracy 0.001

#define Ext_param3 0.0067385254
#define ND_upper_120 1.58868 * pow(10, -8)

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
	float** indignant_acceleration_120_tensor = new float* [6];
	float** indignant_acceleration_500_tensor = new float* [6];

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

	float* newdensity_param_first = new float[6];
	float* density_night_param_first = new float[6];
	float* density_param_second = new float[6];
	float* density_night_param_second = new float[6];

	float* a0_first_level = new float[7]{26.8629, 27.4598, 28.6395, 29.6418, 30.1671, 29.7578, 30.7854};
	float* a1_first_level = new float[7]{ -0.451674, -0.463668, -0.490987, -05.149557, -0.527837, -0.517915, -0.545695 };
	float* a2_first_level = new float[8]{ 0.002903397, 0.002974, 0.00320649, 0.0034126, 0.00353211, 0.00342699, 0.00370328 };
	float* a3_first_level = new float[7]{ -1.06953e-5, -1.0753e-5, -1.1681e-5, -1.25785e-5, -1.30227e-5, -1.24137e-5, -1.37072e-5 };
	float* a4_first_level = new float[7]{ 2.21598e-8, 2.17059e-8, 2.36847e-8, 2.5727e-8, 2.66455e-8, 2.48209e-8, 2.80614e-8 };
	float* a5_first_level = new float[7]{ -2.42941e-11, -230249e-11, -2.51809e-11, -285432e-11, -2.58413e-11, -3.00184e-11 };
	float* a6_first_level = new float[7]{ 1.09926e-14, 1.00123e-14, 1.09536e-14, 1.21091e-14, 1.25009e-14, 1.09383e-14, 1.31142e-14 };

	float* a0_second_level = new float[7]{ 17.8781, -2.54909, -13.9599, -23.3079, -14.7264, -4.912, -5.40952 };
	float* a1_second_level = new float[7]{ -0.132025, 0.0140064, 0.0844951, 0.135141, 0.0713256, 0.0108326, 0.00550749 };
	float* a2_second_level = new float[7]{ 0.000227717, -0.00016946, -0.000328875, -0.000420802, -0.000228015, -8.10546e-5, -3.78851e-5 };
	float* a3_second_level = new float[7]{ -2.2543e-7, 3.27196e-7, 5.05918e-7, 5.73717e-7, 2.8487e-7, 1.15712e-7, 2.4808e-8 };
	float* a4_second_level = new float[7]{ 1.33574e-10, -2.8763e-10, -3.92299e-10, -4.03238e-10, -1.74282e-10, -813296e-11, 4.92183e-12 };
	float* a5_second_level = new float[7]{ -4.50458e-14, 1.22625e-13, 1.52279e-13, 1.42846e-13, 5.08071e-14, 3.04913e-14, -8.45011e-18 };
	float* a6_second_level = new float[7]{ 6.72086e-18, -2.05736e-17, -2.35576e-17, -2.01726e-17, -5.34955e-18, -4.94989e-18, 1.9849e-18 };

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

