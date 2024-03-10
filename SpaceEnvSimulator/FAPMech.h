#pragma once

#include<vector>
#include<iostream>
#include<string>
#include<math.h >


#define H_peregee 450.0
#define H_apogee 340.0
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

class FAPMech
{
public:


	std::vector<float> position_AGSK;
	std::vector<float> position_GSK;
	std::vector<float> position_LBH;

	float transversial_velocity;
	float radial_velocity;
	float velocity_normal;
	std::vector<float> velocity_vector;

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

	float density_param;
	float density_night_param;


public:


	void calculate_anomalies();
	void calculate_AGSK_pos();
	void calculate_velocityes();
	void calculate_GSK_pos();
	float calculate_LBH_params();
	void run_simulation();

};

