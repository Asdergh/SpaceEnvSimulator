#include "FAPMech.h"



FAPMech::FAPMech()
{
	position_AGSK.push_back(0);
	position_AGSK.push_back(0);
	position_AGSK.push_back(0);

	position_GSK.push_back(0);
	position_GSK.push_back(0);
	position_GSK.push_back(0);

	position_LBH.push_back(0);
	position_LBH.push_back(0);
	position_LBH.push_back(0);

	velocity_vector.push_back(0);
	velocity_vector.push_back(0);
	velocity_vector.push_back(0);


	transversial_velocity = 0;
	radial_velocity = 0;
	velocity_normal = 0;

	r_normal_AGSK = 0;
	theta_param = 0;
	ext_param = 0;
	ext_past_param = 0;
	u_param = 0;
	time = 0;

	N_param = 0;
	B_param = 0;
	H_param = 0;
	La_param = 0;
	L_param = 0;
	D_param = 0;

}
void FAPMech::calculate_anomalies()
{
	ext_past_param = M_param;
	int flag = 0;
	while (abs(ext_param - ext_past_param) <= Accuracy || flag == 0)
	{

		ext_param = M_param + Ext_param * sin(ext_past_param);
		ext_past_param = ext_param;
		flag++;

		//out_data_file << "test_ext_param: " << ext_param << std::endl;
		if (flag == 10) {
			break;
		}
	}

	r_normal_AGSK = Semimajor_axis * (1 - Ext_param * cos(ext_param));
	theta_param = 2 * atan(sqrt((1 + Ext_param) / (1 - Ext_param)) * tan(ext_param / 2));
	u_param = theta_param;
}

void FAPMech::calculate_AGSK_pos()
{
	position_AGSK[0] = r_normal_AGSK * (cos(u_param) * cos(Omega_param) - sin(u_param) * sin(Omega_param) * cos(I_param));
	position_AGSK[1] = r_normal_AGSK * (cos(u_param) * sin(Omega_param) + sin(u_param) * cos(Omega_param) * cos(I_param));
	position_AGSK[2] = r_normal_AGSK * (sin(u_param) * sin(I_param));
}

void FAPMech::calculate_velocities()
{
	transversial_velocity = sqrt(Venus_grav_param / Focal_param) * Ext_param * sin(theta_param);
	radial_velocity = sqrt(Venus_grav_param / Focal_param) * (1 + Ext_param * cos(theta_param));
	velocity_normal = sqrt(Venus_grav_param / Focal_param) * sqrt(1 + 2 * Ext_param * cos(theta_param) + pow(Ext_param, 2));
	velocity_vector[0] = -transversial_velocity;
	velocity_vector[1] = -radial_velocity;
	velocity_vector[2] = 0.0;
}

void FAPMech::calculate_GSK_pos()
{
	float ven_angle = Venus_anguler_vel * time;

	position_GSK[0] = (position_AGSK[0] * cos(ven_angle) + position_AGSK[1] * sin(ven_angle));
	position_GSK[1] = (-position_AGSK[0] * sin(ven_angle) + position_AGSK[1] * cos(ven_angle));
	position_GSK[2] = position_AGSK[2];

}

float FAPMech::calculate_LBH_params()
{
	D_param = sqrt(pow(position_GSK[0], 2) + pow(position_GSK[1], 2));

	if (D_param == 0)
	{
		B_param = (3.14 / 2) * position_GSK[2] / abs(position_GSK[2]);
		L_param = 0;
		H_param = position_GSK[2] * sin(B_param) - Semimajor_axis * sqrt(1 - Ext_param3 * pow(sin(B_param), 2));
		return H_param;
	}

	else if (D_param > 0)
	{
		La_param = asin(position_GSK[1] / D_param);

		if (position_GSK[1] < 0 and position_GSK[0] > 0)
		{
			L_param = 2 * 3.14 - La_param;
		}

		else if (position_GSK[1] < 0 and position_GSK[0] < 0)
		{
			L_param = 3.14 + La_param;
		}

		else if (position_GSK[1] > 0 and position_GSK[0] < 0)
		{
			L_param = 3.14 - La_param;
		}

		else if (position_GSK[1] > 0 and position_GSK[0] > 0)
		{
			L_param = La_param;
		}

	}

	if (position_GSK[0] == 0) {
		
		B_param = 0;
		H_param = D_param - Semimajor_axis;
		return H_param;
	
	}

	else if (position_GSK[0] != 0)
	{
		float r_norma_GSK = sqrt(pow(position_GSK[0], 2) + pow(position_GSK[1], 2) + pow(position_GSK[2], 2));
		float c_param = asin(position_GSK[2] / r_norma_GSK);
		float p_param = (Ext_param3 * Semimajor_axis) / (2 * r_norma_GSK);


		float s_past_param = 0;
		float b_param = c_param + s_past_param;
		float s_param = 0;

		while (abs(s_param - s_past_param) < Accuracy)
		{
			s_past_param = s_param;
			b_param = c_param + s_past_param;
			s_param = asin((Focal_param * sin(2 * b_param)) / (sqrt(1 - Ext_param3 * pow(sin(b_param), 2))));

			//out_data_file << "test s_param" << s_param;
		}

		B_param = b_param;
		H_param = D_param * cos(B_param) + position_GSK[2] * sin(B_param) - Semimajor_axis * (sqrt(1 - Ext_param3 * pow(sin(B_param), 2)));
		
		return H_param;
	}


}


void FAPMech::calculate_density()
{
	density_night_param_first[0] = ND_upper_120 * (a0_first_level[0] + H_param * a1_first_level[0] + pow(H_param, 2) * a2_first_level[0] + pow(H_param, 3) * a3_first_level[0] 
		+ pow(H_param, 4) * a4_first_level[0] + pow(H_param, 5) * a5_first_level[0] + pow(H_param, 6) * a6_first_level[0]);
	density_night_param_first[1] = ND_upper_120 * (a0_first_level[1] + H_param * a1_first_level[1] + pow(H_param, 2) * a2_first_level[1] + pow(H_param, 3) * a3_first_level[1]
		+ pow(H_param, 4) * a4_first_level[1] + pow(H_param, 5) * a5_first_level[1] + pow(H_param, 6) * a6_first_level[1]);
	density_night_param_first[2] = ND_upper_120 * (a0_first_level[2] + H_param * a1_first_level[2] + pow(H_param, 2) * a2_first_level[2] + pow(H_param, 3) * a3_first_level[2]
		+ pow(H_param, 4) * a4_first_level[2] + pow(H_param, 5) * a5_first_level[2] + pow(H_param, 6) * a6_first_level[2]);
	density_night_param_first[3] = ND_upper_120 * (a0_first_level[3] + H_param * a1_first_level[3] + pow(H_param, 2) * a2_first_level[3] + pow(H_param, 3) * a3_first_level[3]
		+ pow(H_param, 4) * a4_first_level[3] + pow(H_param, 5) * a5_first_level[3] + pow(H_param, 6) * a6_first_level[3]);
	density_night_param_first[4] = ND_upper_120 * (a0_first_level[4] + H_param * a1_first_level[4] + pow(H_param, 2) * a2_first_level[4] + pow(H_param, 3) * a3_first_level[0]
		+ pow(H_param, 4) * a4_first_level[4] + pow(H_param, 5) * a5_first_level[4] + pow(H_param, 6) * a6_first_level[4]);
	density_night_param_first[5] = ND_upper_120 * (a0_first_level[5] + H_param * a1_first_level[5] + pow(H_param, 2) * a2_first_level[5] + pow(H_param, 3) * a3_first_level[5]
		+ pow(H_param, 4) * a4_first_level[5] + pow(H_param, 5) * a5_first_level[5] + pow(H_param, 6) * a6_first_level[5]);
	density_night_param_first[6] = ND_upper_120 * (a0_first_level[6] + H_param * a1_first_level[6] + pow(H_param, 2) * a2_first_level[6] + pow(H_param, 3) * a3_first_level[6]
		+ pow(H_param, 4) * a4_first_level[6] + pow(H_param, 5) * a5_first_level[6] + pow(H_param, 6) * a6_first_level[6]);


	density_night_param_second[0] = ND_upper_120 * (a0_second_level[0] + H_param * a1_second_level[0] + pow(H_param, 2) * a2_second_level[0] + pow(H_param, 3) * a3_second_level[0]
		+ pow(H_param, 4) * a4_second_level[0] + pow(H_param, 5) * a5_second_level[0] + pow(H_param, 6) * a6_second_level[0]);
	density_night_param_second[1] = ND_upper_120 * (a0_second_level[1] + H_param * a1_second_level[1] + pow(H_param, 2) * a2_second_level[1] + pow(H_param, 3) * a3_second_level[1]
		+ pow(H_param, 4) * a4_second_level[1] + pow(H_param, 5) * a5_second_level[1] + pow(H_param, 6) * a6_second_level[1]);
	density_night_param_second[2] = ND_upper_120 * (a0_second_level[2] + H_param * a1_second_level[2] + pow(H_param, 2) * a2_second_level[2] + pow(H_param, 3) * a3_second_level[2]
		+ pow(H_param, 4) * a4_second_level[2] + pow(H_param, 5) * a5_second_level[2] + pow(H_param, 6) * a6_second_level[2]);
	density_night_param_second[3] = ND_upper_120 * (a0_second_level[3] + H_param * a1_second_level[3] + pow(H_param, 2) * a2_second_level[3] + pow(H_param, 3) * a3_first_level[3]
		+ pow(H_param, 4) * a4_second_level[3] + pow(H_param, 5) * a5_second_level[3] + pow(H_param, 6) * a6_second_level[3]);
	density_night_param_second[4] = ND_upper_120 * (a0_second_level[4] + H_param * a1_second_level[4] + pow(H_param, 2) * a2_second_level[4] + pow(H_param, 3) * a3_second_level[0]
		+ pow(H_param, 4) * a4_second_level[4] + pow(H_param, 5) * a5_second_level[4] + pow(H_param, 6) * a6_second_level[4]);
	density_night_param_second[5] = ND_upper_120 * (a0_second_level[5] + H_param * a1_second_level[5] + pow(H_param, 2) * a2_second_level[5] + pow(H_param, 3) * a3_second_level[5]
		+ pow(H_param, 4) * a4_second_level[5] + pow(H_param, 5) * a5_second_level[5] + pow(H_param, 6) * a6_second_level[5]);
	density_night_param_second[6] = ND_upper_120 * (a0_second_level[6] + H_param * a1_second_level[6] + pow(H_param, 2) * a2_second_level[6] + pow(H_param, 3) * a3_second_level[6]
		+ pow(H_param, 4) * a4_second_level[6] + pow(H_param, 5) * a5_second_level[6] + pow(H_param, 6) * a6_second_level[6]);

}

void FAPMech::run_simulation()
{

	std::fstream out_data_file;
	out_data_file.open("C:\\Users\\1\\Desktop\\MKP_lab_1\\space_mission_params_file.txt", std::ios::trunc | std::ios::out | std::ios::in);

	calculate_anomalies();
	calculate_AGSK_pos();
	calculate_velocities();
	calculate_GSK_pos();
	calculate_LBH_params();



	out_data_file << "pos [x: " << position_AGSK[0] << "y: " << position_AGSK[1] << "z: " << position_AGSK[2] << "]" << "\n";
	out_data_file << "pos [x: " << position_GSK[0] << "y: " << position_GSK[1] << "z: " << position_GSK[2] << "]" << "\n";
	out_data_file << "pos [x: " << position_LBH[0] << "y: " << position_LBH[1] << "z: " << position_LBH[2] << "]" << "\n";
	out_data_file << "transversial_vel: " << transversial_velocity << "\n";
	out_data_file << "rad_vel: " << radial_velocity << "\n";
	out_data_file << "r_normal_AGSK: " << r_normal_AGSK << "\n";
	out_data_file << "theta: " << theta_param << "\n";
	out_data_file << "ext_param" << ext_param << "\n";
	out_data_file << "u_param" << u_param << "\n";
	out_data_file << "time: " << time << "\n";
	out_data_file << "N_param: " << N_param << "\n";
	out_data_file << "B_param: " << B_param << "\n";
	out_data_file << "H_param: " << H_param << "\n";
	out_data_file << "La_param: " << La_param << "\n";
	out_data_file << "L_param: " << L_param << "\n";
	out_data_file << "D_param: " << D_param << "\n";
	//out_data_file << "density: " << density_param << "\n";
	//out_data_file << "night density param: " << density_night_param << "\n";

	out_data_file.close();




}