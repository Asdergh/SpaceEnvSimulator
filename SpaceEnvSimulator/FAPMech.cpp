#include "FAPMech.h"



FAPMech::FAPMech()
{

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
FAPMech::~FAPMech()
{
	delete[] position_GSK;
	delete[] position_AGSK;
	delete[] velocity_vector;
	delete[] density_night_param;

	delete[] a0_first_level;
	delete[] a1_first_level;
	delete[] a2_first_level;
	delete[] a3_first_level;
	delete[] a4_first_level;
	delete[] a5_first_level;
	delete[] a6_first_level;

	delete[] a0_second_level;
	delete[] a1_second_level;
	delete[] a2_second_level;
	delete[] a3_second_level;
	delete[] a4_second_level;
	delete[] a5_second_level;
	delete[] a6_second_level;

	for (int item_index{ 0 }; item_index < 6; item_index++)
	{
		delete[] indignant_acceleration_tensor[item_index];
	}

	delete[] indignant_acceleration_tensor;

	std::cout << "\n" << "data was delted!!!" << "\n";
}
void FAPMech::calculate_anomalies()
{
	ext_past_param = M_param;
	ext_param = M_param + Ext_param * sin(ext_past_param);

	while (fabs(ext_param - ext_past_param) > Accuracy)
	{
		ext_past_param = ext_param;
		ext_param = M_param + Ext_param * sin(ext_past_param);

		//out_data_file << "test_ext_param: " << ext_param << std::endl;

	}

	
	theta_param = (2.0 * atan(sqrt((1.0 + Ext_param) / (1.0 - Ext_param)) * tan(ext_param / 2.0))) * (PI / 180.0);
	u_param = theta_param + omega_param;
	r_normal_AGSK = Focal_param * (1.0 - pow(Ext_param, 2.0)) / (1.0 + Ext_param * cos(theta_param));
}

void FAPMech::calculate_AGSK_pos()
{
	position_AGSK[0] = r_normal_AGSK * (cos(u_param) * cos(Omega_param) - sin(u_param) * sin(Omega_param) * cos(I_param));
	position_AGSK[1] = r_normal_AGSK * (cos(u_param) * sin(Omega_param) + sin(u_param) * cos(Omega_param) * cos(I_param));
	position_AGSK[2] = r_normal_AGSK * (sin(u_param) * sin(I_param));
}

void FAPMech::calculate_velocities()
{
	transversial_velocity = sqrt(Earth_grav_param / Focal_param) * Ext_param * sin(theta_param);
	radial_velocity = sqrt(Earth_grav_param / Focal_param) * (1.0 + Ext_param * cos(theta_param));
	velocity_normal = sqrt(Earth_grav_param / Focal_param) * sqrt(1.0 + 2.0 * Ext_param * cos(theta_param) + pow(Ext_param, 2.0));
	velocity_vector[0] = -transversial_velocity;
	velocity_vector[1] = -radial_velocity;
	velocity_vector[2] = 0.0;
}

void FAPMech::calculate_GSK_pos()
{
	float ven_angle = Earth_anguler_vel * time;

	position_GSK[0] = (position_AGSK[0] * cos(ven_angle) + position_AGSK[1] * sin(ven_angle));
	position_GSK[1] = (-position_AGSK[0] * sin(ven_angle) + position_AGSK[1] * cos(ven_angle));
	position_GSK[2] = position_AGSK[2];

}

float FAPMech::calculate_LBH_params()
{
	
	D_param = sqrt(pow(position_GSK[0], 2.0) + pow(position_GSK[1], 2.0));
	N_param = SM_Earth / sqrt(1.0 - Ext_param3 * pow(sin(B_param), 2.0));

	if (D_param == 0)
	{
		B_param = (PI / 2.0) * position_GSK[2] / fabs(position_GSK[2]);
		L_param = 0.0;
		H_param = position_GSK[2] * sin(B_param) - SM_Earth * sqrt(1.0 - Ext_param3 * pow(sin(B_param), 2.0));

		position_LBH[0] = L_param;
		position_LBH[1] = B_param;
		position_LBH[2] = H_param;

		return H_param;
	}

	else if (D_param > 0.0)
	{
		La_param = asin(position_GSK[1] / D_param);

		if (position_GSK[1] < 0 and position_GSK[0] > 0)
		{
			L_param = 2.0 * PI - La_param;
		}

		else if (position_GSK[1] < 0 and position_GSK[0] < 0)
		{
			L_param = PI + La_param;
		}

		else if (position_GSK[1] > 0 and position_GSK[0] < 0)
		{
			L_param = PI - La_param;
		}

		else if (position_GSK[1] > 0 and position_GSK[0] > 0)
		{
			L_param = La_param;
		}

		if (position_GSK[2] == 0) {

			B_param = 0;
			H_param = D_param - SM_Earth;

			position_LBH[0] = L_param;
			position_LBH[1] = B_param;
			position_LBH[2] = H_param;

			return H_param;
		}

		else
		{
			float r_norma_GSK = sqrt(pow(position_GSK[0], 2.0) + pow(position_GSK[1], 2.0) + pow(position_GSK[2], 2.0));
			float c_param = asin(position_GSK[2] / r_norma_GSK);
			float p_param = (Ext_param3 * SM_Earth) / (2 * r_norma_GSK);


			float s_past_param = 0;
			float b_param = c_param + s_past_param;
			float s_param = 0;

			while (abs(s_param - s_past_param) > Accuracy)
			{
				s_past_param = s_param;
				b_param = c_param + s_past_param;
				s_param = asin((Focal_param * sin((2.0 * b_param) * (PI / 180.0))) / (sqrt(1 - Ext_param3 * pow(sin(b_param * (PI / 180.0)), 2))));

				//out_data_file << "test s_param" << s_param;
			}

			B_param = b_param;
			H_param = D_param * cos(B_param) + position_GSK[2] * sin(B_param) - SM_Earth * (sqrt(1 - Ext_param3 * pow(sin(B_param), 2.0)));

			position_LBH[0] = L_param;
			position_LBH[1] = B_param;
			position_LBH[2] = H_param;

			return H_param;
		}
	}

	

	


}


void FAPMech::calculate_density()
{

	if (H_param < 500)
	{
		for (int den_index{ 0 }; den_index < 7; den_index++)
		{
			density_night_param[den_index] = ND_upper_120 * exp(a0_first_level[den_index] + a1_first_level[den_index] * H_param + a2_first_level[den_index] * pow(H_param, 2)
				+ a3_first_level[den_index] * pow(H_param, 3) + a4_first_level[den_index] * pow(H_param, 4) + a5_first_level[den_index] * pow(H_param, 5) + a6_first_level[den_index] * pow(H_param, 6));
		}

	}
	
	else if (H_param > 500)
	{
		for (int den_index{ 0 }; den_index < 7; den_index++)
		{
			density_night_param[den_index] = ND_upper_120 * exp(a0_second_level[den_index] + a1_second_level[den_index] * H_param + a2_second_level[den_index] * pow(H_param, 2)
				+ a3_second_level[den_index] * pow(H_param, 3) + a4_second_level[den_index] * pow(H_param, 4) + a5_second_level[den_index] * pow(H_param, 5) + a6_second_level[den_index] * pow(H_param, 6));
		}
	}

	

}

void FAPMech::calculate_acceleration()
{
	for (int acc_index{ 0 }; acc_index < 7; acc_index++)
	{
		double simple_acceleration_T = -1.0 * (Sigma_coeff * density_night_param[acc_index] * transversial_velocity);
		double simple_acceleration_S = -1.0 * (Sigma_coeff * density_night_param[acc_index] * radial_velocity);
		double abs_acceleration = sqrt(pow(simple_acceleration_T, 2.0) + pow(simple_acceleration_S, 2.0));
		double* simple_acceleration_vector = new double[3] {simple_acceleration_T, simple_acceleration_S, abs_acceleration};

		indignant_acceleration_tensor[acc_index] = simple_acceleration_vector;
		
		std::cout << indignant_acceleration_tensor[acc_index];
	}
}

void FAPMech::run_simulation()
{

	std::fstream out_data_file;
	out_data_file.open("C:\\Users\\1\\Desktop\\MKP_lab_1\\space_mission_params_file.txt", std::ios::trunc | std::ios::out | std::ios::in);

	calculate_anomalies();
	calculate_AGSK_pos();
	calculate_velocities();
	





	time = 0;
	calculate_GSK_pos();
	H_param = calculate_LBH_params();
	calculate_density();
	calculate_acceleration();





	time = max_time / 2;
	calculate_GSK_pos();
	H_param = calculate_LBH_params();
	calculate_density();
	calculate_acceleration();
	



	out_data_file << "pos AGSK [x: " << position_AGSK[0] << "y: " << position_AGSK[1] << "z: " << position_AGSK[2] << "]" << "\n";
	out_data_file << "pos GSK[x: " << position_GSK[0] << "y: " << position_GSK[1] << "z: " << position_GSK[2] << "]" << "\n";
	out_data_file << "pos LBH[x: " << position_LBH[0] << "y: " << position_LBH[1] << "z: " << position_LBH[2] << "]" << "\n";
	out_data_file << "\n---------------------------------" << "velocity" << "----------------------------------------------\n";
	out_data_file << "transversial_vel: " << transversial_velocity << "\n";
	out_data_file << "rad_vel: " << radial_velocity << "\n";
	out_data_file << "abs vel: " << velocity_normal << "\n";
	out_data_file << "\n---------------------------------" << "trajactory params" << "----------------------------------------------\n";
	out_data_file << "r_normal_AGSK: " << r_normal_AGSK << "\n";
	out_data_file << "theta: " << theta_param << "\n";
	out_data_file << "ext_param: " << ext_param << "\n";
	out_data_file << "u_param: " << u_param << "\n";
	out_data_file << "time: " << time << "\n";
	out_data_file << "\n---------------------------------" << "LBH params" << "----------------------------------------------\n";
	out_data_file << "N_param: " << N_param << "\n";
	out_data_file << "B_param: " << B_param << "\n";
	out_data_file << "H_param: " << H_param << "\n";
	out_data_file << "La_param: " << La_param << "\n";
	out_data_file << "L_param: " << L_param << "\n";
	out_data_file << "D_param: " << D_param << "\n";
	out_data_file << "\n---------------------------------" << "density" << "----------------------------------------------\n";
	out_data_file << "120 density: ";
	for (int den_index{ 0 }; den_index < 6; den_index++) { out_data_file << " " << density_night_param[den_index]; }
	out_data_file << "\n---------------------------------" << "acceleration" << "----------------------------------------------\n";
	out_data_file << "acceleration:\n";
	for (int acc_i{ 0 }; acc_i < 7; acc_i++)
	{
		out_data_file << "\tacc with density number: " << acc_i << "\twith cores: [";
		for (int acc_j{ 0 }; acc_j < 3; acc_j++)
		{
			out_data_file << " " << indignant_acceleration_tensor[acc_i][acc_j];
		}
		out_data_file << "\t]\n";
	}


	out_data_file.close();




}