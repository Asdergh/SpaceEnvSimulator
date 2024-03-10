#include "FAPMech.h"


void FAPMech::calculate_anomalies()
{
	int simulation = 0;
	while (abs(ext_param - ext_past_param) <= Accuracy)
	{
		if (simulation == 0)
		{
			ext_past_param = M_param;
		}

		else
		{
			ext_param = M_param + Ext_param * sin(ext_past_param);
			ext_past_param = ext_param;
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

void FAPMech::calculate_velocityes()
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
		float s_param;

		while (abs(s_param - s_past_param) < Accuracy)
		{
			s_param = asin((Focal_param * sin(2 * b_param)) / (sqrt(1 - Ext_param3 * pow(sin(b_param), 2))));
			s_past_param = s_param;
			b_param = c_param + s_past_param;

		}

		B_param = b_param;
		H_param = D_param * cos(B_param) + position_GSK[2] * sin(B_param) - Semimajor_axis * (sqrt(1 - Ext_param3 * pow(sin(B_param), 2)));
		
		return H_param;
	}


}