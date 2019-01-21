#include <iostream>
#include <math.h>
#include <cmath>
#include <vector>
#include <array>
#include <fstream>

#define check(x) if (isnan(x)  || isinf(x) ) std::cout << x << " diverges" << '\n'

constexpr auto PI = 3.14159265;

std::vector<double> Current_Iteration;
std::vector<double> Previous_Iteration;
double Difference[15];


double lambda;
double lambda_TR, Adv_Ratio;
double Beta_0, Beta_1s, Beta_1c, Beta_p;
double Theta_0, Theta_1s, Theta_1c, Theta_t, Theta_FP, Theta_0_TR;
double alpha_s, phi_s;
double C_T, C_W, C_P, C_H, C_Y;;
double CG_StationLine, CG_WaterLine;
double M_XF, M_YF, YF;		// Moments and Y force
double Weight;
double Drag, Lift;
double rho;
double Cd_0;
double flatPlateArea, Area, Area_TR;
double Velocity;
double v_0, v;

double MR_Radius, MR_Rotational_Speed, MR_Lock, MR_Blade_Twist;
double MR_Precone, MR_Solidity, MR_LC_Slope, MR_Drag_Coeff;
double MR_Hub_Stationline, MR_Hub_Waterline;

double TR_Radius, TR_Rotational_Speed, TR_Solidity, TR_LC_Slope, TR_DP_Ratio;
double TR_Hub_Stationline, TR_Hub_Waterline, TR_Drag_Coeff, TR_Blade_Twist, C_T_TR;

double Fus_Pitch, Fus_Pitch_Alpha;
double Fus_Drag, Fus_Drag_Alpha;

double HS_Stationline, HS_Waterline;
double HS_Area, HS_Aspect_Ratio, HS_CLmax, HS_DP_Ratio, HS_LC_Slope, HS_Deltai;

double X_CG;	// x distance from cg to MR Hub
double Y_CG;	// y distance from cg to MR Hub
double h_CG;	// z distance from cg to MR Hub

double Drag_Fuselage, Drag_HS, Drag_MR, Drag_TR, Lift_HS; // TODO: Add other drags
double Dynamic_Pressure;

double M_Y_Fuselage;

double Adv_Ratio_TR;
double C_Q, Q_main, T_TR;
double lenght_h;	// distance between hs and cg
double alpha_horizontal;  // angle between lenght_h and horizontal
double M_Y_HS;
double del_h_cg_y;
double del_h_cg_x;
double M_Y_TR;
double max_lift_check;

void Get_HelicopterData() {

	Weight = 8000;				// lb
	CG_StationLine = 196 / 12.0;		// feet 
	CG_WaterLine = 73 / 12.0;		// feet

	MR_Radius = 22.0;				// feet
	MR_Rotational_Speed = 32.88;			// rad/s
	MR_Lock = 5.216;
	MR_Blade_Twist = -0.017453;		// rad
	MR_Precone = 0.048;			// rad
	MR_Solidity = 0.0651;
	MR_LC_Slope = 6.281;			// 1/rad
	MR_Drag_Coeff = 0.01;
	MR_Hub_Stationline = 200.0 / 12.0;		// feet
	MR_Hub_Waterline = 152.76 / 12.0;	// feet

	TR_Radius = 4.25;				// feet
	TR_Rotational_Speed = 168.44;			// rad/s
	TR_Blade_Twist = 0.0;				// rad
	TR_Solidity = 0.105;
	TR_LC_Slope = 6.28;				// 1/rad
	TR_Hub_Stationline = 520.7 / 12.0;		// feet
	TR_Hub_Waterline = 118.27 / 12.0;	// feet
	TR_Drag_Coeff = 0.01;
	TR_DP_Ratio = 1.0;

	HS_Stationline = 398.5 / 12.0;
	HS_Waterline = 56.0 / 12.0;
	HS_LC_Slope = 6.28;
	HS_Deltai = 0.03;
	HS_Area = 14.7;
	HS_DP_Ratio = 0.8;
	HS_Aspect_Ratio = 3.0;
	HS_CLmax = 1.2;

	X_CG = MR_Hub_Stationline - CG_StationLine;		// feet
	h_CG = MR_Hub_Waterline - CG_WaterLine;			// feet
	Y_CG = 0.0;  // !!!!! CHANGE 

	flatPlateArea = 20.5;		// sqft
	v = 1.0;
	v_0 = 1.0; // was 1 !!!!!!!!!!!!

	Fus_Pitch = -6.9;
	Fus_Pitch_Alpha = 280.405;
	Fus_Drag = 5.5;
	Fus_Drag_Alpha = -40.1;  // CHANGE

	Cd_0 = 0.01;
}

void Initial_Conditions() {

	// Set angles to zero
	Theta_0 = 0; Theta_1s = 0; Theta_1c = 0;
	Beta_0 = 0; Beta_1s = 0; Beta_1c = 0;

	Theta_0_TR = 0;

	alpha_s = 0; phi_s = 0;

	Adv_Ratio = 0;
	rho = 0.002377;

	// Approximate C_T = C_W
	C_W = Weight / (rho * PI * MR_Radius * MR_Radius  * MR_Rotational_Speed * MR_Rotational_Speed * MR_Radius * MR_Radius);
	C_T = C_W;
	C_H = 0;
	C_Y = 0;

	del_h_cg_y = CG_WaterLine - HS_Waterline;
	del_h_cg_x = HS_Stationline - CG_StationLine;
	lenght_h = sqrt(del_h_cg_x * del_h_cg_x + del_h_cg_y * del_h_cg_y);
	alpha_horizontal = atan(del_h_cg_y / del_h_cg_x);

	// Initialize Tail Rotor Thrust Coefficient
	C_T_TR = 0;

}

void IterateLambda() {
	double error = 584;
	double allowable_error = 0.0000001;
	double lambda1 = 1;
	int max_iterations = 10000000;
	int num_of_iterations = 0;

	while (error > allowable_error)
	{
		lambda = Adv_Ratio * tan(alpha_s) + C_T / (2.0 * sqrt(Adv_Ratio * Adv_Ratio + lambda1 * lambda1));
		error = abs(lambda - lambda1);
		lambda1 -= abs(lambda - lambda1) / 4.0;
		num_of_iterations++;
		if (num_of_iterations > max_iterations) {
			std::cout << "Lambda did not converge in " << max_iterations << " iterations";
			std::cout << "lambda: " << lambda << '\n' << "lamba0: " << lambda1 << std::endl;
			system("Pause");
		}
	}
	Current_Iteration.push_back(lambda);
}

void Find_Theta0() {
	Theta_0 = (2.0 * C_T / (MR_Lock * MR_LC_Slope) - MR_Blade_Twist / 4.0 * (1.0 + Adv_Ratio * Adv_Ratio) + Adv_Ratio / 2.0 * Theta_1s + lambda / 2.0) * 3.0 / (1.0 + 3.0 / 2.0 * Adv_Ratio * Adv_Ratio);
	check(Theta_0);
	Current_Iteration.push_back(Theta_0);
}

void Find_Beta0() {
	Beta_0 = (v_0 * v_0) / (v * v) * MR_Precone + MR_Lock / (v * v) * (Theta_0 / 8.0 * (1.0 + Adv_Ratio * Adv_Ratio) + MR_Blade_Twist / 10.0 * (1.0 + 5.0 / 6.0 * Adv_Ratio * Adv_Ratio) - Adv_Ratio / 6.0 * Theta_1s - lambda / 6.0);
	check(Beta_0);
	Current_Iteration.push_back(Beta_0);
}

void Find_Theta1c() {
	Theta_1c = Beta_1s - (4.0 / 3.0 * Adv_Ratio * Beta_0 - (v * v - 1.0) * 8.0 / MR_Lock * Beta_1c) / (1.0 + Adv_Ratio * Adv_Ratio / 2.0);																		 //4
	check(Theta_1c);
	Current_Iteration.push_back(Theta_1c);
}

void Find_Theta1s() {
	Theta_1s = (8.0 / 3.0 * Adv_Ratio * (Theta_0 + 3.0 / 4.0 * MR_Blade_Twist - 3.0 / 4.0 * lambda) + (v * v - 1) * 8.0 / MR_Lock * Beta_1s) / (1.0 + 3.0 / 2.0 * Adv_Ratio * Adv_Ratio) - Beta_1c;
	check(Theta_1s);
	Current_Iteration.push_back(Theta_1s);
}

void Calculate_Drag() {
	Velocity = Adv_Ratio * MR_Rotational_Speed * MR_Radius;
	Dynamic_Pressure = 0.5 * rho * Velocity * Velocity;
	Drag_Fuselage = (Fus_Drag + Fus_Drag_Alpha * alpha_s * -1.0) * Dynamic_Pressure;
	check(Drag_Fuselage);
	Drag_MR = Cd_0 * Dynamic_Pressure * MR_Radius * MR_Radius * PI;
	Drag_TR = Cd_0 * Dynamic_Pressure * TR_DP_Ratio * TR_Radius * TR_Radius * PI;
	Drag_HS = HS_DP_Ratio * Dynamic_Pressure * HS_Area * (HS_LC_Slope *  HS_LC_Slope * alpha_s * alpha_s / PI / HS_Aspect_Ratio * (1.0 + HS_Deltai) + 0.0065);
	Drag = Drag_Fuselage + Drag_MR + Drag_TR + Drag_HS;
	max_lift_check = HS_LC_Slope * alpha_s;
	if (abs(max_lift_check) >= HS_CLmax) max_lift_check = HS_CLmax;

	Lift_HS = HS_DP_Ratio * Dynamic_Pressure * HS_Area * max_lift_check;
}

void Calculate_Forces_and_Moments() {
	M_Y_Fuselage = (Fus_Pitch - Fus_Pitch_Alpha * alpha_s) * Dynamic_Pressure;
	M_Y_HS = Lift_HS * lenght_h * cos(alpha_s) + Drag_HS * lenght_h * sin(alpha_s);
	M_YF = M_Y_Fuselage + M_Y_HS;
	YF = C_T_TR * rho * PI * TR_Radius * TR_Radius * TR_Rotational_Speed * TR_Rotational_Speed * TR_Radius * TR_Radius;
	M_XF = -1.0 * (YF * TR_Hub_Waterline - CG_WaterLine);
}

void Find_Beta1c() {
	Beta_1c = (X_CG / h_CG - M_YF / Weight / h_CG - C_H / C_T) / (1.0 + (v * v - 1.0) / MR_Lock / (2.0 * h_CG / MR_Radius * C_T / MR_Solidity / MR_LC_Slope));
	check(Beta_1c);
	Current_Iteration.push_back(Beta_1c);
}

void Find_Beta1s() {
	Beta_1s = (M_XF / Weight / h_CG - Y_CG / h_CG - C_Y / C_T) / (1.0 + (v * v - 1.0) / MR_Lock / (2.0 * h_CG / MR_Radius * C_T / MR_Solidity / MR_LC_Slope));
	check(Beta_1s);
	Current_Iteration.push_back(Beta_1s);
}

void Find_Alpha_s() {
	alpha_s = (X_CG / h_CG - M_YF / Weight / h_CG + (v * v - 1.0) / MR_Lock / (2.0 * h_CG / MR_Radius * C_T / MR_Solidity / MR_LC_Slope) * C_H / C_T) / (1.0 - (v * v - 1.0) / MR_Lock / (2.0 * h_CG / MR_Radius * C_T / MR_Lock / MR_LC_Slope)) + Drag / Weight;
	check(alpha_s);
	if (phi_s > 0.7) alpha_s = 0.1;
	if (isnan(alpha_s)) alpha_s = 0.1;
	Current_Iteration.push_back(alpha_s);
}

void Find_Phi_s() {
	phi_s = (Y_CG / h_CG - M_XF / Weight / h_CG - (v * v - 1.0) / MR_Lock / (2.0 * h_CG / MR_Radius * C_T / MR_Solidity / MR_LC_Slope)) / (1.0 - (v * v - 1.0) / MR_Lock / (2.0 * h_CG / MR_Radius * C_T / MR_Lock / MR_LC_Slope)) - YF / Weight;
	check(phi_s);
	if (phi_s > 2) phi_s = 0.1;
	if (isnan(alpha_s)) alpha_s = 0.1;
	Current_Iteration.push_back(phi_s);
}

// Assume C_H = C_H_tpp
void Find_C_H() {
	C_H = MR_Solidity * Cd_0 / 4.0 * Adv_Ratio + (MR_Solidity * MR_LC_Slope / 2.0) * (Adv_Ratio * lambda / 2.0 * (Theta_0 + 0.5 * MR_Blade_Twist) + 1.0 / 6.0 * Theta_1c * Beta_0 - 1.0 / 4.0 * Theta_1s * lambda + 1.0 / 4.0 * Adv_Ratio * Beta_0 * Beta_0);
	check(C_H);
	Current_Iteration.push_back(C_H);
}

//Assume C_Y = C_Y_tpp
void Find_C_Y() {
	C_Y = -1.0 * MR_Solidity * MR_LC_Slope / 2.0 * (3.0 / 4.0 * Adv_Ratio * Beta_0 * (Theta_0 + 2.0 / 3.0 * MR_Blade_Twist) - 1.0 / 4.0 * Theta_1c * lambda - 1.0 / 6.0 * Theta_1s * Beta_0 * (1.0 + 3.0 * Adv_Ratio * Adv_Ratio) - 3.0 / 2.0 * Adv_Ratio * Beta_0 * lambda);
	Current_Iteration.push_back(C_Y); // CHECK DIRECTION - checked
}

void Find_C_T() {
	Weight = 8000;
	Weight += Lift_HS;
	C_W = Weight / (rho * PI * MR_Radius * MR_Radius  * MR_Rotational_Speed * MR_Rotational_Speed * MR_Radius * MR_Radius);
	C_T = C_W / cos(alpha_s) - C_H * tan(alpha_s);
	if (abs(C_T) > 5.0 * C_W) C_T = C_W;
	Current_Iteration.push_back(C_T);
}

void Find_C_P() {		// check f/A										
	C_P = C_T * C_T / (2.0 * sqrt(Adv_Ratio * Adv_Ratio + lambda * lambda)) + MR_Solidity * Cd_0 / 8.0 * (1.0 + 4.6 * Adv_Ratio * Adv_Ratio) + 0.5 * Adv_Ratio * Adv_Ratio * Adv_Ratio * flatPlateArea / (PI * MR_Radius * MR_Radius);  // Check area
	Current_Iteration.push_back(C_P);
}

void Find_C_T_TR() {
	C_Q = C_P;
	Q_main = C_Q * rho * PI * MR_Radius * MR_Radius * MR_Rotational_Speed * MR_Rotational_Speed * MR_Radius * MR_Radius * MR_Radius;
	T_TR = Q_main / (TR_Hub_Stationline - CG_StationLine);
	C_T_TR = T_TR / (rho * PI * TR_Radius * TR_Radius * TR_Rotational_Speed * TR_Rotational_Speed * TR_Radius * TR_Radius * TR_Radius);
	Current_Iteration.push_back(C_T_TR);
}

void Iterate_TR() {
	double errortr = INT_MAX;
	double allowable_errortr = 0.00000001;
	double lambdatr2 = 1.0;
	int max_iterations = 1000000;
	int num_of_iterationstr = 0;
	lambda_TR = 5;
	while (errortr > allowable_errortr)	// CHECK
	{
		Adv_Ratio_TR = Adv_Ratio * (MR_Rotational_Speed * MR_Radius / TR_Rotational_Speed / TR_Radius);
		lambda_TR = C_T_TR / (2.0 * sqrt(Adv_Ratio_TR * Adv_Ratio_TR + lambdatr2 * lambdatr2));
		errortr = abs(lambda_TR - lambdatr2);
		lambdatr2 -= abs(lambda_TR - lambdatr2) / 4.0;
		num_of_iterationstr++;
		if (num_of_iterationstr > max_iterations) {
			std::cout << "LambdaTR did not converge in " << max_iterations << " iterations";
			std::cout << "lambda: " << lambda_TR << '\n' << "lamba0: " << lambdatr2 << std::endl;
			system("Pause");
		}
	}
	Current_Iteration.push_back(lambda_TR);
}

void Find_Theta_0_TR() { //CHECK - checked
	Theta_0_TR = (2.0 * C_T_TR / (TR_Solidity * TR_LC_Slope) + lambda_TR / 2.0) * 3.0 / (1.0 + 1.5 * Adv_Ratio_TR * Adv_Ratio_TR);
	Current_Iteration.push_back(Theta_0_TR);
}


int main() {

	std::ofstream t0;
	std::ofstream t1s;
	std::ofstream t1c;
	std::ofstream b0;
	std::ofstream b1s;
	std::ofstream b1c;
	std::ofstream t0tr;
	std::ofstream ct;
	std::ofstream cp;
	std::ofstream phi;
	std::ofstream alpha;

	t0.open("Theta0.txt");
	t1c.open("Theta1c.txt");
	t1s.open("Theta1s.txt");
	b0.open("Beta0.txt");
	b1s.open("Beta1s.txt");
	b1c.open("Beta1c.txt");
	t0tr.open("Theta0TR.txt");
	ct.open("CT.txt");
	cp.open("CP.txt");
	alpha.open("alpha.txt");
	phi.open("phi.txt");



	Get_HelicopterData();
	Initial_Conditions();


	for (Adv_Ratio = 0; Adv_Ratio < 0.41; Adv_Ratio += 0.01)
	{
		double allowable_difference{ 0.0000000001 };
		double residual{ 0 };
		bool SolutionConverged = false;
		bool FirstIteration = true;
		int number_of_iterations{ 0 };

		while (!SolutionConverged)
		{

			IterateLambda();
			Find_Theta0();
			Find_Beta0();
			Find_Theta1c();
			Find_Theta1s();
			Find_Beta1c();
			Find_Beta1s();
			Calculate_Drag();
			Calculate_Forces_and_Moments();
			Find_Alpha_s();
			Find_Phi_s();
			Find_C_H();
			Find_C_Y();
			Find_C_T();
			Find_C_P();
			Find_C_T_TR();
			Iterate_TR();
			Find_Theta_0_TR();

			if (FirstIteration) {
				Previous_Iteration = Current_Iteration;
				FirstIteration = false;
			}
			else {
				for (unsigned i = 0; i < Current_Iteration.size(); i++)
				{
					double difference = Current_Iteration[i] - Previous_Iteration[i];
					Difference[i] = abs(difference);
				}
				double max = (double)INT_MIN;
				for (double val : Difference) {
					if (max < val) max = val;
				}
				number_of_iterations++;
				if (max < allowable_difference) {
					residual = max;
					SolutionConverged = true;
				}
			}
			Previous_Iteration = Current_Iteration;
			lambda = 0;
			Current_Iteration.clear();
		}

		std::cout << "# of iterations " << number_of_iterations << std::endl << "Maximum error = " << residual << std::endl;

		std::cout << "Advance Ratio = " << Adv_Ratio << std::endl;
		std::cout << "Theta0 = " << Theta_0 * 180.0 / PI << "\t" << "Theta1c = " << Theta_1c * 180.0 / PI << "\t" << "Theta1s = " << Theta_1s * 180.0 / PI << "\t" << "Theta0_TR = " << Theta_0_TR * 180.0 / PI << std::endl;
		std::cout << "Beta0 = " << Beta_0 * 180.0 / PI << "\t" << "Beta1c = " << Beta_1c * 180.0 / PI << "\t" << "Beta1s = " << Beta_1s * 180.0 / PI << "\t" << "Alpha = " << alpha_s * 180.0 / PI << "\n\n\n";

		t0 << Theta_0 * 180.0 / PI << "\t" << Adv_Ratio << std::endl;
		t1s << Theta_1s * 180.0 / PI << "\t" << Adv_Ratio << std::endl;
		t1c << Theta_1c * 180.0 / PI << "\t" << Adv_Ratio << std::endl;
		t0tr << Theta_0_TR * 180.0 / PI << "\t" << Adv_Ratio << std::endl;
		b0 << Beta_0 * 180.0 / PI << "\t" << Adv_Ratio << std::endl;
		b1s << Beta_1s * 180.0 / PI << "\t" << Adv_Ratio << std::endl;
		b1c << Beta_1c * 180.0 / PI << "\t" << Adv_Ratio << std::endl;
		alpha << alpha_s * -180.0 / PI << "\t" << Adv_Ratio << std::endl;
		phi << phi_s * 180.0 / PI << "\t" << Adv_Ratio << std::endl;
		cp << C_P << "\t" << Adv_Ratio << std::endl;
		ct << C_T << "\t" << Adv_Ratio << std::endl;
	}

	system("Pause");
}


