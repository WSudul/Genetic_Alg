

#include <array>
#include <iostream>
#include <string.h>
#include <string>
#include <fstream>
#include <vector>
#include <sstream>
#include <cmath>

class ReservoirModel
{
public:
	//explicit ReservoirModel(QWidget *parent = 0);
	ReservoirModel();
	
	void setData(std::vector<double> &p, std::vector<double> &g_p, std::vector<double> &t, std::vector<double> &w_p);
	

	~ReservoirModel();
	int fy_gas();
	void non_vol(std::array<int, 4>& intParams, std::array<float, 6> &floatParams);

	double zpar(double p_ri);
	double z_func(double t_pr, double p_pr, double y);
	double z_func_der(double t_pr, double p_pr, double y);
	double doCalculcations(std::array<int, 4> &ar1, std::array<float, 6> &ar2);
private:
	double a_lin(double x, double y, double x2, double xy, int n);
	double b_lin(double x, double y, double x2, double xy, int n);
	double B_w(double p, double t);
	double t_pc();
	double p_pc();
	double pitzer_pc();
	void SRK_parameters(double t_pr, double p_pr);
	void RK_parameters(double t_pr, double p_pr);
	double E_WaA();
	//void on_pushButton_clicked();
	double m, acc_factor, alfa, A, B;

	//void on_pushButton_2_clicked();

	
	std::vector<double> p, g_p, t, w_p; //vector with data!


	void temperature();

	//void on_pushButton_4_clicked();

	//void on_checkBox_toggled(bool checked);



	//void on_pushButton_5_clicked();

	//void on_q_DoubleSpinBox_valueChanged(double arg1);

	//void on_time_SpinBox_valueChanged(int arg1);

	//void on_time_DoubleSpinBox_valueChanged(double arg1);

	double U(double r_o, double h, double f, double por, double c);
	double dim_t(double k, double t, double por, double lep, double c_t, double r);
	double WeD(double t, double r_ed);
	double Bg(double z, double p, double t);
	double F_HO(double Gp, double z, double p, double W_p);
	double Efw(double pr, double c_w, double s_w, double c_f, double zi, double pi, double ti);
	//Ui::ReservoirModel *ui;
	void funkcja1();
	//void funkcja2();
	
	double funkcja3(std::array<int, 4>& intParams, std::array<float, 6> &floatParams);

	//"global" values

	double h;// = 20;
	double por;// = 0.2;
	double k;// = 15 * 0.986923*pow(0.1, 15);
	double mi_w;// = 1 * 0.001;
	double c_w;// = 0.5*pow(0.1, 9);
	double c_f;// = 0.6*pow(0.1, 9);
	double f;// = 0.75;
	double r_res;// = 300;
	double r_aquifer;// = 700;
	double t_r;// = 273.15;
	double r_ed;// = r_aquifer / r_res;
	double s_w;// = 0.23;
	double pi = 3.1415926535897;


	const long year = 31536000; //czas trwania roku w sekundach
	double p_c[15] = {
		4600200.00,
		4883900.00,
		4250344.828,
		3648965.517,
		3800000.0,
		3390000.0,
		3372413.793,
		3013793.103,
		2740000.0,
		1290000.0,
		3390000.0,
		5040000.0,
		7376500.00,
		8940000.00,
		22062068.97,
	};

	double t_c[15] = {
		190.6,
		305.4,
		369.5000,
		407.8333,
		425.0000,
		460.0000,
		469.1111,
		506.8889,
		540.0000,
		33.00000,
		126.0000,
		154.0000,
		304.2,
		373.2,
		646.8889,
	};

	double pitzer[15] = {
		0.011,    //    Methan, CH4
		0.098,    //    Ethan, C2H6
		0.149,    //    Propan, C3H8
		0.177,    //    i-Butane, C4H10
		0.197,    //    n-Butane, C4H10
		0.226,    //    i-Pentane C5H12
		0.251,    //    n-Pentane C5H12
		0.304,    //    Hexane C6H14
		0.346,    //    Heptane C7H16
		-0.215,    //    Hydogen, H2
		0.037,    //    Nitrogen, N2
		0.020,    //    Oxygen, O2
		0.224,    //    Carbon dioxid, CO2
		0.096,    //    Hydrogensulfid, H2S
		0.344,    //    Dihydrogenoksid, H2O
	};

	double y_gas[15] = {
		0.7,
		0.10,
		0.1,
		0.05,
		0.05,
		0,
		0.0,
		0,
		0,
		0,
		0,
		0,
		0,
		0,
		0,

	};

};



