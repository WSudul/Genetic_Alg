#include "stdafx.h"
#include "ReservoirModel.h"









double ReservoirModel::a_lin(double x,double y,double x2, double xy, int n)
{
    return (n*xy-x*y)/(n*x2-x*x);
}

double ReservoirModel::b_lin(double x,double y,double x2, double xy, int n)
{
    double n_new=n;
    return (1/n_new*(y-((n_new*xy-x*y)/(n_new*x2-x*x))*x));
}




double ReservoirModel::B_w(double p,double t) //McCain, W.D. Jr.: McCain, W.D. Jr. 1990. The Properties of Petroleum Fluids, second edition. Tulsa, Oklahoma: PennWell Books.
{
    float t_fah=(t-32)/1.8;
    float p_psi=14.5037738*p;
    double deltaV_wp=-0.0010001+1.33391*pow(10,-4)*t_fah+5.50654*pow(10,-7)*pow(t_fah,2);
    double deltaV_wT=-1.95301*pow(10,-9)*p_psi*t_fah-1.72834*pow(10,-13)*pow(p_psi,2)*t_fah-3.58922*pow(10,-7)*p_psi
            -2.25341*pow(10,-10)*pow(p_psi,2);

    return (1+deltaV_wp)*(1+deltaV_wT);

}




double ReservoirModel::t_pc()
{
    double sum=0;
    for(int i=0;i<15;i++)
    {
        sum=sum+t_c[i]*y_gas[i];
    }
    return sum;
}

double ReservoirModel::p_pc()
{
    double sum=0;
    for(int i=0;i<15;i++)
    {
        sum=sum+p_c[i]/(pow(10,5))*y_gas[i];
    }
    return sum;

}

double ReservoirModel::pitzer_pc()
{
    double sum=0;
    for(int i=0;i<15;i++)
    {
        sum=sum+pitzer[i]*y_gas[i];
    }
    return sum;

}



void ReservoirModel::SRK_parameters(double t_pr,double p_pr)
{
    acc_factor=pitzer_pc();
    m=0.48508+1.55171*acc_factor-0.15613*pow(acc_factor,2);
    alfa=pow((1+m*(1-pow(t_pr,0.5))),2);
    A=alfa*0.42747*p_pr/(pow(t_pr,2));
    B=0.08664*p_pr/t_pr;

}
void ReservoirModel::RK_parameters(double t_pr, double p_pr)
{
    A=0.42748*p_pr/(pow(t_pr,2.5));
    B=0.08664*p_pr/t_pr;
    //std::cout<<"RK A="<<A<<"\tRK B=|"<<B;

}

double ReservoirModel::E_WaA() //H2S, CO2 correction factor
{
    double A=y_gas[12]+y_gas[13];
    double B=y_gas[13];
    double E_F=120*(pow(A,0.9)-pow(A,1.6))+15*(pow(B,0.5)-pow(B,4.0));
    /*std::cout<<"E_WA RESULTS!";
    std::cout<<"\t\tA="<<A<<" B="<<B;
    std::cout<<"\t\tE_F="<<E_F;*/
    return E_F*(5.0/9.0);
}


double ReservoirModel::U(double r_o,double h,double f,double por,double c)
{

    //std::cout<<"U="<<(2*pi*f*por*h*c*pow(r_o,2));
    return (2*pi*f*por*h*c*pow(r_o,2));

}

double ReservoirModel::dim_t(double k,double t,double por,double lep,double c_t,double r)  //c_t - dotyczy aquifer = c_w+c_f
{
    return (k*t)/(por*lep*c_t*pow(r,2));
}



double ReservoirModel::WeD(double t,double r_ed)
{
    double const_t=k/(por*mi_w*(c_w+c_f)*pow(r_res,2));
    double t_d=const_t*t*year;
    double t_dkr=0.4*pow((r_ed-1),2);
    double result;
    //std::cout<<"t_d="<<t_d<<" dla t="<<t;
    if(t_d>t_dkr*100) //TRUE=ogranczony
    {
        double J=pow(r_ed,4)*log(r_ed)/(pow(r_ed,2)-1.0)+0.25*(1-3*pow(r_ed,2));
        result=0.5*pow((r_ed-1),2)*(1-exp(-2*t_d/J));
    }
    else //NIEOGRANICZONY
    {
            if(t_d<0.01)
        {
            result=2*pow(t_d/pi,0.5);
        }
        else if(t_d<=200)
        {
            result=(1.2838*sqrt(t_d)+1.19327*t_d+0.269872*pow(t_d,1.5)
                    +0.00855294*pow(t_d,2))/(1+0.616599*sqrt(t_d)+0.0413008*t_d);
        }
        else
        {
            result=(-4.29881+2.02566*t_d)/(log(t_d));
        }
    }
    //std::cout<<"WeD(t_d)="<<result;
    return result;


}

double ReservoirModel::Bg(double z,double p,double t)
{
    double B_g=(z*t*1.01325/(273.15*p));
    //std::cout<<"B_g("<<p<<")="<<B_g;
    return B_g;

}

double ReservoirModel::F_HO(double Gp,double z,double p,double W_p)
{

    return Gp*Bg(z,p,t_r)+W_p*B_w(p,t_r);

}

double ReservoirModel::Efw(double pr,double c_w,double s_w,double c_f,double zi,double pi,double ti)
{
    if(pi==pr)
    {
        return 0;
    }
    return Bg(zi,pi,ti)*(c_w*s_w+c_f)/(1-s_w)*(pi-pr)*pow(10,5);
};




/*!
#TODO fix input
*/
void ReservoirModel::funkcja1() //otwieranie plikow i sprawdzenie czy mamy taka sama ilosc danych
{
    int n=5;  //minimalna liczba danych

    int sTotal1=0;
    int sTotal2=0;
    int sTotal3=0;
    int x=0;
    /*std::std::fstream file("p_data.txt");
    std::fstream file2("g_p_data.txt");
    std::fstream file3("time_data.txt");

    file.open(QIODevice::ReadOnly);
    std::string line[1000];
    QTextStream in(&file);
    while(!in.atEnd())
    {
        line[sTotal1]=in.readLine();
            sTotal1++;
    }*/



    //file.close();



//    file2.open(QIODevice::ReadOnly);
//    std::string line2[1000];
//    QTextStream in2(&file2);
//            while(!in2.atEnd())
//            {
//                line2[sTotal2]=in2.readLine();
//                    sTotal2++;
//            }
//    file2.close();
//
//    if(sTotal1<=n && sTotal2<=n)
//    {
//        x=0;
//        std::cout<<x<<" z true dla 1.";
//    }
//    else
//    {
//        x=1;
//        std::cout<<x<<" z false dla 1.";
//    };
//std::cout<<sTotal1;
//std::cout<<sTotal2;


 /*   if(sTotal1!=sTotal2)
    {
        x= 0;
        std::cout<<x<<" z true dla 2.";
    }
    else
    {
     x=1;
     std::cout<<x<<" z false dla 2.";
     }


    if(ui->checkBox->isChecked()==false)
    {
        file3.open(QIODevice::ReadOnly);
        std::string line3[1000];
        QTextStream in3(&file3);
                while(!in3.atEnd())
                {
                    line3[sTotal3]=in3.readLine();
                        sTotal3++;
                }
        file3.close();*/

    /*    if(x==1)
        {
            if(sTotal1!=sTotal3)
            {
                x=0;
            }

            std::cout<<x<<" z true dla 2.";
        }
        else
        {
         x=1;
         std::cout<<x<<" z false dla 2.";
         }


    }*/


//
//std::cout<<x<<"po warunkach";
//
//        if(x==1)
//        {
//
//            QMessageBox::information(this,tr("Wczytanie danych"),tr("Wczytanie danych powiodło sie"));
//        }
//        else
//        {
//            QMessageBox::critical(this,tr("Wczytanie danych"),tr("Wczytanie danych NIE powiodło sie"));
//        }
}


/*
	#TODO fix Radiobutons to bool/false values
*/
double ReservoirModel::z_func(double t_pr,double p_pr,double y)
{
    //if(ui->Z_HaY_RadioButton->isChecked())
    if(true)
	{
    double t=1/t_pr;

    return (-1)*(0.0615*t*exp(-1.2*(1-t)*(1-t)))*p_pr+((y+y*y+pow(y,3)-pow(y,4))/pow((1-y),3))
                 -(14.76*t-9.76*pow(t,2)+4.58*pow(t,3))*pow(y,2)+(90.7*t-242.2*pow(t,2)+42.4*pow(t,3))*pow(y,(2.18+2.82*t));
    }
    else if(false/*ui->Z_SRK_RadioButton->isChecked()*/)
    {

        SRK_parameters(t_pr,p_pr);
        return pow(y,3)-pow(y,2)+(A-B-B*B)*y -A*B;
    }
    else if(false /*ui->Z_RK_RadioButton->isChecked()*/)
    {
        RK_parameters(t_pr,p_pr);
        return pow(y,3)-pow(y,2)+(A-B-B*B)*y -A*B;
    }

}

double ReservoirModel::z_func_der(double t_pr,double p_pr,double y)
{
    if(true/*ui->Z_HaY_RadioButton->isChecked()*/)
    {
        double t=1/t_pr;
        return ((1+4*y+4*y*y-4*pow(y,3)+pow(y,4))/pow((1-y),4))
        -(29.52*t-19.52*pow(t,2)+9.16*pow(t,3))*y+(2.18+2.82*t)*(90.7*t-242.2* pow(t,2)+42.4* pow(t,3))*pow(y,(1.18+2.82*t));
    }
    else if(false/*ui->Z_SRK_RadioButton->isChecked()*/)
    {
        SRK_parameters(t_pr,p_pr);
        return 3*pow(y,2)-2*y+(A-B-B*B);
    }
    else if(false/*ui->Z_RK_RadioButton->isChecked()*/)
    {
        RK_parameters(t_pr,p_pr);
        return 3*pow(y,2)-2*y+(A-B-B*B);
    }
}

double ReservoirModel::zpar(double p_ri)
{
    double t_pr=t_r/t_pc();
    double p_pr=p_ri/p_pc();
    if(false/*ui->Z_RK_RadioButton->isChecked()*/)
    {
        RK_parameters(t_pr,p_pr);
    }
    else if(false/*ui->Z_SRK_RadioButton->isChecked()*/)
    {
        SRK_parameters(t_pr,p_pr);
    }



    double z_temp;
    double z=0;
    double y=0.01;



    //std::cout<<t_pr<<" and p ppr "<<p_pr<<std::endl;

    if(true/*ui->Z_HaY_RadioButton->isChecked()*/)
    {
        double t_pc_tmp=t_pc()-E_WaA();
        double p_pc_tmp=p_pc()*t_pc_tmp/(t_pc()-y_gas[13]*(1-y_gas[13])*E_WaA());

        t_pr=t_r/t_pc_tmp;
        p_pr=p_ri/p_pc_tmp;

        double t=1/t_pr;
        //std::cout<<"P_PC="<<p_pc();

        do
        {
            z_temp=z;
            y=y-z_func(t_pr,p_pr,y)/z_func_der(t_pr,p_pr,y);
            z=(0.0615*t*exp(-1.2*(1-t)*(1-t)))*p_pr/y;
        }while(fabs(z_temp-z)>0.0001);


        return z;
    }
    else
    {
		//punkty startowe, i zakres poszukiwania rozwiazania
        double a1=0.3;
        double b1=1.2;
        double x1=(a1+b1)/2;
        do
        {
                if (z_func(t_pr,p_pr,a1)*z_func(t_pr,p_pr,x1)<0) //sprawdzenie czy f(x1) znajduje sie powyzej/ponizej osi OX
                {

                    b1=x1;  //przypisanie b1 wartoœci x1 gdy f(x1)<0
                    x1=(a1+b1)/2;  //wykonanie kolejnego kroku
                }
                else
                {

                    a1=x1;  //przypisanie a1 wattoœci x1 gdy f(x1)>0
                    x1=(a1+b1)/2;   //wykonanie kolejnego kroku
                }

        }while(fabs(a1-b1)>0.00001);
        return x1;  //zwrocenie wartoœci przez funkcje
    }




}

/*!
	#TODO file input needs to be fixed
	un_i gives ammount of lines in file!

*/
//void funkcja2()
//{
//
//
//
//    std::fstream file("p_data.txt");
//
//
//
//    /*file.open(QIODevice::ReadOnly);
//    QTextStream in(&file);*/
//    std::vector<double> p;
//    int un_i = 0;
//
//
//
//        while (!in.atEnd())
//        {
//
//
//
//
//            std::string str=in.readLine();
//            double dstr=str.toDouble();
//            if(un_i==0)
//            {
//                p[0]=dstr;
//            }
//            else
//            {
//               p.append(dstr);
//            }
//
//            std::cout<<p[un_i]<<"p"<<un_i;
//
//            ++un_i;
//
//        }
//        std::cout<<"CLOSED FILE";
//        file.flush();
//        file.close();
//
//        std::fstream file2("g_p_data.txt");
//
//        file2.open(QIODevice::ReadOnly);
//        QTextStream in2(&file2);
//        std::vector<double> g_p(1);
//        un_i = 0;
//
//        while (!in2.atEnd())
//        {
//
//            std::string str2=in2.readLine();
//            double dstr2=str2.toDouble();
//            if(un_i==0)
//            {
//                g_p[0]=dstr2;
//            }
//            else
//            {
//                g_p.append(dstr2);
//            }
//
//            std::cout<<g_p[un_i]<<"g_p"<<un_i;
//
//            ++un_i;
//        }
//        file2.flush();
//        file2.close();
//        std::cout<<"CLOSED FILE";
//    std::vector<double> z(un_i+1);
//    std::cout<<"Created dynamic vector z[]-size"<<un_i;
//
//    for(int j=0;j<un_i;j++)
//    {
//        double p_el=p[j];
//        std::cout<<"j="<<j<<" p="<<p_el;
//        z[j]=zpar(p_el);
//        std::cout<<"z="<<z[j];
//    }
//
//    std::fstream file3("z_data.txt");
//    file3.open(QIODevice::ReadWrite | QIODevice::Truncate);
//    QTextStream out3(&file3);
//    for(int j=0;j<un_i;j++)
//    {
//        out3<<z[j]<<"\n";
//        std::cout<<"wpisane z["<<j<<"]";
//    }
//
//    double sum_x=0,sum_y=0,sum_x2=0,sum_xy=0,sum_y2=0;
//
//
//
//    for (int i=0;i<un_i;i++)
//    {
//
//        sum_x+=g_p[i];
//        sum_y+=p[i]/z[i];
//        sum_x2+=pow(g_p[i],2);
//        sum_y2+=pow(p[i]/z[i],2);
//        sum_xy+=g_p[i]*p[i]/z[i];
//        std::cout<<i;
//        std::cout<<"sum of x="<<sum_x<<"\taand y="<<sum_y<< std::endl;
//        std::cout<<"xy\t"<<sum_xy<<"\tx2 "<<sum_x2<< endl;
//    }
//
//
//    double a_par=a_lin(sum_x,sum_y,sum_x2,sum_xy,un_i);
//    double b_par=b_lin(sum_x,sum_y,sum_x2,sum_xy,un_i);
//    double G=(-b_par)/a_par;
//    std::cout<<"y="<<a_par<<"x+"<<b_par<< std::endl;
//    std::cout<<"G="<<G<< std::endl;
//
//    std::cout<<a_par;
//    std::cout<<b_par;
//
//    for(int k=0;k<un_i;k++)
//    {
//        std::cout<<p[k]<<"\t"<<g_p[k]<<"\t"<<z[k]<<"\t"<<p[k]/z[k]<<endl;
//    }
//
//
//
//
//    //plot p/z
//    std::vector<double> pz_plot_y;
//    for(int j=0;j<un_i;j++)
//    {
//
//        pz_plot_y[j]=p[j]/z[j];
//
//
//    }
//
//    //010000100110110001100101011100110111001101100101011001000010000001100010011001010010000001110100011010000110010100100000010011010110000101100011011010000110100101101110011001010010000001010011011100000110100101110010011010010111010000100000011101000110100001100001011101000010000001110111011010010110110001101100001000000110001101101111011011010111000001110101011101000110010100100000011101000110100001101001011100110010000001100011011011110110010001100101000011010000101001000010011011000110010101110011011100110110010101100100001000000110001001100101001000000111010001101000011001010010000001101110011000010110110101100101001000000110111101100110001000000100111101101101011011100110100101110011011100110110100101100001011010000000110100001010001100000011000100110000001100000011000100110001001100010011000100110000001100010011000100110000001100010011000100110000001100010011000000110001001100010011000000110001001100000011000000110001001100000011000100110001001100000011000100110001001100010011000000110000001100010011000100110000001100010011000000110000001100010011000000110001001100010011000100110000001100000011000100110001001100000011000100110001001100000011000100110000001100000011000100110000001100010011000100110000001100000011000000110000001100010011000000110001001100010011000000110001001100000011000000110000
//
//
//
//
//
//
//    double x_0=-b_par/a_par;
//    std::vector<double> x_linia, y_linia;
//
//    std::cout<<"y="<<a_par<<"x+"<<b_par;
//    std::cout<<"x_0="<<x_0;
//    for(int i=0;i<10;i++)
//    {
//        x_linia[i]=(x_0)*0.111*i;
//        y_linia[i]=a_par*x_linia[i]+b_par;
//        std::cout<<x_linia[i]<<"\t"<<y_linia[i];
//    }
//
//
//   /* ui->pz_plot_CustomPlot->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectPlottables);
//
//    ui->pz_plot_CustomPlot->addGraph();
//    ui->pz_plot_CustomPlot->graph(0)->setData(g_p,pz_plot_y);
//    ui->pz_plot_CustomPlot->xAxis->setLabel("G_p [mld n m^3]");
//    ui->pz_plot_CustomPlot->yAxis->setLabel("p/z [bar]");
//    ui->pz_plot_CustomPlot->xAxis->setRange(0,x_0*1.1);
//    ui->pz_plot_CustomPlot->yAxis->setRange(0,pz_plot_y[un_i-1]*1.1);
//    ui->pz_plot_CustomPlot->addGraph();
//    ui->pz_plot_CustomPlot->graph(1)->setData(x_linia,y_linia);
//
//    ui->pz_plot_CustomPlot->graph(0)->setLineStyle(QCPGraph::lsNone);
//    ui->pz_plot_CustomPlot->graph(0)->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, 4));*/
//
//
//
//
//    //ui->pz_plot_CustomPlot->replot();
//
//
//
//    //std::string msgG;
//    //msgG="liczba punktów: ";
//    //msgG+=std::string::number(un_i);
//    //msgG+="\nG=";
//    //msgG+=std::string::number(G, 'f', 5);
//    //msgG+=std::string("%1").arg(G_units());
//    //msgG+=std::string("\nWzór prostej:");
//    //msgG+=std::string("\ny=%1x+%2").arg(std::string::number(a_par),std::string::number(b_par));
//    //double S_a=pow((un_i/(un_i-2)*(sum_y2-a_par*sum_xy-b_par*sum_y)/(un_i*sum_x2-pow(sum_x,2))),0.5);
//    //double S_b=pow(((pow(S_a,2)*sum_x2)/un_i),0.5);
//    //std::cout<<"S_b="<<S_b;
//    //msgG+=std::string("\nS_a -odch. st. dla parametru a\n S_b -odch. st. dla parametru b\n");
//
//    //msgG+=std::string("\nS_a=%1\tS_b=%2").arg(std::string::number(S_a, 'f', 3),std::string::number(S_b, 'f', 3));
//    //msgG+=std::string("\nRozpiętość danych:\np=%1 bar\tG_p=%2 mld m^3").arg((std::string::number(p[0]-p[un_i-1], 'f', 3)),(std::string::number(g_p[un_i-1]-g_p[0],'f',3)));
//    //ui->label->setText(msgG);
//
//}



/*!
	This is the candidate for calcfitness()
	Non-volumetric ; needs to be fixed !
*/
double ReservoirModel::funkcja3(std::array<int,4>& intParams,std::array<float,6> &floatParams) //non-vol obliczenia
{
	
	//intParams {h,k,r_res,r_aquifer}
	//floatPrams {por,u_w,c_w,c_f,s_w}
	

    non_vol(intParams, floatParams); //wprowadza zmiany w parametrach


                            //p tworzy dynamiczna tablice p g_p oraz z
                            // czesc z obliczeniem z jest niepotrzebna narazie
    //std::fstream file("p_data.txt");
	


    //file.open(QIODevice::ReadOnly);
    //QTextStream in(&file);
	
	int un_i = p.size();

	//std::cout << "un_i" << un_i << std::endl;

	//#TODO data input from files?
        //while (!in.atEnd())
        //{




        //    std::string str=in.readLine();
        //    double dstr=str.toDouble();
        //    if(un_i==0)
        //    {
        //        p[0]=dstr;
        //    }
        //    else
        //    {
        //       p.append(dstr);
        //    }

        //    std::cout<<p[un_i]<<"p"<<un_i;

        //    ++un_i;

        //}
        //std::cout<<"CLOSED FILE";
        //file.flush();
        //file.close();

        //std::fstream file2("g_p_data.txt");

        //file2.open(QIODevice::ReadOnly);
        //QTextStream in2(&file2);
        //std::vector<double> g_p(1);
        //un_i = 0;


        //while (!in2.atEnd())
        //{

        //    std::string str2=in2.readLine();
        //    double dstr2=str2.toDouble();
        //    if(un_i==0)
        //    {
        //        g_p[0]=dstr2;
        //    }
        //    else
        //    {
        //        g_p.append(dstr2);
        //    }

        //    std::cout<<g_p[un_i]<<"g_p"<<un_i;

        //    ++un_i;
        //}
        //file2.flush();
        //file2.close();
        //std::cout<<"CLOSED FILE";

        //std::fstream filetd("time_data.txt");

        //filetd.open(QIODevice::ReadOnly);
        //QTextStream intd(&filetd);
        //std::vector<double> t(1);
        //un_i = 0;


        //while (!intd.atEnd())
        //{

        //    std::string strtd=intd.readLine();
        //    double dstrtd=strtd.toDouble();
        //    if(un_i==0)
        //    {
        //        t[0]=dstrtd;
        //    }
        //    else
        //    {
        //        t.append(dstrtd);
        //    }

        //    std::cout<<t[un_i]<<"t"<<un_i;

        //    ++un_i;
        //}
        //filetd.flush();
        //filetd.close();
        //std::cout<<"CLOSED FILE";


        //std::fstream file4("w_p_data.txt");
        //file4.open(QIODevice::ReadOnly);
        //QTextStream int4(&file4);
        //std::vector<double> w_p(1);
        //un_i = 0;


        //while (!int4.atEnd())
        //{

        //    std::string strtd=int4.readLine();
        //    double dstrtd=strtd.toDouble();
        //    if(un_i==0)
        //    {
        //        w_p[0]=dstrtd;
        //    }
        //    else
        //    {
        //        w_p.append(dstrtd);
        //    }

        //    std::cout<<w_p[un_i]<<"w_p"<<un_i;

        //    ++un_i;
        //}
        //file4.flush();
        //file4.close();
        //std::cout<<"CLOSED FILE";


	//SOME MATH, LEAVE IT
    std::vector<double> z;
	z.resize(un_i);
    //std::cout<<"Created dynamic vector z[]-size"<<un_i;

    std::vector<double> p_delta;
	p_delta.resize(un_i);
    //std::cout<<"Created dynamic vector p_delta[]-size"<<un_i;

    for(int j=0;j<un_i;j++)
    {
        if(j==0)
        {
            p_delta[j]=0;
  //          std::cout<<"p_delta przy j==0";

        }
        else if(j==1)
        {
            p_delta[j]=(p[0]-p[j])/2;
  //          std::cout<<"p_delta przy j==1";
        }
        else if(j==2)
        {
            p_delta[j]=(p[0]-p[j])/2;
 //           std::cout<<"p_delta przy j==2";
        }
        else
        {
            p_delta[j]=(p[j-2]-p[j])/2;
  //          std::cout<<"wynik="<<p[j-2]<<"-"<<p[j];
        }
   //     std::cout<<"p_delta"<<j<<" ="<<p_delta[j];
    }



    for(int j=0;j<un_i;j++)
    {
        w_p[j]=0;
    }


    //Obliczanie We
	//#TODO Add control to calculation
    std::vector<double> W_e;
	W_e.resize(un_i);
    if(true/*ui->Hurst_We_RadioButton->isChecked()*/)
    {
        for(int j=0;j<un_i;j++)
        {
            if(j==0)
            {
                W_e[j]=0;
               std::cout<<"W_e"<<j<<" ="<<W_e[j] << std::endl;
            }
            else
            {
                double sum_aquifer=0;
                for(int i=1;i<=j;i++)
                {
                    sum_aquifer=sum_aquifer+p_delta[i]*pow(10,5)*WeD(t[j]-t[i-1],r_ed);
                  //  std::cout<<"sum_aquifer="<<sum_aquifer;
                }

                W_e[j]=U(r_res,h,f,por,(c_w+c_f))*sum_aquifer;
             //   std::cout<<"W_e["<<j<<"]="<<W_e[j];
            }

        }


    }
    else if(false/*ui->Fetkovich_We_RadioButton->isChecked()*/)
    {
        std::cout<<"WYBRANO FETKOVICHA";
        double J_Fetk=2*pi*k*h*f/(mi_w*(log(r_aquifer/r_res)-3.0/4));
        double W_i=pi*(pow(r_aquifer,2)-pow(r_res,2))*h*por*f;
        double W_ei=(c_f+c_w)*W_i*p[0]*pow(10,5);
        std::vector<double> p_sr(un_i), p_a(un_i); //Pa !
        double deltaW_e;
        std::vector<double> deltaW_e_wetor(un_i);

		std::cout << "J=" << J_Fetk << "\tW_ei=" << W_ei << std::endl;
        std::cout<<"ln(ra/r_res)-3/4="<<(log(r_aquifer/r_res)-3.0/4)<<"\tr_aq/r_res="<<r_aquifer/r_res << std::endl;
        std::cout<<"gora + u_w="<<2*pi*k*h*f/(mi_w) << std::endl;

        p_sr[0]=p[0]*pow(10,5); //Pa
        W_e[0]=0;
        p_a[0]=p[0]*pow(10,5); //Pa
        std::cout<<"p_sr["<<0<<"]="<<p_sr[0] << std::endl;
        for(int i=1;i<un_i;i++)
        {
            p_sr[i]=(p[i-1]+p[i])/2*pow(10,5); //Pa
            std::cout<<"p_sr["<<i<<"]="<<p_sr[i];


        }


        for(int i=1;i<un_i;i++)
        {
           p_a[i-1]=p[0]*(1-W_e[i-1]/W_ei)*pow(10,5);


               deltaW_e=(W_ei/p[0]/pow(10,5))*(p_a[i-1]-p_sr[i])*(1-exp(-J_Fetk*p[0]*pow(10,5)*(t[i]-t[i-1])*year/W_ei));
               std::cout<<"deltaW_e"<<deltaW_e;
               deltaW_e_wetor[i]=deltaW_e;
               for(int j=0;j<=i;j++)
               {
                   W_e[i]=W_e[i]+deltaW_e_wetor[j];
               }

           std::cout<<"W_e["<<i<<"]="<<W_e[i]<<"p_a["<<i-1<<"]="<<p_a[i-1]<<"\tFetkovich";
        }








    }

    double We_sum=0;

	//OUTPUT TO FILE

	//std::fstream fileW("W_e_data.txt");
    //fileW.open(QIODevice::ReadWrite | QIODevice::Truncate);
    //QTextStream outW(&fileW);
    //for(int j=0;j<un_i;j++) //plus suma
    //{
    //    We_sum+=W_e[j];
    //    outW<<W_e[j]<<"\n";
    //    std::cout<<"wpisane W_e["<<j<<"]";
    //}

    //fileW.close();

    // oblicza z na podstawie tablicy p
    //// TRZEBA ZMIENIC (?)
    for(int j=0;j<un_i;j++)
    {
        double p_el=p[j];
    //    std::cout<<"j="<<j<<" p="<<p_el;
        z[j]=zpar(p_el);
    //    std::cout<<"z="<<z[j];
    }

  //  std::cout<<"un_i="<<un_i;


    for(int i=0;i<un_i;i++)
    {
  //      std::cout<<"z["<<i<<"]="<<z[i];
    }

	//OUTPUT TO FILE

    /*std::fstream file3("z_data.txt");
    file3.open(QIODevice::ReadWrite | QIODevice::Truncate);
    QTextStream out3(&file3);
    for(int j=0;j<un_i;j++)
    {
        out3<<z[j]<<"\n";
        std::cout<<"wpisane z["<<j<<"]";
    }
    file3.close();*/

    ////SUMA Doplywow We
 //   std::cout<<"Suma We="<<We_sum;


    //Obliczenia do wykresu

    double Bg_i=Bg(z[0],p[0],t_r);
    // F/E_g
    std::vector<double> F;
    std::vector<double> E_g;
	std::vector<double> E_fw;
    std::vector<double> drive_plot_y, pz_plot_y, match_plot_x, error_y,error_x;
	F.resize(un_i);
	E_g.resize(un_i);
	E_fw.resize(un_i);
	drive_plot_y.resize(un_i);
	match_plot_x.resize(un_i);

   /* std::fstream fileF("F_results.txt");
    std::fstream fileE_g("E_g_results.txt");
    std::fstream fileE_fw("E_fw_results.txt");
*/
    //fileF.open(QIODevice::ReadWrite | QIODevice::Truncate);
    //fileE_g.open(QIODevice::ReadWrite | QIODevice::Truncate);
    //fileE_fw.open(QIODevice::ReadWrite | QIODevice::Truncate);

    //QTextStream outFres(&fileF);
    //QTextStream outEgres(&fileE_g);
    //QTextStream outEfwres(&fileE_fw);



    for(int j=0;j<un_i;j++)
    {
        F[j]=F_HO(g_p[j],z[j],p[j],w_p[j]);
        E_g[j]=Bg(z[j],p[j],t_r)-Bg_i;
        E_fw[j]=Efw(p[j],c_w,s_w,c_f,z[0],p[0],t_r);
        //outFres<<F[j]<<"\n";
        //outEgres<<E_g[j]<<"\n";\
        //outEfwres<<E_fw[j]<<"\n";

		//#TODO THIS DATA NEEDED FOR FITNESS
        drive_plot_y[j]=F[j]/(E_g[j]+E_fw[j]);
        //pz_plot_y[j]=p[j]/z[j];
        match_plot_x[j]=W_e[j]*B_w(p[j],t_r)/(E_g[j]+E_fw[j])/pow(10,9); 
        //error_y[j]=drive_plot_y[j]*0.001;
        //error_x[j]=match_plot_x[j]*0.020;
      // std::cout<<"\nWyniki\n"<<F[j]<<"\t"<<E_g[j]<<"\t"<<E_fw[j]<<"\t"<<"\nplot\n"<<drive_plot_y[j]<<"\t"<<match_plot_x[j]<<"\\ Bg "<<Bg(z[j],p[j],t_r) <<"j="<<j << std::endl;
    }

	double sum_x = 0, sum_y = 0, sum_x2 = 0, sum_xy = 0, sum_y2 = 0;
	for (int i = 1; i<un_i; i++)
	{
		//10e-6 due to billion m^3 to m^3 rate
		sum_x += match_plot_x[i];
		sum_y += drive_plot_y[i];
		sum_x2 += pow(match_plot_x[i], 2);
		sum_y2 += pow(drive_plot_y[i], 2);
		sum_xy += match_plot_x[i] * drive_plot_y[i];
		//std::cout << i;
		//std::cout << "sum of x=" << sum_x << "\taand y=" << sum_y << std::endl;
		//std::cout << "xy\t" << sum_xy << "\tx2 " << sum_x2 << std::endl;
	}
	  // std::cout<<sum_x<<"\t"<<sum_y<<"\t"<<sum_x2<<"\t"<<sum_xy<<std::endl;

	double a_par = a_lin(sum_x, sum_y, sum_x2, sum_xy, un_i);
	std::cout << "\n\na_par=" << a_par << std::endl;
	
	//#TODO Sanity check for Nan - done temporary
	if (isnan(a_par))
		a_par = 0;

	return (a_par);
	//double b_par = b_lin(sum_x, sum_y, sum_x2, sum_xy, un_i);




    //fileF.close();
    //fileE_g.close();
    //fileE_fw.close();


	//GRAPH DRAWING

    /*ui->drive_plot_CustomPlot->addGraph();
    ui->drive_plot_CustomPlot->graph(0)->setData(g_p, drive_plot_y);

    ui->drive_plot_CustomPlot->xAxis->setLabel("G_p [mld nm^3]");
    ui->drive_plot_CustomPlot->yAxis->setLabel("F/(E_g+E_fw)[mld m^3]");

    ui->drive_plot_CustomPlot->xAxis->setRange(0, g_p[un_i-1]*1.2);
    ui->drive_plot_CustomPlot->yAxis->setRange(0, drive_plot_y[un_i-1]*2);

    ui->drive_plot_CustomPlot->graph(0)->setLineStyle(QCPGraph::lsLine);
    ui->drive_plot_CustomPlot->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectPlottables);
    ui->drive_plot_CustomPlot->replot();
*/


    //p/z plot


	//SQUARE METHOD FOR LINEAR REGRESSION 

   // double sum_x=0,sum_y=0,sum_x2=0,sum_xy=0,sum_y2=0;





//    for (int i=0;i<un_i;i++)
//    {
//
//        sum_x+=g_p[i];
//        sum_y+=p[i]/z[i];
//        sum_x2+=pow(g_p[i],2);
//        sum_y2+=pow(p[i]/z[i],2);
//        sum_xy+=g_p[i]*p[i]/z[i];
//        std::cout<<i;
//        std::cout<<"sum of x="<<sum_x<<"\taand y="<<sum_y<< std::endl;
//        std::cout<<"xy\t"<<sum_xy<<"\tx2 "<<sum_x2<< std::endl;
//    }
////    cout<<sum_x<<"\t"<<sum_y<<"\t"<<sum_x2<<"\t"<<sum_xy<<endl;
//
//    double a_par=a_lin(sum_x,sum_y,sum_x2,sum_xy,un_i); 
//    double b_par=b_lin(sum_x,sum_y,sum_x2,sum_xy,un_i);
//    //y=ax+b
//    //x_0=-b/a
//    double x_0=-b_par/a_par;
//    std::vector<double> x_linia(10), y_linia(10);
//
//    std::cout<<"y="<<a_par<<"x+"<<b_par;
//    std::cout<<"x_0="<<x_0;
//    for(int i=0;i<10;i++)
//    {
//        x_linia[i]=(x_0)*0.111*i;
//        y_linia[i]=a_par*x_linia[i]+b_par;
//    }



   /* ui->pz_plot_CustomPlot->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectPlottables);

    ui->pz_plot_CustomPlot->addGraph();
    ui->pz_plot_CustomPlot->graph(0)->setData(g_p,pz_plot_y);
    ui->pz_plot_CustomPlot->xAxis->setLabel("G_p [mld nm^3]");
    ui->pz_plot_CustomPlot->yAxis->setLabel("p/z [bar]");
    ui->pz_plot_CustomPlot->xAxis->setRange(0,x_0*1.1);
    ui->pz_plot_CustomPlot->yAxis->setRange(0,pz_plot_y[un_i-1]*1.1);
    ui->pz_plot_CustomPlot->addGraph();
    ui->pz_plot_CustomPlot->graph(1)->setData(x_linia,y_linia);
    ui->pz_plot_CustomPlot->graph(1)->rescaleAxes();
    ui->pz_plot_CustomPlot->graph(0)->setLineStyle(QCPGraph::lsNone);
    ui->pz_plot_CustomPlot->graph(0)->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, 4));

    ui->pz_plot_CustomPlot->replot();*/




    //match plot


    //straigh lines
    //std::vector<double> x1(2),y1(2),x2(2),y2(2),x3(2),y3(2),x4(2),y4(2);


    //    x1[0]=match_plot_x[1];
    //    x1[1]=match_plot_x[un_i-1];
    //    y1[0]=x1[0]+drive_plot_y[1];
    //    y1[1]=x1[1]+drive_plot_y[1];

    //    x2[0]=match_plot_x[1];
    //    x2[1]=match_plot_x[un_i-1];
    //    y2[0]=x2[0]+drive_plot_y[1]*1.2;
    //    y2[1]=x2[1]+drive_plot_y[1]*1.2;

    //    x3[0]=match_plot_x[1];
    //    x3[1]=match_plot_x[un_i-1];
    //    y3[0]=x3[0]+drive_plot_y[1]*0.8;
    //    y3[1]=x3[1]+drive_plot_y[1]*0.8;



   /* ui->match_plot_CustomPlot->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectPlottables);

    ui->match_plot_CustomPlot->addGraph();
    ui->match_plot_CustomPlot->graph(0)->setData(match_plot_x,drive_plot_y);
    ui->match_plot_CustomPlot->xAxis->setLabel("We*Bw/(Eg+Efw [mld m^3])");
    ui->match_plot_CustomPlot->yAxis->setLabel("F/(Eg+Efw) [mld m^3]");
    ui->match_plot_CustomPlot->graph(0)->rescaleAxes();
    ui->match_plot_CustomPlot->graph(0)->setErrorType(QCPGraph::etBoth);
    ui->match_plot_CustomPlot->graph(0)->setDataBothError(match_plot_x, drive_plot_y,error_x, error_y);
    ui->match_plot_CustomPlot->graph(0)->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, 4));
    ui->match_plot_CustomPlot->graph(0)->setPen(QPen(Qt::red));*/

    /*QPen blueDotPen;
    blueDotPen.setColor(QColor(30, 40, 255, 150));
    blueDotPen.setStyle(Qt::DotLine);
    blueDotPen.setWidthF(2);

    ui->match_plot_CustomPlot->addGraph();
    ui->match_plot_CustomPlot->graph(1)->setData(x1,y1);
    ui->match_plot_CustomPlot->graph(1)->setPen(blueDotPen);
    ui->match_plot_CustomPlot->addGraph();
    ui->match_plot_CustomPlot->graph(2)->setData(x2,y2);
    ui->match_plot_CustomPlot->graph(2)->setPen(blueDotPen);
    ui->match_plot_CustomPlot->addGraph();
    ui->match_plot_CustomPlot->graph(3)->setData(x3,y3);
    ui->match_plot_CustomPlot->graph(3)->setPen(blueDotPen);

    ui->match_plot_CustomPlot->replot();*/



    //std::string non_volmsg="Obliczenia dla warunków niewolumetrycznych:\nZakończone\n";



    //non_volmsg+=std::string("Rozpiętość danych:\np=%1 bar\tG_p=%2 mld m^3").arg((std::string::number(p[0]-p[un_i-1], 'f', 3)),(std::string::number(g_p[un_i-1]-g_p[0],'f',3)));

    //ui->label->setText(non_volmsg);



}

/*!
	Reads data from file to textboxes
*/
//void ReservoirModel::on_pushButton_clicked()
//{
//
//   std::fstream file("p_data.txt");
//   if(!file.open(QIODevice::ReadWrite))
//           QMessageBox::information(0,"Error opening file",file.errorString());
//   QTextStream in(&file);
//   ui->textEdit->setText(in.readAll());
//
//
//
//    std::fstream file2("g_p_data.txt");
//    if(!file2.open(QIODevice::ReadWrite))
//            QMessageBox::information(0,"Error opening file",file2.errorString());
//    QTextStream in2(&file2);
//    ui->textEdit_2->setText(in2.readAll());
//
//    if(ui->checkBox->isChecked()==false)
//    {
//        std::fstream file3("time_data.txt");
//        if(!file3.open(QIODevice::ReadWrite))
//                QMessageBox::information(0,"Error opening file",file3.errorString());
//        QTextStream in3(&file3);
//        ui->textEdit_3->setText(in3.readAll());
//
//
//        std::fstream file4("w_p_data.txt");
//        if(!file4.open(QIODevice::ReadWrite))
//                QMessageBox::information(0,"Error opening file",file4.errorString());
//        QTextStream in4(&file4);
//        ui->textEdit_4->setText(in4.readAll());
//    }
//
//
//
//}

/*!
	Saves data to file 
*/
//void ReservoirModel::on_pushButton_2_clicked()
//{
//    funkcja1();
//
//
//    std::fstream outfile;
//    outfile.setFileName("p_data.txt");
//    outfile.open(QIODevice::ReadWrite | QIODevice::Truncate);
//    QTextStream out(&outfile);
//    out << ui->textEdit->toPlainText() << "\n";
//
//    outfile.close();
//
//    std::fstream outfile2;
//    outfile2.setFileName("g_p_data.txt");
//    outfile2.open(QIODevice::ReadWrite | QIODevice::Truncate);
//    QTextStream out2(&outfile2);
//    out2 << ui->textEdit_2->toPlainText() << endl;
//
//    outfile2.close();
//
//    if(ui->checkBox->isChecked()==false) //non vol
//    {
//        std::fstream outfile3;
//        outfile3.setFileName("time_data.txt");
//        outfile3.open(QIODevice::ReadWrite | QIODevice::Truncate);
//        QTextStream out3(&outfile3);
//        out3 << ui->textEdit_3->toPlainText() << endl;
//
//        std::fstream outfile4;
//        outfile4.setFileName("w_p_data.txt");
//        outfile4.open(QIODevice::ReadWrite | QIODevice::Truncate);
//        QTextStream out4(&outfile4);
//        out4 << ui->textEdit_4->toPlainText() << endl;
//    }
//
//
//
//
//
//
//}

double ReservoirModel::doCalculcations(std::array<int, 4> &ar1,std::array<float, 6> &ar2)
{
    /*ui->z_textBrowser->setText("");
    ui->W_e_textBrowser->setText("");
    ui->label->setText("");*/
    //if(fy_gas()==!100)
    //{
    //    QMessageBox ::information(0,"Obliczenia przerwane","Skład gazu nieprawidłowy");
    //}
        ;
    temperature();
    if(false/*ui->checkBox->isChecked()==true*/) //true = war wolumetryczne
    {

        //funkcja2();

        //std::fstream file3("z_data.txt");
        //if(!file3.open(QIODevice::ReadOnly))
        //        QMessageBox ::information(0,"Error opening file",file3.errorString());
        //QTextStream in3(&file3);
        //ui->z_textBrowser->setText(in3.readAll());

        //file3.close();
    }
    else
    {
		//#TODO ARRAYS WITH DATA

        return funkcja3(ar1,ar2);

		//EXTRACTING RESULTS FROM FILE (funkcja3() saved those to file)

        /*std::fstream fileW("W_e_data.txt");
        if(!fileW.open(QIODevice::ReadOnly))
                QMessageBox ::information(0,"Error opening file",fileW.errorString());
        QTextStream inW(&fileW);
        ui->W_e_textBrowser->setText(inW.readAll());

        fileW.close();

        std::fstream file3("z_data.txt");
        if(!file3.open(QIODevice::ReadOnly))
                QMessageBox ::information(0,"Error opening file",file3.errorString());
        QTextStream in3(&file3);
        ui->z_textBrowser->setText(in3.readAll());

        file3.close();

        std::fstream fileF("F_results.txt");
        std::fstream fileE_g("E_g_results.txt");
        std::fstream fileE_fw("E_fw_results.txt");

        if(!fileF.open(QIODevice::ReadOnly))
                QMessageBox ::information(0,"Error opening file",fileF.errorString());

        if(!fileE_g.open(QIODevice::ReadOnly))
                QMessageBox ::information(0,"Error opening file",fileE_g.errorString());

        if(!fileE_fw.open(QIODevice::ReadOnly))
                QMessageBox ::information(0,"Error opening file",fileE_fw.errorString());


        QTextStream inp1(&fileF);
        ui->F_textBrowser->setText(inp1.readAll());
        fileF.close();
        std::cout<<"Wypisano F";

        QTextStream inp2(&fileE_g);
        ui->E_g_textBrowser->setText(inp2.readAll());
        fileE_g.close();
        std::cout<<"Wypisano Eg";

        QTextStream inp3(&fileE_fw);
        ui->E_fw_textBrowser->setText(inp3.readAll());
        fileE_fw.close();
        std::cout<<"Wypisano Efw";*/


    }




}

/*
	sets global t_r value 
	t_r - reservoir
	#TODO set ssetter
*/
void ReservoirModel::temperature() //do globalnej zmiennej temp_zlozowej
{
    //t_r=(ui->temperatureSpinBox->value());
	t_r = 330; //<-sample value!
}


/*
	#TODO this where input parameters are
*/
void ReservoirModel::non_vol(std::array<int, 4>& intParams, std::array<float, 6> &floatParams) //wprowadzenie zmian w parametrach
{

	h = intParams.at(0);
	k = intParams.at(1)*0.986923*pow(10, -15);
	r_res = intParams.at(2);
	r_aquifer = intParams.at(3);
	
	por = floatParams.at(0);
	mi_w = floatParams.at(1)*0.001;
	c_w = floatParams.at(2)*pow(10, -9);
	c_f = floatParams.at(3)*pow(10, -9);
	f = floatParams.at(4);
	s_w = floatParams.at(5);
	
    //h=ui->hDoubleSpinBox->value();
    //por=ui->porosityDoubleSpinBox->value();
    //k=ui->perm_absDoubleSpinBox->value()*0.986923*pow(10,-15);
    //mi_w=ui->water_viscDoubleSpinBox->value()*0.001;
    //c_w=ui->compr_waterDoubleSpinBox->value()*pow(10,-9);
    //c_f=ui->compr_formationDoubleSpinBox->value()*pow(10,-9);
    //f=ui->fDoubleSpinBox->value();
    //r_res=ui->r_reservoirDoubleSpinBox->value();
    //r_aquifer=ui->r_aquiferDoubleSpinBox->value();
    //r_ed=r_aquifer/r_res;
    //s_w=ui->s_wDoubleSpinBox_2->value(); //dunno why _2 - KIEDYS skasowalem i dodalem ponownie wszystkie SpinBoxy


}




ReservoirModel::ReservoirModel():
	h(20),por(0.20),k(15 * 0.986923*pow(0.1, 15)),mi_w(1.0),c_w(0.5e-9),c_f(0.6e-9),f(0.75),r_res(300.0),r_aquifer(700.0),t_r(300),r_ed(r_aquifer/r_res),s_w(0.23)
{
}

void ReservoirModel::setData(std::vector<double> &pressure, std::vector<double> &gas_p, std::vector<double> &time, std::vector<double> &water_p)
{
	if (p.size()==g_p.size() && t.size()==w_p.size() && p.size()==t.size())
	{
		p= pressure;
		g_p = gas_p;
		w_p = water_p;
		t = time;
		std::cout << "assigned vectors ,size=" << pressure.size() << "\t p.size=" << p.size()<<std::endl;
	}
}

ReservoirModel::~ReservoirModel()
{
}

/*only checks if gas properties are ok*/
int ReservoirModel::fy_gas()
{
    //y_gas[0]=ui->y_gas_1->value();
    //y_gas[1]=ui->y_gas_2->value();
    //y_gas[2]=ui->y_gas_3->value();
    //y_gas[3]=ui->y_gas_4->value();
    //y_gas[4]=ui->y_gas_5->value();
    //y_gas[5]=ui->y_gas_6->value();
    //y_gas[6]=ui->y_gas_7->value();
    //y_gas[7]=ui->y_gas_8->value();
    //y_gas[8]=ui->y_gas_9->value();
    //y_gas[9]=ui->y_gas_10->value();
    //y_gas[10]=ui->y_gas_11->value();
    //y_gas[11]=ui->y_gas_12->value();
    //y_gas[12]=ui->y_gas_13->value();
    //y_gas[13]=ui->y_gas_14->value();
    //y_gas[14]=ui->y_gas_15->value();



    //int sum=0;
    //for(int k=0;k<15;k++)
    //{

    //     sum=sum+round(y_gas[k]*100);
    //    std::cout<<"y_gas= "<<sum;
    //}
    //if(sum!=100)
    //{
    //    //QMessageBox::critical(this,tr("Wczytanie danych"),tr("Suma udzialow rozna od 1"));
    //}
	int sum = 100;

    return sum;

}



//void ReservoirModel::on_checkBox_toggled(bool checked)
//{
//    std::string chkres;
//    bool chk= ui->checkBox->checkState();
//    if(chk==true)
//    {
//        chkres="warunki wolumetryczne";
//    }
//    else
//    {
//         chkres="warunki niewolumetryczne";
//    }
//    QMessageBox ::information(0,tr("Wybrano"),chkres);
//
//
//}


/*
	#TODO Data/file input
	un_i for looping needs to be known!
	Probably this method is NOT needed!
*/
//void ReservoirModel::on_pushButton_5_clicked()
//{
//    fy_gas();
//    non_vol();
//    const int day=24*3600;
//    const double U=2*3.1415*r_res*r_res*h*por*(c_w+c_f)*f;
//    double G_i=ui->G_i_DoubleSpinBox->value()*pow(10,9);
//    int TOL=ui->tol_DoubleSpinBox->value();
//    long time=ui->time_DoubleSpinBox->value()*year; //SpinBox - years to sec
//    double W_e=WeD(time/year,r_ed);
//    double q=ui->q_DoubleSpinBox->value()*1000/(3600*24);
//    double G_p=q*time;
//    double p_initial;
//
//
//    std::vector<double> p_temp(20),We(20), B_g(20);
//
//    std::fstream file("p_data.txt");
//    file.open(QIODevice::ReadOnly);
//    QTextStream in(&file);
//    std::string str=in.readLine();
//    double dstr=str.toDouble();
//    p_initial=dstr;
//    file.close();
//    std::cout<<"p_initial="<<p_initial;
//
//    double z_initial=zpar(p_initial);
//    double pz_initial=p_initial/z_initial;
//    std::cout<<"pz="<<pz_initial;
//    p_temp[0]=pz_initial*(1-G_p/G_i);
//    std::cout<<"1-Gp/G="<<(1-G_p/G_i);
//    std::cout<<"p_temp="<<p_temp[0]<<"\t"<<0;
//    We[0]=0;
//    We[1]=U*WeD(time/year,r_ed)*(p_initial-p_temp[0]);
//    double Bgi=Bg(z_initial,p_initial,t_r);
//    B_g[0]=Bgi;
//    float difference=0;
//    int i=0;
//    do
//    {
//        i++;
//        p_temp[i]=pz_initial*(1-G_p/G_i)/(1-We[i]*(fabs(Bgi-B_g[i-1]))/G_i)*zpar(p_temp[i-1]);
//        std::cout<<"p_temp="<<p_temp[i]<<"\t"<<i;
//        We[i]=U*WeD(time/year,r_ed)*(p_initial-p_temp[i]);
//        B_g[i]=Bg(zpar(p_temp[i]),p_temp[i],t_r);
//        difference=fabs(p_temp[i]-p_temp[i-1]);
//        std::cout<<"difference="<<difference;
//
//
//    }while(difference>TOL 	&& i<19);
//    std::cout<<p_temp[i];
//    ui->p_res_TextLabel->setText(std::string("p=%1 bar").arg(std::string::number(p_temp[i])));
//
//    if(ui->Fetkovich_We_RadioButton->isChecked())
//    {
//        std::fstream file("p_data.txt");
//        std::fstream file2("g_p_data.txt");
//        std::fstream file3("z_data.txt");
//        std::fstream file4("time_data.txt");
//        std::fstream file5("W_e_data.txt");
//		std::vector<double> p;
//		std::vector<double> g_p;
//		std::vector<double> z;
//		std::vector<double> t;
//		std::vector<double> W_e;
//
//        int un_i = 0;
//		
//		//FILE INPUT
//
//        /*file.open(QIODevice::ReadOnly);
//        QTextStream in(&file);
//        
//
//        file2.open(QIODevice::ReadOnly);
//        QTextStream in2(&file2);
//        
//
//        file3.open(QIODevice::ReadOnly);
//        QTextStream in3(&file3);
//        
//
//        file4.open(QIODevice::ReadOnly);
//        QTextStream in4(&file4);
//        
//
//        file5.open(QIODevice::ReadOnly);
//        QTextStream in5(&file5);
//        
//
//
//
//            while (!in.atEnd())
//            {
//
//
//
//
//                std::string str=in.readLine();
//                double dstr=str.toDouble();
//                std::string str2=in2.readLine();
//                double dstr2=str2.toDouble();
//                std::string str3=in3.readLine();
//                double dstr3=str3.toDouble();
//                std::string str4=in4.readLine();
//                double dstr4=str4.toDouble();
//                std::string str5=in5.readLine();
//                double dstr5=str5.toDouble();
//
//                if(un_i==0)
//                {
//                    p[0]=dstr;
//                    g_p[0]=dstr2;
//                    z[0]=dstr3;
//                    t[0]=dstr4;
//                    W_e[0]=dstr5;
//                }
//                else
//                {
//                   p.append(dstr);
//                   g_p.append(dstr2);
//                   z.append(dstr3);
//                   t.append(dstr4);
//                   t.append(dstr5);
//                }
//
//                std::cout<<p[un_i]<<"p"<<un_i;
//
//                ++un_i;
//
//            }
//
//            file.flush();
//            file.close();
//            std::cout<<"CLOSED FILE";
//            file2.flush();
//            file2.close();
//            std::cout<<"CLOSED FILE";
//            file3.flush();
//            file3.close();
//            std::cout<<"CLOSED FILE";
//            file4.flush();
//            file4.close();
//            std::cout<<"CLOSED FILE";
//            file5.flush();
//            file5.close();
//            std::cout<<"CLOSED FILE";
//*/
//
//            double G_i=ui->G_i_DoubleSpinBox->value()*pow(10,9);
//            int TOL=ui->tol_DoubleSpinBox->value();
//            long time=ui->time_DoubleSpinBox->value()*year; //SpinBox - years to sec
//
//            double q=ui->q_DoubleSpinBox->value()*1000/(3600*24);
//            double G_p=q*time;
//
//
//
//            std::cout<<"WYBRANO FETKOVICHA un_i="<<un_i;
//            double J_Fetk=2*pi*k*h*f/(mi_w*(log(r_aquifer/r_res)-3.0/4));
//            double W_i=pi*(pow(r_aquifer,2)-pow(r_res,2))*h*por*f;
//            double W_ei=(c_f+c_w)*W_i*p[0]*pow(10,5);
//            std::vector<double> p_sr, p_a ; //Pa !
//            double deltaW_e;
//			std::vector<double> deltaW_e_wetor;
//
//			
//
//
//            p_sr[0]=p[0]*pow(10,5); //Pa
//            W_e[0]=0;
//            p_a[0]=p[0]*pow(10,5); //Pa
//            std::cout<<"p_sr["<<0<<"]="<<p_sr[0];
//            for(int i=1;i<un_i;i++)
//            {
//                p_sr[i]=(p[i-1]+p[i])/2*pow(10,5); //Pa
//                std::cout<<"p_sr["<<i<<"]="<<p_sr[i];
//            }
//
//
//
//            for(int i=1;i<un_i;i++)
//            {
//               p_a[i-1]=p[0]*(1-W_e[i-1]/W_ei)*pow(10,5);
//               std::cout<<"p_a["<<i-1<<"]="<<p_a[i-1];
//
//
//            }
//
//
//            double p_initial=p[0];
//            double Bgi=Bg(z_initial,p_initial,t_r);
//            double z_initial=zpar(p_initial);
//            double pz_initial=p_initial/z_initial;
//            std::vector<double> p_temp;
//			p_temp.reserve(un_i);
//
//            if(time/year>t[un_i])
//            {
//                float difference=0;
//                int i=0;
//                do
//                {
//                    i++;
//                    p_temp[un_i]=pz_initial*(1-(g_p[un_i-1]+G_p)/G_i)/(1-W_e[un_i-1]*Bgi/G_i);
//                    deltaW_e=W_ei/p_initial*(p_a[un_i-1] + p_temp[un_i])*0.5*(1-exp(-J_Fetk*p_initial*(time/year-t[un_i-1])/W_ei));
//                    difference=fabs(p_temp[i]-p_temp[i-1]);
//                    std::cout<<"difference="<<difference;
//
//                }while(difference>TOL 	&& i<19);
//            }
//			std::cout << "KURWA JEST! p=" << p_temp[i] << std::endl;
//            //ui->p_res_TextLabel->setText(std::string("p=%1 bar").arg(std::string::number(p_temp[i])));
//		
//
//
//
//
//
//
//
//
//
//        std::cout<<"J="<<J_Fetk<<"\tW_ei="<<W_ei;
//        std::cout<<"ln(ra/r_res)-3/4="<<(log(r_aquifer/r_res)-3.0/4)<<"\tr_aq/r_res="<<r_aquifer/r_res;
//        std::cout<<"gora + u_w="<<2*pi*k*h*f/(mi_w);
//		std::cout<<std::endl;
//
//
//    }
//
//
//}

/*!
	#TODO
	q5,t5  to be extracted
*/

//void ReservoirModel::on_q_DoubleSpinBox_valueChanged(double arg1)
//{
//    double q5=ui->q_DoubleSpinBox->value()*(1000.0/(3600*24));
//    int t5=ui->time_DoubleSpinBox->value()*year;
//
//    ui->G_p_LabelText->setText(std::string("G_p=%1 [mld m^3]").arg(std::string::number(q5*t5*pow(0.1,9))));
//}
//
//
//void ReservoirModel::on_time_DoubleSpinBox_valueChanged(double arg1)
//{
//    double q5=ui->q_DoubleSpinBox->value()*(1000.0/(3600*24));
//    int t5=ui->time_DoubleSpinBox->value()*year;
//
//    ui->G_p_LabelText->setText(std::string("G_p=%1 [mld m^3]").arg(std::string::number(q5*t5*pow(0.1,9))));
//
//}



