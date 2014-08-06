#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>

#include "dynamic_topography.h"
// #include "integration.h"
#include "log.h"

// #include <windows.h>


#include <iostream>
#include <sstream>
#include <locale>

using namespace std;


void print_eng_usage()
{
	cout << "This Software is for calculating dynamic topography between two points\n\n";
	cout << "USAGE:\tintegral_dt.exe <vp_out_file> <boundary_points_list> <dt_out_file>\n\n"

		 << "\t<input_file>\t\tVecPlotter output file\n"
		 	<< "\t    string format:\n"
				<< "\t\tgeo longitude start\n"
				<< "\t\tgeo latitude start\n"
				<< "\t\tgeo longitude end\n"
				<< "\t\tgeo latitude end\n"
				<< "\t\tpixel longitude start\n"
				<< "\t\tpixel latitude start\n"
				<< "\t\tpixel longitude end\n"
				<< "\t\tpixel latitude end\n"				
				<< "\t\tcorrelation\n"
				<< "\t\tvelocity of movement\n"
				<< "\t\ta priori error\n"

		 << "\t<boundary_points_list>\tList of cut definitions\n"
			<< "\t    string format:\n"
				<< "\t\tgeo longitude start\n"
				<< "\t\tgeo latitude start\n"
				<< "\t\tgeo longitude end\n"
				<< "\t\tgeo latitude end\n"
				<< "\t\tcut width (radius), [degrees]\tor -1 by default\n"
				<< "\t\tinterpolation interval (diameter), [degrees]\tor -1 by default\n"
				<< "\t\tweight coefficient\tor -1 by default\n"

		 << "\t<dt_out_file>\t\tfile for result output\n"
			<< "\t    string format:\n"
				<< "\t\tgeo longitude start\n"
				<< "\t\tgeo latitude start\n"
				<< "\t\tgeo longitude end\n"
				<< "\t\tgeo latitude end\n"
				<< "\t\tDynamic Topography, [meters]\n"
				<< "\t\tcut width (radius), [degrees]\n"
				<< "\t\tinterpolation interval (diameter), [degrees]\n"
				<< "\t\tweight coefficient\n"
				<< "\t\tintegration error (K|In-I2n|=dDT)\n"
				<< "\t\tinterpolation accuracy\n"
				<< "\t\tintegration error\n"
				<< "\t\tmean-square deviation\n"
				<< "\t\ta priori error\n"
				<< "\t\tcut length, [meters]\n"
				<< "\t\tDT coefficient\n"
				<< "\t\tintegration step size [meters]\n"
				<< "\t\tintegration step count\n"
				<< "\t\tvector count\n"
		 << "\n";
}

void print_usage()
{
	// char buf[1024], rbuf[1024];
	cout <<	 "Программное средство для расчета динамической топографии океана (моря) между двумя точками по заданному полю перемещений\n\n";
	// CharToOemBuff(buf, rbuf, sizeof(buf));
	// cout << rbuf;
	cout << "USAGE:\tintegral_dt.exe <vp_out_file> <boundary_points_list> <dt_out_file>\n\n"

		 << "\t<input_file>\t\tVecPlotter output file\n"
		 	<< "\t    string format:\n"
				<< "\t\tgeo longitude start\n"
				<< "\t\tgeo latitude start\n"
				<< "\t\tgeo longitude end\n"
				<< "\t\tgeo latitude end\n"
				<< "\t\tpixel longitude start\n"
				<< "\t\tpixel latitude start\n"
				<< "\t\tpixel longitude end\n"
				<< "\t\tpixel latitude end\n"				
				<< "\t\tcorrelation\n"
				<< "\t\tvelocity of movement\n"
				<< "\t\ta priori error\n"

		 << "\t<boundary_points_list>\tList of cut definitions\n"
			<< "\t    string format:\n"
				<< "\t\tgeo longitude start\n"
				<< "\t\tgeo latitude start\n"
				<< "\t\tgeo longitude end\n"
				<< "\t\tgeo latitude end\n"
				<< "\t\tcut width (radius), [degrees]\n"
				<< "\t\tinterpolation interval (diameter), [degrees]\n"
				<< "\t\tweight coefficient\n"

		 << "\t<dt_out_file>\t\tfile for result output\n"
			<< "\t    string format:\n"
				<< "\t\tgeo longitude start\n"
				<< "\t\tgeo latitude start\n"
				<< "\t\tgeo longitude end\n"
				<< "\t\tgeo latitude end\n"
				<< "\t\tDynamic Topography, [meters]\n"
				<< "\t\tcut width (radius), [degrees]\n"
				<< "\t\tinterpolation interval (diameter), [degrees]\n"
				<< "\t\tweight coefficient\n"
				<< "\t\tintegration error (K|In-I2n|=dDT)\n"
				<< "\t\tinterpolation accuracy\n"
				<< "\t\tintegration error\n"
				<< "\t\tmean-square deviation\n"
				<< "\t\ta priori error\n"
				<< "\t\tcut length, [meters]\n"
				<< "\t\tDT coefficient\n"
				<< "\t\tintegration step size [meters]\n"
				<< "\t\tintegration step count\n"
				<< "\t\tvector count\n"
		 << "\n";
}

char* move_points_file;
char* station_points_file;
char* out_file;
// char* output_log = (char *)"log.txt";
// char* itg_log = (char *)"itg_log.txt";

bool file_exists(const char *fname)
{
	return ifstream(fname) != NULL;
}

bool analize_options(int argc, char** argv)
{
	if (argc > 3)
	{
		move_points_file = argv[1];
		station_points_file = argv[2];
		out_file = argv[3];
		if (strcmp(move_points_file, station_points_file) && file_exists(move_points_file) && file_exists(station_points_file))
			return true;
	}
	// else
	// 	print_usage();
	return false;
}

void read_cuts(char *file_name, vector <scut> &cut)
{
	fstream fcut;
	fcut.open(file_name);
	double gsx, gsy, gex, gey, cut_width, itp_diameter, weight_coef;
	
	const char *dm = " \t";

	while(fcut >> gsx >> gsy >> gex >> gey >> cut_width >> itp_diameter >> weight_coef)
	{
		vec v(point(gsx, gsy), point(gex, gey));
		cut.push_back(scut(v, cut_width, itp_diameter, weight_coef));
	}
	fcut.close();
}

void read_movement_field(char *file_name, vector <movement> &mvn)
{
	fstream fmoves;
	fmoves.open(file_name);
	double gsx, gsy, gex, gey, crl, vlc, err;
	int psx, psy, pex, pey;
	while(fmoves >> gsx >> gsy >> gex >> gey >> psx >> psy >> pex >> pey >> crl >> vlc >> err)
	{
		vec v(point(gsx, gsy), point(gex, gey));
		mvn.push_back(movement(v, vlc, err));
		// getline(fmoves, s);
	}
	fmoves.close();
}

/*wstring AnsiToWide(const string& in_sAnsi)
{
    wstring wsWide;
    wsWide.resize(in_sAnsi.length(), 0);
    MultiByteToWideChar(
        1251, 
        0, 
        &in_sAnsi[0], 
        (int)in_sAnsi.length(), 
        &wsWide[0], 
        (int)wsWide.length());

    return wsWide;
}*/

int main(int argc, char** argv)
{

	if (!analize_options(argc, argv))
	{
		cout << "Incorrect arguments!\n\n";
		print_eng_usage();
	}

	ofstream fres;	

	vector <movement> mvn;
	vector <scut> station;

	read_cuts(station_points_file, station);
	
	read_movement_field(move_points_file, mvn);


	// flog.open(output_log);
	// fitg.open(itg_log);

	fres.open(out_file);

	DynamicTopography dyn_tpg(mvn);

	for (size_t i = 0; i < station.size(); ++i)
	{

		dyn_tpg.set_cut(station[i]);
		dyn_tpg.set_file_index(i + 1);

		struct dt_result dt_res = dyn_tpg.take();

		fres << dt_res.toDTFormat() << endl;
		
	}

	fres.close();

	// flog.close();
	// fitg.close();

	return 0;
}