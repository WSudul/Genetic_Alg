#pragma once
#include "Selection.h"
#include "array"


/*!
	Small class for storing reservoir parameters
	intParams:
		-h		average thickness
		-k		average permeability of reservoir rock
		-r_res	average radius of reservoir
		-r_aqu	average radius of aquifier (r_aqu>r_res)

	floatParams:
		-fi		average porosity
		-u_w	water viscosity in reservoir
		-c_w	water compressibility coefficient (1/pressure unit)
		-c_f	formation compressibility coefficient (1/pressure unit)
		-f		fraction of aquifier encirclement
		-s_w	water saturation
*/
class ReservoirParams
{
public:
	ReservoirParams();
	~ReservoirParams();

	//int param1;
	//int param2;

	std::array<int, 4>	intParams;
	std::array<float, 6>	floatParams;


};

