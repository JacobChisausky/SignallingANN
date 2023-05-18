/*
#ifndef PARAMETERS_H
#define PARAMETERS_H

#include "json.hpp"

struct parameters
{
	int	replicates	= 1;
	int	k	= 2;
	int	seed	= 123456789;
	double	N	= 1000.0;
	int	G	= 200000;
//	double	c	= 10.0;
//	double d = 5.0;
//	double p = 1.0;  // Use 0.5 for additive fitness

	int s_levels = 2;
	int q_levels = 2;

	double c1 = 0.5;
	double c0 = 1.5;

	double	init_ann_range	= 1.0;
	double	mut_rate_ann_S	= 0.01;
	double	mut_rate_ann_R	= 0.01;
	double	mut_step_ann_S	= 0.01;
	double	mut_step_ann_R	= 0.01;
//	bool	send_0_first	= true;
//	int	s_max	= 10;
	int	interactionPartners	= 10;
//	int fitnessFunction = 0;	//0 = additive. 1 = multiplicative
	bool	complexInit	= true;
//	bool	nullReceivers	= false;
//	bool	nullSenders	= false;
//	int	nullHonestBeginG	= 100000;
	int	Report_annVar	= 10;
	int	Report_annVar_N	= 100;
	bool	Report_annInit	= true;
	bool	recordFittestANNs = true;
	std::string	dataFileName	= "test_json";
	std::string	dataFileFolder	= "/.";
};

void from_json(const nlohmann::json& j, parameters& p);


#endif // SIM_PARAMETERS_H
*/
