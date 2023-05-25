
#ifndef PARAMETERS_H
#define PARAMETERS_H

#include "json.hpp"

struct parameters
{
	int initOption = 1;
	int	maxTrainingTime	=	100000;
	double	targetAccuracy	=	0.98;
	int	replicates	=	1;
	int	s_levels	=	2;
	int	q_levels	=	2;
	int	r_levels	=	2;
	double	cMax	=	2.0;
	double	cMin	=	0.0;
	double	m	=	0.50;
	int	seed	=	12345678;
	double	N	=	1000;
	int	G	=	2500;
	double	mut_rate_ann_S	=	0.01;
	double	mut_rate_ann_R	=	0.01;
	double	mut_step_ann_S	=	0.01;
	double	mut_step_ann_R	=	0.01;
	int	tries_max	=	10;
	int	interactionPartners	=	10;
	int	k	=	2;
	double	init_ann_range	=	1;
	bool	complexInit	=	1;
	int	Report_annVar	=	5;
	int	Report_annVar_N	=	100;
	bool	Report_annInit	=	1;
	bool	recordFittestANNs	=	1;
	bool	nullReceivers	=	false;
	bool	nullSenders	=	false;
	std::string	dataFileName	=	"test_correct_c_noNull";
	std::string	dataFileFolder	=	"C:/Users/owner/eclipse-workspace/SignallingANN/data";
};

void from_json(const nlohmann::json& j, parameters& p);


#endif // SIM_PARAMETERS_H

