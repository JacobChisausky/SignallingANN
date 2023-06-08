
#include "parameters.h"

void from_json(const nlohmann::json& j, parameters& t)
{
	NLOHMANN_JSON_FROM(initOption);
	NLOHMANN_JSON_FROM(maxTrainingTime);
	NLOHMANN_JSON_FROM(targetAccuracy);
	NLOHMANN_JSON_FROM(replicates);
	NLOHMANN_JSON_FROM(s_levels);
	NLOHMANN_JSON_FROM(q_levels);
	NLOHMANN_JSON_FROM(r_levels);
	NLOHMANN_JSON_FROM(cMax);
	NLOHMANN_JSON_FROM(cMin);
	NLOHMANN_JSON_FROM(m);
	NLOHMANN_JSON_FROM(seed);
	NLOHMANN_JSON_FROM(N);
	NLOHMANN_JSON_FROM(G);
	NLOHMANN_JSON_FROM(p_q);
	NLOHMANN_JSON_FROM(p_s);
	NLOHMANN_JSON_FROM(p_rS);
	NLOHMANN_JSON_FROM(p_rR);
	NLOHMANN_JSON_FROM(p_m);
	NLOHMANN_JSON_FROM(try_0_first_S);
	NLOHMANN_JSON_FROM(try_0_first_R);
	NLOHMANN_JSON_FROM(mut_rate_ann_S);
	NLOHMANN_JSON_FROM(mut_rate_ann_R);
	NLOHMANN_JSON_FROM(mut_step_ann_S);
	NLOHMANN_JSON_FROM(mut_step_ann_R);
	NLOHMANN_JSON_FROM(tries_max);
	NLOHMANN_JSON_FROM(interactionPartners);
	NLOHMANN_JSON_FROM(k);
	NLOHMANN_JSON_FROM(init_ann_range);
	NLOHMANN_JSON_FROM(complexInit);
	NLOHMANN_JSON_FROM(Report_annVar);
	NLOHMANN_JSON_FROM(Report_annVar_N);
	NLOHMANN_JSON_FROM(Report_annInit);
	NLOHMANN_JSON_FROM(recordFittestANNs);
	NLOHMANN_JSON_FROM(nullReceivers);
	NLOHMANN_JSON_FROM(nullSenders);
	NLOHMANN_JSON_FROM(dataFileName);
	NLOHMANN_JSON_FROM(dataFileFolder);
}
