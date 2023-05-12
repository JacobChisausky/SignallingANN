#include "parameters.h"

void from_json(const nlohmann::json& j, parameters& t)
{
	NLOHMANN_JSON_FROM(replicates);
	NLOHMANN_JSON_FROM(k);
	NLOHMANN_JSON_FROM(seed);
	NLOHMANN_JSON_FROM(N);
	NLOHMANN_JSON_FROM(G);
	NLOHMANN_JSON_FROM(c);
	NLOHMANN_JSON_FROM(p);
	NLOHMANN_JSON_FROM(init_ann_range);
	NLOHMANN_JSON_FROM(mut_rate_ann_S);
	NLOHMANN_JSON_FROM(mut_rate_ann_R);
	NLOHMANN_JSON_FROM(mut_step_ann_S);
	NLOHMANN_JSON_FROM(mut_step_ann_R);
	NLOHMANN_JSON_FROM(send_0_first);
	NLOHMANN_JSON_FROM(s_max);
	NLOHMANN_JSON_FROM(interactionPartners);
	NLOHMANN_JSON_FROM(fitnessFunction);
	NLOHMANN_JSON_FROM(complexInit);
	NLOHMANN_JSON_FROM(nullReceivers);
	NLOHMANN_JSON_FROM(nullSenders);
	NLOHMANN_JSON_FROM(nullHonestBeginG);
	NLOHMANN_JSON_FROM(Report_annVar);
	NLOHMANN_JSON_FROM(Report_annVar_N);
	NLOHMANN_JSON_FROM(Report_annInit);
	NLOHMANN_JSON_FROM(dataFileName);
	NLOHMANN_JSON_FROM(dataFileFolder);
}
