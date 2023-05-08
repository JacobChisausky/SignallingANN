//============================================================================
// Name        : SignallingANN.cpp
// Author      : Jacob Chisausky
// Description : Signalling ANN
//============================================================================

#include <iostream>
#include <iostream>
#include <vector>
#include <algorithm>
#include <random>
#include <ctime>
#include <fstream>
#include <string>

#include "parameters.h"
#include "agents.h"
//using namespace std;

//To talk to GR about: should we introduce x new individuals each generation. I worry that mutations are too gradual - and we might want completely new individuals every once in a while

//These cost and benefit functions, and calculation of fitness, are what I really want to go over with Graeme

//Cost function
//double cost_function(double s, double q, double coeff){			//Add options to input too...
//For now - just linear
//	return s - q;
//}

//Benefit function for senders
double sender_fitness_function(bool response, double s, double q, double c, int interactionPartners){
	//This will be ran after each interaction and summed for all interaction. Each value is divided by interactionPartners to keep consistency across replicates - max fitness = 1, min = 0

	//This is multiplicative: fitness = benefits * survival probability
	//Where survival probability is 1 - costs

	//Benefits per interaction: 1 if receive response, 0 if not
	//Survival prob: 1 - s^1+cq
	//Follows Biernaskie et al 2014
	//c is a coefficient that determines the effect of q on signal costs.
	//Problem? If no responses are given, fitness is 0.

	double fitness = ( (int(response) * (1 - std::pow(s,(1+ (c*q))) ) ) / double(interactionPartners) );

	return fitness;

	//Write a version which ignores senders and just is grafen-esuqe with reprod success = s/sbar ?

}

//Benefit function for receivers
double receiver_benefit_function(bool response, double q){
	//If response sent: payoff = q
	//If response not sent: payoff = 1-q
	//So payoffs are always positive - but it is better to ignore q < 0.5 and respond to q > 0.5

	//Response sent			//Response not sent (response = 0)
	return (double(response)*q  +  -1*(double(response) - 1.0)*(1-q) );
}

//Benefit function for selecting null receivers where Pr = s
double receiver_benefit_function_PrEqualsS(double Pr, double s){
	return (1.0 - std::abs(Pr-s));
}

//json version
int main(int argc, char* argv[]){
	// getting params from the parameter file via json

	std::cout << argv[1] << std::endl;

	nlohmann::json json_in;
	std::ifstream is(argv[1]);   //assumes that the file name is given as a parameter in the command line
	is >> json_in;
	parameters params = json_in.get<parameters>();

	int replicates = params.replicates;
	int k = params.k;
	//constexpr int k = k1;
	int seed = params.seed;
	double N = params.N;
	int G = params.G;
	double c = params.c;
	double init_ann_range = params.init_ann_range;
	double mut_rate_ann_S = params.mut_rate_ann_S;
	double mut_rate_ann_R = params.mut_rate_ann_R;
	double mut_step_ann_S = params.mut_step_ann_S;
	double mut_step_ann_R = params.mut_step_ann_R;
	bool send_0_first = params.send_0_first;
	int s_max = params.s_max;
	int interactionPartners = params.interactionPartners;
	bool complexInit = params.complexInit;
	bool nullReceivers = params.nullReceivers;
	bool nullSenders = params.nullSenders;
	int nullHonestBeginG = params.nullHonestBeginG;
	int Report_annVar = params.Report_annVar;
	int Report_annVar_N = params.Report_annVar_N;
	bool Report_annInit = params.Report_annInit;
	std::string dataFileName = params.dataFileName;
	std::string dataFileFolder = params.dataFileFolder;

	//const std::string dataFileFolder = "./";

	/*
int main() {

	const int replicates = 1;

	const int k = 2;

	int seed = time(0);		//This won't work if we send many jobs at once but for now it's okay
	const double N = 1000;	//There will be N receivers and N senders. Stored as double for calculations
	const int G = 200000;

	const double c = 10.0;		//This determines the relationship between sender quality and signal cost. Higher = stronger reduction of cost with quality. 0 = same cost for all signallers.
	const double init_ann_range = 1.0;	//ANN stats will be initialized randomly + or - this value

	const double mut_rate_ann_S = 0.01;
	double mut_rate_ann_R = 0.01;

	const double mut_step_ann_S = 0.01;
	const double mut_step_ann_R = 0.05;	//Try .05 and 0.01

	const bool send_0_first = true; //if true, a signal of s=0.0 is always tried first - so that 'no signal' is an acceptable strategy. What does GR think of this? Not implemented yet.
	const int s_max = 10;			//The number of signals strengths tried before a signal of 0 is sent by default

	const int interactionPartners = 10;	//How many interactions per generation does each signaller and receiver engage in?

	const bool complexInit = true;		//If true, initial ANNs will not be allowed to be completely flat

	const bool nullReceivers = false;	//If true: receivers all respond with Pr = s. Their ANNs never mutate. This is for testing senders.
	const bool nullSenders = false;		//If true: senders only send s = q. Bypass ANNs. This is for testing receivers. Sender ANNs will drift.

	const int nullHonestBeginG = 210000;	//If >0, let everyone evolve against null senders and null receivers for this many generations. Then let things evolve. This should start it at an honest equilibrium.
	//This may not be working - in earlier tests it only took <10 000 gens to get receivers honest.
	//The way this has to work: Sender send normally but receiver behavior is null (Pr=s)
	//Receivers get cues from 'fake' senders where s=q.
	//After nullHonestBeginG generations, senders and receivers actually start interacting

	const int Report_annVar = 1; //Export this many generations of ANN variable data. 0 for no report, 1 for only last generation
	const int Report_annVar_N = 1000;	//How many individuals do you want to report ann stats for? The highest fitness individuals will be reported
	const bool Report_annInit = true;	//Do you want initial generation ann data?

	const std::string dataFileName = "make_honest_weights";

	const std::string dataFileFolder = "C:/Users/owner/Documents/S4/Simulation_ANN";
	//const std::string dataFileFolder = "./";
	 */


	//_________End parameter input

	//Get ANNs from .csv


	//Reports
	int ReportFreq_annVar = 0;
	if (Report_annVar > 0){
		ReportFreq_annVar = floor(G/Report_annVar);
	}

	//Time
	const std::time_t now = std::time(nullptr) ; // get the current time point
	const std::tm calendar_time = *std::localtime( std::addressof(now) ) ;

	tm *ltm = localtime(&now);
	int yday = ltm->tm_yday;

	std::string hr = std::to_string(calendar_time.tm_hour);
	if (hr.length() == 1){
		hr = "0" + hr;
	}

	std::string min = std::to_string(calendar_time.tm_min);
	if (min.length() == 1){
		min = "0" + min;
	}

	std::string sec = std::to_string(calendar_time.tm_sec);
	if (sec.length() == 1){
		sec = "0" + sec;
	}

	//For file names
	std::string strTime = std::to_string(yday) + "_" +  hr + "_" + min + "_" + sec;
	std::string str1 = dataFileFolder + "/" + strTime + "_annVars_" +  dataFileName + ".csv";
	std::string str2 = dataFileFolder + "/" + strTime + "_params_" + dataFileName + ".csv";
	//std::string str3 = dataFileFolder + "/" + strTime + "_summaryStats_" + dataFileName + ".csv";

	//Files to write to
	std::ofstream annVars;
	//	std::ofstream dataLog;
	std::ofstream paramFile;
	//	std::ofstream summaryStats;

	annVars.open(str1);
	//		dataLog.open(str1);
	paramFile.open(str2);
	//Prepare output files
	annVars << "rep,gen,indType,indNum,quality,fitness,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38";
	//		dataLog << "rep,gen,ind,indType,sendType,strategy,alphaBeta,fitness";
	paramFile << "replicates,k,seed,N,G,c,init_ann_range,mut_rate_ann_S,mut_rate_ann_R,mut_step_ann_S,mut_step_ann_R,s_max,send_0_first,interactionPartners,complexInit,nullReceivers,nullSenders,nullHonestBeginG,Report_annVar,Report_annVar_N,Report_annInit,dataFileName,dataFileFolder";
	paramFile << "\n"<<std::to_string(replicates)<<","<<std::to_string(k)<<","<<std::to_string(seed)<<","<<std::to_string(N)<<","<<std::to_string(G)<<","<<std::to_string(c)<<","<<std::to_string(init_ann_range)<<","<<std::to_string(mut_rate_ann_S)<<","<<std::to_string(mut_rate_ann_R)<<","<<std::to_string(mut_step_ann_S)<<","<<std::to_string(mut_step_ann_R)<<","<<std::to_string(s_max)<<","<<std::to_string(send_0_first)<<","<<std::to_string(interactionPartners)<<","<<std::to_string(complexInit)<<","<<std::to_string(nullReceivers)<<","<<std::to_string(nullSenders)<<","<<std::to_string(nullHonestBeginG)<<","<<std::to_string(Report_annVar)<<","<<std::to_string(Report_annVar_N)<<","<<std::to_string(Report_annInit)<<","<<dataFileName<<","<<dataFileFolder;

	//For k-selection tournament
	//std::array<int, k> tourn_arr;

	std::vector<int> tourn_arr;
	for (int i = 0; i < k; i++){
		tourn_arr.push_back(0);
	}

	//Disable mutations if using nullReceivers
	if (nullReceivers == true){
		mut_rate_ann_R = 0.000;
	}

	auto rng = std::default_random_engine {seed};
	std::uniform_real_distribution<double> prob(0,1);
	std::uniform_real_distribution<double> init_ann_dist(-std::abs(init_ann_range),std::abs(init_ann_range));
	std::uniform_int_distribution<int> randN(0,N-1);
	auto ann_mu_S = std::bernoulli_distribution(std::abs(mut_rate_ann_S)); //distribution for ann mutation chance - sender
	auto ann_mu_R = std::bernoulli_distribution(std::abs(mut_rate_ann_R)); //distribution for ann mutation chance - receivers

	auto ann_mu_size_S = std::cauchy_distribution<double>(0.0, std::abs(mut_step_ann_S));		//Size of ann mutations - senders
	auto ann_mu_size_R = std::cauchy_distribution<double>(0.0, std::abs(mut_step_ann_R));		//Size of ann mutations - receivers

	//Quality distribution
	//Only flat for now..
	std::uniform_real_distribution<double> qual_dist(0,1);
	//Just use prob(rng) to assign quality

	//Vector used for drawing random numbers from 0 to N-1 without replacement
	std::vector<int> nullVecS;
	for (int i = 0; i < N; i++){
		nullVecS.push_back(i);
	}

	//Vector used for drawing random numbers from 0 to N-1 without replacement
	std::vector<int> nullVecR;
	for (int i = 0; i < N; i++){
		nullVecR.push_back(i);
	}

	//Parameterized relationships - cost, benefit, etc...
	//1. Signal costs
	//2. Sender benefits (start flat)
	//3. quality distribution (start flat)
	//4. Receiver benefits?


	//Start replicate loop
	for (int rep = 1; rep <= replicates; rep++){

		//Population vectors
		std::vector<Sender> SenderPopulation;
		std::vector<Receiver> ReceiverPopulation;

		//Offspring vectors which will store offspring generation before it replaces parental generation
		std::vector<Sender> SenderOffspring;
		std::vector<Receiver> ReceiverOffspring;

		//Initialize Population

		for (int i = 0; i < N; i++){	//senders
			std::array<double, 39> sender_ann; //Create ann for sender
			if (complexInit == false){
				for (int j = 0; j < sender_ann.size(); j++){
					sender_ann[j] = init_ann_dist(rng);				//Create a random neural network for senders
				}
			} else if (complexInit==true){
				bool annOkay = false;
				while (annOkay == false){
					for (int j = 0; j < sender_ann.size(); j++){
						sender_ann[j] = init_ann_dist(rng);				//Create a random neural network for senders
					}
					//Now test the network and if it is complex enough, set annOkay to true
					annOkay = annS_test(10, sender_ann);
				}
			}

			//For now - flat quality distribution only
			double qualityInit = qual_dist(rng);

			SenderPopulation.push_back(Sender(sender_ann,qualityInit));		//Create sender agent  with that network
			SenderOffspring	.push_back(Sender(sender_ann,qualityInit));		//Offspring = parents in first generation.

		}

		if (nullReceivers == false){
			for (int i = 0; i < N; i++){	//receivers

				std::array<double, 34> receiver_ann; //Create ann for receiver
				if (complexInit == false){
					for (int j = 0; j < receiver_ann.size(); j++){		//Create a random neural network for receivers
						receiver_ann[j] = init_ann_dist(rng);
					}
				} else if (complexInit==true){
					bool annOkay = false;
					while (annOkay == false){
						for (int j = 0; j < receiver_ann.size(); j++){		//Create a random neural network for receivers
							receiver_ann[j] = init_ann_dist(rng);
						}
						//Now test the network and if it is complex enough, set annOkay to true
						annOkay = annR_test(10, receiver_ann);
					}
				}

				ReceiverPopulation.push_back(Receiver(receiver_ann));		//Create sender agent  with that network
				ReceiverOffspring.push_back(Receiver(receiver_ann));		//Offspring = parents in first generation.
			}
		} else if (nullReceivers == true){
			//Disable mutations

			//This ANN will give output = input	   //= {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20 , 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33};
			std::array<double, 34> receiver_ann_NULL = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  1,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0 ,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0};
			for (int i = 0; i < N; i++){

				ReceiverPopulation.push_back(Receiver(receiver_ann_NULL));		//Create sender agent  with that network
				ReceiverOffspring.push_back(Receiver(receiver_ann_NULL));		//Offspring = parents in first generation.
			}
		}

		//Now SenderPopulation and ReceiverPopulation are populated with agents of random networks and fitness 0.

		//Write initial ANNs to file - unless 0 reports are requested
		if (Report_annVar > 0 & Report_annInit == true){
			for (int i = 0; i < Report_annVar_N; i++){
				annVars<<"\n"<<rep<<","<<0<<",Sender,"<<i<<","<<SenderPopulation[i].get_quality()<<","<<0;
				for (int m = 0; m <= 38; m++){
					annVars<<","<<SenderPopulation[i].get_ann(m);
				}
			}
			for (int i = 0; i < Report_annVar_N; i++){
				annVars<<"\n"<<rep<<","<<0<<",Receiver,"<<i<<","<<0<<","<<0;
				for (int m = 0; m <= 38; m++){
					if (m <= 33){
						annVars<<","<<ReceiverPopulation[i].get_ann(m);
					} else {
						annVars<<","<<0;
					}
				}
			}
		}


		//Start Generation Loop
		for (int g = 0; g <= G; g++){

			//Assign qualities to senders in offspring loop, not here.

			//Randomize order of sender and receiver populations
			//Then pair sender 1 with receiver 1, 2 with 2, etc
			//Randomize order again and repeat interactionPartner # of times

			for (int i = 0; i < interactionPartners; i++){

				//These vectors are shuffled to randomly make S-R pairs
				std::shuffle(std::begin(nullVecS), std::end(nullVecS), rng);
				std::shuffle(std::begin(nullVecR), std::end(nullVecR), rng);

				for (int j = 0; j < N; j++){

					Sender S_cur = SenderPopulation[nullVecS[j]];
					double q_cur = S_cur.get_quality();
					Receiver R_cur = ReceiverPopulation[nullVecR[j]];
					double s = 0;

					//Determine strength of signal sent by signaller
					//Share for all individuals per generation? No - this could produce a generation with mostly strong signals, another with mostly weak, etc.
					//vector<double> s_list();
					if (nullSenders == true){	//Null sender behavior ONLY FOR RECEIVERS xxx
						//They are not actually coevolving yet
						s = q_cur;
					} else {
						//First - try s = 0
						if (send_0_first == true & prob(rng) > S_cur.annS_output(0.0, q_cur)){    //If so, then signal of 0 is NOT sent. Go on to pick more signals. If not (prob < ann output), then keep signal = 0;
							int s_num = 0;
							while (s_num < s_max){
								s_num++;
								double s_try = prob(rng);
								if (prob(rng) < S_cur.annS_output(s_try, q_cur)){
									s = s_try;
									s_num = s_max;
								}
							}	//If we don't pick a signal after s_max tries, keep at s = 0
						} else { //don't send 0 first
							int s_num = 0;
							while (s_num < s_max){
								s_num++;
								double s_try = prob(rng);
								if (prob(rng) < S_cur.annS_output(s_try, q_cur)){
									s = s_try;
									s_num = s_max;
								}
							}	//If we don't pick a signal after s_max tries, keep at s = 0
						}
					}

					//Determine cost of signal
					//double cost_S = cost_function(s, q_cur, cost_function_coeff);

					//Fitness = benefit - cost...? additive not multiplicative - is that okay? Talk to GR about it.

					//Null receiver behavior. Receiver gets payoff from response. Sender gets payoff from nullResponse
					double payoffReceiver = -1.0;
					bool response = 0;
					bool nullResponse = 0;
					if (g < nullHonestBeginG){
						double Pr = R_cur.annR_output(s);	//Use different benefit function for now to make Pr = s for receivers
						payoffReceiver = receiver_benefit_function_PrEqualsS(Pr, s);

						//old
						//if (prob(rng) < s){    //Does receiver respond to signal? Chance is equal to s for NULL RECEIVERS
						//	nullResponse = 1;		//null response is response for SENDER
						//}

						//double sig = std::pow((1/(2+c*q_cur)),(1/(1+c*q_cur)));
						//old
						//if (prob(rng) < R_cur.annR_output(s)){
						//	response = 1;			//response is response for receiver fitness
						//}
					} else {
						//Now check receiver ANN - do they respond to signal?
						//For real evolution - not honest start
						if (prob(rng) < R_cur.annR_output(s)){    //If so - respond to signal
							response = 1;
						}
						payoffReceiver = receiver_benefit_function(response, q_cur);
					}


					//Benefit to sender
					//double benefit_S = sender_benefit_function(response, q_cur, sender_benefit_option);

					//Payoffs to senders and receivers
					double payoffSender = sender_fitness_function(response, s, q_cur, c, interactionPartners);

					if (g < nullHonestBeginG){
						//payoffReceiver = receiver_benefit_function(response, q_cur); //I moved this up
						payoffSender = sender_fitness_function(nullResponse, s, q_cur, c, interactionPartners);
					}

					//Fitness changes
					SenderPopulation[nullVecS[j]].change_fitness(payoffSender);
					ReceiverPopulation[nullVecR[j]].change_fitness(payoffReceiver);

				}//End loop for this individual

			}//End loop for interaction partners
			//Now all fitnesses have been determined.

			//k-tournament selection
			//First, pick k individuals (k=2)
			//Pick *with replacement*
			//Determine highest fitness one
			//If tie, pick random
			//Winner reproduces
			//Repeat
			//randN(rng) produces random int from 0 to N-1

			for (int n = 0; n < N; n++){	//Reproduction, mutation

				tourn_arr[0] = randN(rng);
				double maxFit = SenderPopulation[tourn_arr[0]].get_fitness();
				int maxFitInd = 0;

				for (int k_cur = 1; k_cur < k; k_cur++){
					tourn_arr[k_cur] = randN(rng);

					while (tourn_arr[1] == tourn_arr[0]){ //This will only remove duplicates for k = 2 but since that is what we plan on using, this is fine
						tourn_arr[1] = randN(rng);
					}

					if (SenderPopulation[tourn_arr[k_cur]].get_fitness() > maxFit){
						maxFit = SenderPopulation[tourn_arr[k_cur]].get_fitness();
						maxFitInd = k_cur;
					}
				}

				SenderOffspring[n] = SenderPopulation[tourn_arr[maxFitInd]];
				SenderOffspring[n].reset_fitness();
				SenderOffspring[n].set_quality(qual_dist(rng));

				//Mutation for SenderOffspring[n]
				for (int m = 0; m <= 38; m++){
					if (ann_mu_S(rng)){	//Mutation to this weight occurs
						//which ann variable to mutate, size of mutation
						SenderOffspring[n].ann_mutate(m, ann_mu_size_S(rng));
					}
				}

				//A new sender is born!

				//Now the same for receivers
				tourn_arr[0] = randN(rng);
				maxFit = ReceiverPopulation[tourn_arr[0]].get_fitness();
				maxFitInd = 0;

				for (int k_cur = 1; k_cur < k; k_cur++){
					tourn_arr[k_cur] = randN(rng);

					while (tourn_arr[1] == tourn_arr[0]){ //This will only remove duplicates for k = 2 but since that is what we plan on using, this is fine
						tourn_arr[1] = randN(rng);
					}

					if (ReceiverPopulation[tourn_arr[k_cur]].get_fitness() > maxFit){
						maxFit = ReceiverPopulation[tourn_arr[k_cur]].get_fitness();
						maxFitInd = k_cur;
					}
				}

				ReceiverOffspring[n] = ReceiverPopulation[tourn_arr[maxFitInd]];
				ReceiverOffspring[n].reset_fitness();

				//Mutation for SenderOffspring[n]
				for (int m = 0; m <= 33; m++){
					if (ann_mu_R(rng)){	//Mutation to this weight occurs
						//which ann variable to mutate, size of mutation
						ReceiverOffspring[n].ann_mutate(m, ann_mu_size_R(rng));
					}
				}
			}

			if (Report_annVar > 0 & g > 0 & g%ReportFreq_annVar == 0){	//Report ANN stats for fittest individuals

				std::vector<int> fittestSenders;	//Find fittest individuals in population by populating this vector with their indices
				fittestSenders.push_back(0);
				double minFit = SenderPopulation[0].get_fitness();
				int minFitNum = 0;	//position in vector

				//Now populate with Report_annVar_N individuals, no replacement
				for (int i = 1; i < Report_annVar_N; i++){
					fittestSenders.push_back(i);
					if (SenderPopulation[i].get_fitness() < minFit){
						minFit = SenderPopulation[i].get_fitness();
						minFitNum = i;
					}
				}

				for (int i = Report_annVar_N; i < N; i++){
					if (SenderPopulation[i].get_fitness() > minFit){
						fittestSenders[minFitNum] = i;
						//Find new least fit individual
						minFit = SenderPopulation[fittestSenders[0]].get_fitness();
						minFitNum = 0;
						for (int j = 1; j < Report_annVar_N; j++){
							if (SenderPopulation[fittestSenders[j]].get_fitness() < minFit){
								minFit = SenderPopulation[fittestSenders[j]].get_fitness();
								minFitNum = j;
							}
						}
					}
				}
				//By the end of this, fittestSenders will contain the indices of the Report_annVar_N fittest senders in the population
				//Do the same for receivers
				std::vector<int> fittestReceivers;
				fittestReceivers.push_back(0);
				minFit = ReceiverPopulation[0].get_fitness();
				minFitNum = 0;	//position in vector

				//Now populate with Report_annVar_N individuals, no replacement
				for (int i = 1; i < Report_annVar_N; i++){
					fittestReceivers.push_back(i);
					if (ReceiverPopulation[i].get_fitness() < minFit){
						minFit = ReceiverPopulation[i].get_fitness();
						minFitNum = i;
					}
				}

				for (int i = Report_annVar_N; i < N; i++){
					if (ReceiverPopulation[i].get_fitness() > minFit){
						fittestReceivers[minFitNum] = i;  //replace old least fit individual
						//Find new least fit individual
						minFit = ReceiverPopulation[fittestReceivers[0]].get_fitness();
						minFitNum = 0;
						for (int j = 1; j < Report_annVar_N; j++){
							if (ReceiverPopulation[fittestReceivers[j]].get_fitness() < minFit){
								minFit = ReceiverPopulation[fittestReceivers[j]].get_fitness();
								minFitNum = j;
							}
						}
					}
				}


				for (int i = 0; i < Report_annVar_N; i++){ // fittest individuals as found above
					annVars<<"\n"<<rep<<","<<g<<",Sender,"<<fittestSenders[i]<<","<<SenderPopulation[fittestSenders[i]].get_quality()<<","<<SenderPopulation[fittestSenders[i]].get_fitness();
					for (int m = 0; m <= 38; m++){
						annVars<<","<<SenderPopulation[fittestSenders[i]].get_ann(m);
					}
				}
				for (int i = 0; i < Report_annVar_N; i++){
					annVars<<"\n"<<rep<<","<<g<<",Receiver,"<<fittestReceivers[i]<<","<<0<<","<<ReceiverPopulation[fittestReceivers[i]].get_fitness();
					for (int m = 0; m <= 38; m++){
						if (m <= 33){
							annVars<<","<<ReceiverPopulation[fittestReceivers[i]].get_ann(m);
						} else {
							annVars<<","<<0;
						}
					}
				}
			}


			//Replace Parents with Offspring
			SenderPopulation.swap(SenderOffspring);
			ReceiverPopulation.swap(ReceiverOffspring);


		}//End Generation Loop

		//For testing...
		std::cout << "\n-----------------\n";
		for (int i = 0; i < 10; i++){
			for (int j = 0; j <=33; j++){
				std::cout << ReceiverPopulation[i].get_ann(j) << " ";
			}
			std::cout << "\n";
		}

	}//End replicate loop

	//	dataLog.close();
	paramFile.close();
	annVars.close();
	//	summaryStats.close();





	std::cout << "\nDone";

	return 0;
}
