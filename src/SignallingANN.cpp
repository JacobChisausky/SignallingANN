//============================================================================
// Name        : SignallingANN.cpp
// Author      : Jacob Chisausky
// Description : Signalling ANN - Discrete
//============================================================================

//Some open questions:
//How sensitive are (any) results to cost functions?
//Difficult to test..
//What if responding is competetive? I.e., receivers can only respond to one sender
//I expect that this would have very different dynamics.

#include <iostream>
#include <iostream>
#include <vector>
#include <algorithm>
#include <random>
#include <ctime>
#include <fstream>
#include <string>

//Temporary for non-json testing
//#include "parameters.h"
#include "agents.h"


//Cost function - since discrete model
double cost_function(double s, double q, double c0, double c1){
	//c1 = the cost of a s = 1 signal to q = 1
	//c0 = the cost of a s = 1 signal to q = 0
	double cost = (1.0 - q)*(s*c0 - s*c1) + (s*c1);
	return cost;
}

double sender_benefit_function(double r){
	//Can easily implement differential benefits here (b1 and b2)
	//For r_levels = 2: this is fine.
	//As we increase r_levels, we may need to tweak this.
	//Otherwise, sender benefits will be low relative to costs.
	//Adding a parameter to increase benefits (like 2.0*r) might be the easiest solution.

	return r;
}

//Receiver benefit function
//r has levels.
//It is best to respond with high r to high q
//It is best to respond with low r to low q
//Or - is r = q best?
//In this case with 2 r levels:
//r = 1 to q = 1 is good (A1 to T1 in Zollman et al terms)
//r =0 to q = 0 is good (A2 to T2)
//Fitness = 1-std(r-q)

double receiver_benefit_function(double r, double q){
	double benefit = 1.0 - std::abs(r-q);
	return benefit;
}


int main(){

	//To add:
	int initOption = 1;
	int maxTrainingTime = 100000;
	double targetAccuracy = 0.9;

	int replicates = 1;

	int s_levels = 2;		//Number of s values to use
	int q_levels = 2;		//Number of q values to use
	int r_levels = 2;		//Number of r values to use. For now only 2. Talk to GR about this.

	double c0 = 1.5;
	double c1 = 0.5;

	int seed = 12345678;
	double N = 1000;
	int G = 5000;

	double mut_rate_ann_S = 0.01;
	double mut_rate_ann_R = 0.01;
	double mut_step_ann_S = 0.01;
	double mut_step_ann_R = 0.01;

	//bool send_0_first = params.send_0_first;
	//int s_max = params.s_max;

	int interactionPartners = 10;
	int k = 2;

	double init_ann_range = 1;
	bool complexInit = 1;

	int Report_annVar = 10;
	int Report_annVar_N = 100;
	bool Report_annInit = 1;
	bool recordFittestANNs = 1;

	std::string dataFileName = "testCPP";
	std::string dataFileFolder = "C:/Users/owner/eclipse-workspace/SignallingANN/data";


	/*
//json version
int main(int argc, char* argv[]){
	// getting params from the parameter file via json

	std::cout << argv[1] << std::endl;

	nlohmann::json json_in;
	std::ifstream is(argv[1]);   //assumes that the file name is given as a parameter in the command line
	is >> json_in;
	parameters params = json_in.get<parameters>();

	int replicates = params.replicates;

	int s_levels = params.s_levels;		//Number of s values to use
	int q_levels = params.q_levels;		//Number of q values to use
	//int r_levels = params.r_levels;		//Number of r values to use. For now only 2. Talk to GR about this.

	double c0 = params.c0;
	double c1 = params.c1;

	//constexpr int k = k1;
	int seed = params.seed;
	double N = params.N;
	int G = params.G;

	//double c = params.c;
	//double d = params.d;
	//double p = params.p;
//	int fitnessFunction = params.fitnessFunction; //0 = additive. 1 = multiplicative

	double mut_rate_ann_S = params.mut_rate_ann_S;
	double mut_rate_ann_R = params.mut_rate_ann_R;
	double mut_step_ann_S = params.mut_step_ann_S;
	double mut_step_ann_R = params.mut_step_ann_R;

	//bool send_0_first = params.send_0_first;
	//int s_max = params.s_max;

	int interactionPartners = params.interactionPartners;
	int k = params.k;

	double init_ann_range = params.init_ann_range;
	bool complexInit = params.complexInit;

	//Replace these with 'start' options to run before the real sim begins instead of current initialization
	//bool nullReceivers = params.nullReceivers;
	//bool nullSenders = params.nullSenders;
	//int nullHonestBeginG = params.nullHonestBeginG;

	int Report_annVar = params.Report_annVar;
	int Report_annVar_N = params.Report_annVar_N;
	bool Report_annInit = params.Report_annInit;
	bool recordFittestANNs = params.recordFittestANNs;

	std::string dataFileName = params.dataFileName;
	std::string dataFileFolder = params.dataFileFolder;
	 */
	//_________End parameter input

	//Get ANNs from .csv eventually for honest start?


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
	paramFile << "replicates,s_levels,q_levels,c0,c1,seed,N,G,mut_rate_ann_S,mut_rate_ann_R,mut_step_ann_S,mut_step_ann_R,interactionPartners,k,init_ann_range,complexInit,Report_annVar,Report_annVar_N,Report_annInit,recordFittestANNs,dataFileName,dataFileFolder";
	paramFile << "\n"<< std::to_string(replicates) << "," << std::to_string(s_levels) << "," << std::to_string(q_levels) << "," << std::to_string(c0)<<","<<std::to_string(c1)<<","<<std::to_string(seed)<<","<<std::to_string(N)<<","<<std::to_string(G)<<","<<std::to_string(mut_rate_ann_S)<<","<<std::to_string(mut_rate_ann_R)<<","<<std::to_string(mut_step_ann_S)<<","<<std::to_string(mut_step_ann_R)<<","<<std::to_string(interactionPartners)<<","<<std::to_string(k)<<","<<std::to_string(init_ann_range)<<","<<std::to_string(complexInit)<<","<<std::to_string(Report_annVar)<<","<<std::to_string(Report_annVar_N)<<","<<std::to_string(Report_annInit)<<","<<std::to_string(recordFittestANNs)<<","<<dataFileName<<","<<dataFileFolder;

	//Determine what values of s and q are used
	//These are evenly spaced from 0 to 1
	if (s_levels < 2 || q_levels < 2 || r_levels < 2){
		std::cout << "s_levels and q_levels and r_levels must be > 1";
		return 97;
	}


	std::vector<double> s_vals;
	double s_increment = 1.0/double(s_levels-1);
	for (double val = 0.0; val <= 1.0; val += s_increment){
		s_vals.push_back(val);
	}
	if (s_vals.size() != s_levels){
		return 88;
	}

	std::vector<double> q_vals;
	double q_increment = 1.0/double(q_levels-1);
	for (double val = 0.0; val <= 1.0; val += q_increment){
		q_vals.push_back(val);
	}
	if (q_vals.size() != q_levels){
		return 88;
	}

	std::vector<double> r_vals;
	double r_increment = 1.0/double(r_levels-1);
	for (double val = 0.0; val <= 1.0; val += r_increment){
		r_vals.push_back(val);
	}
	if (r_vals.size() != r_levels){
		return 88;
	}


	//For k-selection tournament
	std::vector<int> tourn_arr;
	for (int i = 0; i < k; i++){
		tourn_arr.push_back(0);
	}

	//Disable mutations if using nullReceivers
	//if (nullReceivers == true){
	//	mut_rate_ann_R = 0.000;
	//}

	auto rng = std::default_random_engine {seed};
	std::uniform_real_distribution<double> prob(0,1);
	std::uniform_real_distribution<double> init_ann_dist(-std::abs(init_ann_range),std::abs(init_ann_range));
	std::uniform_int_distribution<int> randN(0,N-1);
	std::uniform_int_distribution<int> randQ(0,(q_levels-1));
	std::uniform_int_distribution<int> randS(0,(s_levels-1));
	std::uniform_int_distribution<int> randR(0,(r_levels-1));
	auto ann_mu_S = std::bernoulli_distribution(std::abs(mut_rate_ann_S)); //distribution for ann mutation chance - sender
	auto ann_mu_R = std::bernoulli_distribution(std::abs(mut_rate_ann_R)); //distribution for ann mutation chance - receivers

	auto ann_mu_size_S = std::cauchy_distribution<double>(0.0, std::abs(mut_step_ann_S));		//Size of ann mutations - senders
	auto ann_mu_size_R = std::cauchy_distribution<double>(0.0, std::abs(mut_step_ann_R));		//Size of ann mutations - receivers

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

	//Start replicate loop
	for (int rep = 1; rep <= replicates; rep++){

		//Population vectors
		std::vector<Sender> SenderPopulation;
		std::vector<Receiver> ReceiverPopulation;

		//Offspring vectors which will store offspring generation before it replaces parental generation
		std::vector<Sender> SenderOffspring;
		std::vector<Receiver> ReceiverOffspring;

		//Initialize Population
		//Begin with random ANNs

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

			//Flat quality distribution - discrete

			double qualityInit = q_vals[randQ(rng)];

			SenderPopulation.push_back(Sender(sender_ann,qualityInit));		//Create sender agent  with that network
			SenderOffspring.push_back(Sender(sender_ann,qualityInit));		//Offspring = parents in first generation.
		}

		for (int i = 0; i < N; i++){	//receivers
			std::array<double, 39> receiver_ann; //Create ann for receiver
			if (complexInit == false){
				for (int j = 0; j < receiver_ann.size(); j++){
					receiver_ann[j] = init_ann_dist(rng);				//Create a random neural network for senders
				}
			} else if (complexInit==true){
				bool annOkay = false;
				while (annOkay == false){
					for (int j = 0; j < receiver_ann.size(); j++){
						receiver_ann[j] = init_ann_dist(rng);				//Create a random neural network for senders
					}
					//Now test the network and if it is complex enough, set annOkay to true
					annOkay = annR_test(10, receiver_ann);
				}
			}
			ReceiverPopulation.push_back(Receiver(receiver_ann));		//Create sender agent  with that network
			ReceiverOffspring.push_back(Receiver(receiver_ann));		//Offspring = parents in first generation.
		}

		//Old receiver ANNs
		//if (nullReceivers == false){
		/*	for (int i = 0; i < N; i++){	//receivers

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
		}*/
		/*} else if (nullReceivers == true){
			//Disable mutations

			//This ANN will give output = input	   //= {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20 , 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33};
			std::array<double, 34> receiver_ann_NULL = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0,  1,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0 ,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0};
			for (int i = 0; i < N; i++){

				ReceiverPopulation.push_back(Receiver(receiver_ann_NULL));		//Create sender agent  with that network
				ReceiverOffspring.push_back(Receiver(receiver_ann_NULL));		//Offspring = parents in first generation.
			}
		}*/

		//Now SenderPopulation and ReceiverPopulation are populated with agents of random networks and fitness 0.

		//Initialize option: Honest - take the ANNs and train them to prespecified behavior.

		int genAcheived = maxTrainingTime;
		if (initOption == 1){	//This is for s_levels = 2, q_levels = 2. This is a Zollman Honest Equilibrium
			//What if we init at hybrid equilibrium - q = 0, prob sending =

			//Signaller: s = q.
			//A1 with Pr = q

			double maxDeviation = N;
			double targetDeviation = (1.0-targetAccuracy)*maxDeviation;
			bool targetAcheived = false;

			//Train signallers
			for (int trainingGen = 1; trainingGen <= maxTrainingTime; trainingGen++){
				if (targetAcheived == true){
					break;
				}

				//Total deviation from perfect performance for receivers and senders
				double totalDeviationR = 0.0;
				double totalDeviationS = 0.0;

				for (int i = 0; i < N; i++){

					//Senders
					Sender S_cur = SenderPopulation[i];
					double q_rand = q_vals[randQ(rng)]; //random q

					//Try a random q. If s = q, reward.


					double s = 0;
					bool s_selected = false;
					for (int tries = 0; tries < 10; tries++){
						//Shuffle s_vals
						std::shuffle(std::begin(s_vals), std::end(s_vals), rng);
						//Try each value in vector s_vals. This is sampling without replacement.
						for (int s_try = 0; s_try < s_levels; s_try++){

							if (prob(rng) < S_cur.annS_output(s_vals[s_try], q_rand)){
								s = s_vals[s_try];
								s_selected = true;
								break;
							}
						}
						if (s_selected == true){
							break;
						}
					}

					//s is now selected.
					//Fitness = 1 - abs(q_rand - s)
					double fitnessS = 1.0 - std::abs(q_rand - s);
					SenderPopulation[i].change_fitness(fitnessS);

					//Receivers

					Receiver R_cur = ReceiverPopulation[i];
					double s_rand = s_vals[randS(rng)]; //random s

					//Try a random s. Output r = s is best

					double r = 0;
					bool r_selected = false;
					for (int tries = 0; tries < 10; tries++){
						//Shuffle r_vals
						std::shuffle(std::begin(r_vals), std::end(r_vals), rng);
						//Try each value in vector s_vals. This is sampling without replacement.
						for (int r_try = 0; r_try < r_levels; r_try++){

							if (prob(rng) < R_cur.annR_output(s_rand, r_vals[r_try])){
								r = r_vals[r_try];
								r_selected = true;
								break;
							}
						}
						if (r_selected == true){
							break;
						}
					}

					double fitnessR = receiver_benefit_function(r, s_rand);	//Giving s and NOT q to the benefit function. At honest state, s = q.
					ReceiverPopulation[i].change_fitness(fitnessR);

					//Total deviation for receivers and senders
					totalDeviationR += (1.0 - fitnessR);
					totalDeviationS += (1.0 - fitnessS);

					//Selection: k-selection tournament. Set k to one nonparameter value after testing
					for (int n = 0; n < N; n++){	//Reproduction, mutation

						tourn_arr[0] = randN(rng);
						double maxFit = SenderPopulation[tourn_arr[0]].get_fitness();
						int maxFitInd = 0;

						for (int k_cur = 1; k_cur < k; k_cur++){
							tourn_arr[k_cur] = randN(rng);

							while (tourn_arr[1] == tourn_arr[0]){ //This will only remove duplicates for k = 2
								tourn_arr[1] = randN(rng);
							}

							if (SenderPopulation[tourn_arr[k_cur]].get_fitness() > maxFit){
								maxFit = SenderPopulation[tourn_arr[k_cur]].get_fitness();
								maxFitInd = k_cur;
							}
						}

						SenderOffspring[n] = SenderPopulation[tourn_arr[maxFitInd]];
						SenderOffspring[n].reset_fitness();
						//		SenderOffspring[n].set_quality(q_vals[randQ(rng)]); quality doesn't matter since it is randomly assigned above

						//Mutation for SenderOffspring[n]
						for (int m = 0; m <= 38; m++){
							if (ann_mu_S(rng)){	//Mutation to this weight occurs
								//which ann variable to mutate, size of mutation
								SenderOffspring[n].ann_mutate(m, ann_mu_size_S(rng));
							}
						}

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
						for (int m = 0; m <= 38; m++){
							if (ann_mu_R(rng)){	//Mutation to this weight occurs
								//which ann variable to mutate, size of mutation
								ReceiverOffspring[n].ann_mutate(m, ann_mu_size_R(rng));
							}
						}
					}

				}//end indvidual loop
				//Check total deviation - break loop if it is low enough
				//Max deviation = 1 * N * N
				//Total deviation holds how 'wrong' the networks are. N*N = worst, 0 = best
				//We want to acheive less than targetAccuracy % deviation.
				if (totalDeviationR < targetDeviation && totalDeviationS < targetDeviation){
					genAcheived = trainingGen;
					targetAcheived = true;
				}
				//Replace Parents with Offspring
				SenderPopulation.swap(SenderOffspring);
				ReceiverPopulation.swap(ReceiverOffspring);

				//std::cout << "r " << totalDeviationR << "\ns " << totalDeviationS <<"\n";
			}//end trainingGen loop
		}
		//	return 2;

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
					annVars<<","<<ReceiverPopulation[i].get_ann(m);
				}
			}
		}

		//Start Generation Loop
		for (int g = 1; g <= G; g++){

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
					//		if (nullSenders == true){	//Null sender behavior ONLY FOR RECEIVERS
					//They are not actually coevolving yet
					//This probably needs to be scrapped
					//			s = q_cur;
					//		} else {

					bool s_selected = false;
					for (int tries = 0; tries++; tries < 10){
						//Shuffle s_vals
						std::shuffle(std::begin(s_vals), std::end(s_vals), rng);
						//Try each value in vector s_vals. This is sampling without replacement.
						for (int s_try = 0; s_try++; s_try < s_levels){
							if (s_vals[s_try] < S_cur.annS_output(s_vals[s_try], q_cur)){
								s = s_vals[s_try];
								s_selected = true;
								break;
							}
						}
						if (s_selected == true){
							break;
						}
						//If we select one, break
						//If not, repeat. For now, try 10 times
						//If still nothing, send 0
					}

					//Receivers

					double r = 0;
					bool r_selected = false;
					for (int tries = 0; tries < 10; tries++){
						//Shuffle r_vals
						std::shuffle(std::begin(r_vals), std::end(r_vals), rng);
						//Try each value in vector s_vals. This is sampling without replacement.
						for (int r_try = 0; r_try < r_levels; r_try++){
							if (prob(rng) < R_cur.annR_output(s, r_vals[r_try])){
								r = r_vals[r_try];
								r_selected = true;
								break;
							}
						}
						if (r_selected == true){
							break;
						}
					}

					//Modify sender fitness
					double senderCost = cost_function(s, q_cur, c0, c1);
					double senderBenefit = sender_benefit_function(r);
					double senderPayoff = senderBenefit - senderCost;
					SenderPopulation[nullVecS[j]].change_fitness(senderPayoff);

					// Modify receiver fitness
					ReceiverPopulation[nullVecR[j]].change_fitness(receiver_benefit_function(r, q_cur));

					//Old - before discrete implementation
					/*	if (send_0_first == true & prob(rng) > S_cur.annS_output(0.0, q_cur)){    //If so, then signal of 0 is NOT sent. Go on to pick more signals. If not (prob < ann output), then keep signal = 0;
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
						}*/


					//Old fitness functions
					/*					double payoffReceiver = -1.0;
					bool response = 0;
					bool nullResponse = 0;
					if (g <= nullHonestBeginG){ //Null receiver behavior. Receiver gets payoff from response. Sender gets payoff from nullResponse
						//Just send random s - not actual sender s. q_cur is random
						if (prob(rng) < p*s){    //Does NULL receiver respond to signal? //Parameter p which modifies this. Used for additive fitness; if p = 1.0, then it is always best for senders to send s = 1 during null phase (additive fitness)
							nullResponse = 1;		//null response is response for SENDER
						}
						if (p == 1.0){ //p should be 0.5 if you want to do threshold additive fitness
							double Pr = R_cur.annR_output(q_cur);	//Use different benefit function for now to make Pr = s for receivers
							payoffReceiver = receiver_benefit_function_PrEqualsS(Pr, q_cur);	//Fitness function doesn't matter for this - this is just to get receiver phenotype we want
						} else { // If q != 1.0, we let receivers evolve as if s = q
							if (prob(rng) < R_cur.annR_output(q_cur)){    //If so - respond to signal
								response = 1;
							}
							payoffReceiver = receiver_benefit_function(response, q_cur, d);
						}
					} else {	//Real receiver behavior
						//Now check receiver ANN - do they respond to signal?
						//For real evolution - not honest start
						if (prob(rng) < R_cur.annR_output(s)){    //If so - respond to signal
							response = 1;
						}
						payoffReceiver = receiver_benefit_function(response, q_cur, d);
					}

					//Receiver fitness
					ReceiverPopulation[nullVecR[j]].change_fitness(payoffReceiver);

					//Benefit to sender
					//Now - fitnessFunction (additive = 0, multiplicative = 1) matters

					if (fitnessFunction == 0){ //Additive fitness
						double payoffSender = 0.0;
						if (g <= nullHonestBeginG){ //Null behavior
							payoffSender = sender_fitness_function_Additive(nullResponse, s, q_cur, c, interactionPartners);
						} else { //not null
							payoffSender = sender_fitness_function_Additive(response, s, q_cur, c, interactionPartners);
						}
						//Fitness change
						SenderPopulation[nullVecS[j]].change_fitness(payoffSender);

					} else { //Multiplicative fitness
						//Up to the last interaction, just sum benefits and costs. Then calculate fitness in last interaction
						if (g > nullHonestBeginG){ //not null
							SenderPopulation[nullVecS[j]].incrementBenefit(response);
							SenderPopulation[nullVecS[j]].incrementCost(q_cur, s);

							if (i == (interactionPartners - 1)) {
								//Calculate fitness
								SenderPopulation[nullVecS[j]].setMultFitness(interactionPartners);
							}

						} else { //Null behavior. This doesn't produce as nice of a phenotype as I would like.
							SenderPopulation[nullVecS[j]].incrementBenefit(nullResponse);
							SenderPopulation[nullVecS[j]].incrementCost(q_cur, s);

							if (i == (interactionPartners - 1)) {
								//Calculate fitness
								SenderPopulation[nullVecS[j]].setMultFitness(interactionPartners);
							}
						}
					}
					 */
				}//End loop for this individual

			}//End loop for interaction partners
			//Now all fitnesses have been determined.

			//k-tournament selection

			for (int n = 0; n < N; n++){	//Reproduction, mutation

				tourn_arr[0] = randN(rng);
				double maxFit = SenderPopulation[tourn_arr[0]].get_fitness();
				int maxFitInd = 0;

				for (int k_cur = 1; k_cur < k; k_cur++){
					tourn_arr[k_cur] = randN(rng);

					while (tourn_arr[1] == tourn_arr[0]){ //This will only remove duplicates for k = 2
						tourn_arr[1] = randN(rng);
					}

					if (SenderPopulation[tourn_arr[k_cur]].get_fitness() > maxFit){
						maxFit = SenderPopulation[tourn_arr[k_cur]].get_fitness();
						maxFitInd = k_cur;
					}
				}

				SenderOffspring[n] = SenderPopulation[tourn_arr[maxFitInd]];
				//if (fitnessFunction == 0){
				SenderOffspring[n].reset_fitness();
				//} else {
				//	SenderOffspring[n].reset_fitnessAndCostBenefit();
				//}

				//Determine quality - discrete
				SenderOffspring[n].set_quality(q_vals[randQ(rng)]);

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

				//Mutation for ReceiverOffspring[n]
				for (int m = 0; m <= 38; m++){
					if (ann_mu_R(rng)){	//Mutation to this weight occurs
						//which ann variable to mutate, size of mutation
						ReceiverOffspring[n].ann_mutate(m, ann_mu_size_R(rng));
					}
				}
			}

			if (Report_annVar > 0 & g > 0 & g%ReportFreq_annVar == 0){	//Report ANN stats for fittest individuals

				if (recordFittestANNs == true){
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
						annVars<<"\n"<<rep<<","<<g<<",Receiver,"<<i<<","<<0<<","<<ReceiverPopulation[i].get_fitness();
						for (int m = 0; m <= 38; m++){
							annVars<<","<<ReceiverPopulation[i].get_ann(m);
						}
					}
				} else { //recordFittestANNs==false
					for (int i = 0; i < Report_annVar_N; i++){
						annVars<<"\n"<<rep<<","<<g<<",Sender,"<<i<<","<<SenderPopulation[i].get_quality()<<","<<SenderPopulation[i].get_fitness();
						for (int m = 0; m <= 38; m++){
							annVars<<","<<SenderPopulation[i].get_ann(m);
						}
					}
					for (int i = 0; i < Report_annVar_N; i++){
						annVars<<"\n"<<rep<<","<<g<<",Receiver,"<<i<<","<<0<<","<<ReceiverPopulation[i].get_fitness();
						for (int m = 0; m <= 38; m++){
							annVars<<","<<ReceiverPopulation[i].get_ann(m);
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
		for (int i = 0; i < 5; i++){
			for (int j = 0; j <=33; j++){
				std::cout << ReceiverPopulation[i].get_ann(j) << " ";
			}
			std::cout << "\n";
		}

	}//End replicate loop

	//	dataLog.close();
	paramFile.close();
	annVars.close();


	std::cout << "\nDone";

	return 0;
}
