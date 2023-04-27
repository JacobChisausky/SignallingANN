//============================================================================
// Name        : SignallingANN.cpp
// Author      : Jacob Chisausky
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include <iostream>
#include <vector>
#include <algorithm>
#include <random>
#include <ctime>
#include <fstream>

#include "agents.h"
//using namespace std;

//These cost and benefit functions, and calculation of fitness, are what I really want to go over with Graeme

//Cost function
double cost_function(double s, double q){			//Add options to input too...
	//For now - just linear
	return q - s;
}

//Benefit function for senders
double sender_benefit_function(bool response, double q){
	//For now - flat (quality doesn't impact benefit)
	//Receiver response = 1 fitness. No response = 0.
	return double(response);
}

//Benefit function for receivers
double receiver_benefit_function(bool response, double q){
	//If response sent: payoff = q
	//IF response not sent: payoff = 1-q
	//So payoffs are always positive - but it is better to ignore q < 0.5 and respond to q > 0.5

	//Response sent			//Response not sent (response = 0)
	return (double(response)*q  +  -1*(double(response) - 1.0)*(1-q) );
}

int main() {

	const int k = 2;

	int seed = 123456789;
	double N = 500;	//There will be N receivers and N senders. Stored as double for calculations
	int G = 1000;

	double init_ann_range = 1.0;	//ANN stats will be initialized randomly + or - this value

	double mut_rate_ann_S = 0.01;
	double mut_rate_ann_R = 0.01;

	double mut_step_ann_R = 0.001;
	double mut_step_ann_S = 0.001;


	int s_max = 10;			//The number of signals strengths tried before a signal of 0 is sent by default

	int interactionPartners = 10;

	int replicates = 1;

	int coutReport = 0;

	int Report_annVar = 4; //Export this many generations of ANN variable data. 0 for no report, 1 for only last generation
	int Report_annVar_N = 25;	//How many individuals do you want to report ann stats for? The highest fitness individuals will be reported
	bool Report_annInit = true;	//Do you want initial generation ann data?

	int reportFreq = 50; //Export data every this many generations... change that

	std::string dataFileName = "ann_model";
	std::string dataFileFolder = "C:/Users/owner/Documents/S4/Simulation_ANN";





	//_________End parameter input

	//Reports
	int ReportFreq_annVar = 0;
	if (Report_annVar > 0){
		ReportFreq_annVar = floor(N/Report_annVar);
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
	std::string str3 = dataFileFolder + "/" + strTime + "_summaryStats_" + dataFileName + ".csv";

	//Files to write to
	std::ofstream annVars;
	//	std::ofstream dataLog;
	//	std::ofstream params;
	//	std::ofstream summaryStats;

	annVars.open(str1);
	//		dataLog.open(str1);
	//		params.open(str2);
	//Prepare output files
	annVars << "rep,gen,indType,indNum,quality,fitness,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38";
	//		dataLog << "rep,gen,ind,indType,sendType,strategy,alphaBeta,fitness";
	//		params << "seed,N,G,c1,c2,v1,v2,w1,w2,w3,w4,m,interactionPartners,mutRateAlpha,mutRateBeta,mutRateStrategySender,mutRateStrategyReceiver,mutStepAlpha,mutStepBeta,initializationType,initStrategySender,initStrategyReceiver,initAlpha,initBeta,replicates,alphaBetaMutation,cauchyDist";
	//		params << "\n" << to_string(seed) << "," << to_string(N) <<","<< to_string(G) <<","<< to_string(c1) <<","<< to_string(c2) <<","<< to_string(v1)<<","<<to_string(v2)<<","<<to_string(w1)<<","<<to_string(w2)<<","<<to_string(w3)<<","<<to_string(w4)<<","<<to_string(m)<<","<<to_string(interactionPartners)<<","<<to_string(mutRateAlpha)<<","<<to_string(mutRateBeta)<<","<<to_string(mutRateStrategySender)<<","<<to_string(mutRateStrategyReceiver)<<","<<to_string(mutStepAlpha)<<","<<to_string(mutStepBeta)<<","<<initializationType<<","<<to_string(initStrategySender)<<","<<to_string(initStrategyReceiver)<<","<<to_string(initAlpha)<<","<<to_string(initBeta)<<","<<to_string(replicates)<<","<<alphaBetaMutation<<","<<to_string(cauchyDist),"\n";
	//		if (computeMeansInCpp == true){
	//			summaryStats << "rep,gen,indType,stratNum,stratType,meanAlphaBeta,meanFit,expAlphaBeta,seed,N,G,c1,c2,v1,v2,w1,w2,w3,w4,m,interactionPartners,mutRateAlpha,mutRateBeta,mutRateStrategySender,mutRateStrategyReceiver,mutStepAlpha,mutStepBeta,initializationType,initStrategySender,initStrategyReceiver,initAlpha,initBeta,replicates,alphaBetaMutation,cauchyDist";
	//		}

	//For k-selection tournament
	std::array<int, k> tourn_arr;

	//Random number generators and distributions
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

		for (int i = 0; i < N; i++){
			std::array<double, 39> sender_ann; //Create ann for sender
			for (int j = 0; j < size(sender_ann); j++){
				sender_ann[j] = init_ann_dist(rng);				//Create a random neural network for senders
			}
			//For now - flat quality distribution only
			double qualityInit = qual_dist(rng);

			SenderPopulation.push_back(Sender(sender_ann,qualityInit));		//Create sender agent  with that network
			SenderOffspring	.push_back(Sender(sender_ann,qualityInit));		//Offspring = parents in first generation.

			std::array<double, 34> receiver_ann; //Create ann for receiver
			for (int j = 0; j < size(receiver_ann); j++){		//Create a random neural network for receivers
				receiver_ann[j] = init_ann_dist(rng);
			}
			ReceiverPopulation.push_back(Receiver(receiver_ann));		//Create sender agent  with that network
			ReceiverOffspring.push_back(Receiver(receiver_ann));		//Offspring = parents in first generation.
		}

		//Now SenderPopulation and ReceiverPopulation are populated with agents of random networks and fitness 0.

		//Testing purposes...
		for (int i = 0; i < 10; i++){
			for (int j = 0; j <=33; j++){
				std::cout << ReceiverPopulation[i].get_ann(j) << " ";
			}
			std::cout << "\n";
		}

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
		for (int g = 0; g < G; g++){

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

					//First - try s = 0
					if (prob(rng) > S_cur.annS_output(0.0, q_cur)){    //If so, then signal of 0 is NOT sent. Go on to pick more signals. If not (prob < ann output), then keep signal = 0;
						int s_num = 0;
						while (s_num < s_max){
							s_num++;
							double s_try = prob(rng);
							if (prob(rng) < S_cur.annS_output(s_try, q_cur)){
								s = s_try;
								s_num = s_max;
							}
						}	//If we don't pick a signal after s_max tries, keep at s = 0
					}	//By now, s is (0,1)

					//Determine cost of signal
					double cost_S = cost_function(s, q_cur);

					//Fitness = benefit - cost...? additive not multiplicative - is that okay? Talk to GR about it.

					//Now check receiver ANN - do they respond to signal?
					bool response = 0;
					if (prob(rng) < R_cur.annR_output(s)){    //If so - respond to signal
						response = 1;
					}

					//Benefit to sender
					double benefit_S = sender_benefit_function(response, q_cur);

					//Payoffs to senders and receivers
					double payoffReceiver = receiver_benefit_function(response, q_cur);
					double payoffSender = benefit_S - cost_S;

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

			if (Report_annVar > 0 & g%ReportFreq_annVar == 0){

				std::vector<int> fittestSenders;
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

			/*
				if (g%reportFreq == 0){
					if (coutReport == 1){
						for (int i = 0; i < 20; i++){
							cout << SenderPopulation[i].Strategy << " " << SenderPopulation[i].Alpha << " " << SenderPopulation[i].fitness << " | ";
						}
						cout << endl;
						for (int i = 0; i < 20; i++){
							cout << ReceiverPopulation[i].Strategy << " " << ReceiverPopulation[i].Beta << " " << ReceiverPopulation[i].fitness << " | ";
						}
						cout << "\n--------------------------\n";
					}

					if (computeMeansInCpp == false){

						for (int i = 0; i < N; i++){
							dataLog << "\n" << rep << "," << g << "," << i << ",Sender," << SenderPopulation[i].Type << "," << SenderPopulation[i].Strategy << "," << SenderPopulation[i].Alpha << "," << SenderPopulation[i].fitness;
						}

						for (int i = 0; i < N; i++){
							dataLog << "\n" << rep << "," << g << "," << i << ",Receiver,na," << ReceiverPopulation[i].Strategy << "," << ReceiverPopulation[i].Beta << "," << ReceiverPopulation[i].fitness;
						}
					} else {
						//Compute means
						// rep, gen, indType, stratNum, stratType, meanAlphaBeta, meanFit, expAlphaBeta, fileNumber
						int numSS1 = 0;
						int numSS2 = 0;
						int numSS3 = 0;
						int numRS1 = 0;
						int numRS2 = 0;
						int numRS3 = 0;

						long double totalAlpha = 0.0;
						long double totalBeta = 0.0;

						long double totalFitSS1 = 0.0;
						long double totalFitSS2 = 0.0;
						long double totalFitSS3 = 0.0;

						long double totalFitRS1 = 0.0;
						long double totalFitRS2 = 0.0;
						long double totalFitRS3 = 0.0;

						for (int i = 0; i < N; i++){
							if (SenderPopulation[i].Strategy == 1){
								numSS1 += 1;
								totalFitSS1 += SenderPopulation[i].fitness;
							} else if (SenderPopulation[i].Strategy == 2){
								numSS2 += 1;
								totalFitSS2 += SenderPopulation[i].fitness;
							} else {
								numSS3 += 1;
								totalFitSS3 += SenderPopulation[i].fitness;
							}
							if (ReceiverPopulation[i].Strategy == 1){
								numRS1 += 1;
								totalFitRS1 += ReceiverPopulation[i].fitness;
							} else if (ReceiverPopulation[i].Strategy == 2){
								numRS2 += 1;
								totalFitRS2 += ReceiverPopulation[i].fitness;
							} else {
								numRS3 += 1;
								totalFitRS3 += ReceiverPopulation[i].fitness;
							}

							totalAlpha += SenderPopulation[i].Alpha;
							totalBeta += ReceiverPopulation[i].Beta;

							//totalFitS += SenderPopulation[i].fitness;
							//totalFitR += ReceiverPopulation[i].fitness;
						}

						double meanAlpha = totalAlpha/N;
						double meanBeta = totalBeta/N;

						//double meanFitS = totalFitS/N;
						//double meanFitR = totalFitR/N;

						double meanFitSS1 = totalFitSS1/double(numSS1);
						double meanFitSS2 = totalFitSS2/double(numSS2);
						double meanFitSS3 = totalFitSS3/double(numSS3);

						double meanFitRS1 = totalFitRS1/double(numRS1);
						double meanFitRS2 = totalFitRS2/double(numRS2);
						double meanFitRS3 = totalFitRS3/double(numRS3);


						//                      rep          gen    indType    stratNum    stratType   meanAlphaBeta         meanFit         expAlphaBeta | Params |
				//		summaryStats << "\n" << rep << "," << g << ",Sender," << numSS1 << ",strat1," << meanAlpha << "," << meanFitSS1 << "," << expAlpha  <<","<< seed<<","<<N<<","<<G<<","<<c1<<","<<c2<<","<<v1<<","<<v2<<","<<w1<<","<<w2<<","<<w3<<","<<w4<<","<<m<<","<<interactionPartners<<","<<mutRateAlpha<<","<<mutRateBeta<<","<<mutRateStrategySender<<","<<mutRateStrategyReceiver<<","<<mutStepAlpha<<","<<mutStepBeta<<","<<initializationType<<","<<initStrategySender<<","<<initStrategyReceiver<<","<<initAlpha<<","<<initBeta<<","<<replicates<<","<<alphaBetaMutation<<","<<cauchyDist;
				//		summaryStats << "\n" << rep << "," << g << ",Sender," << numSS2 << ",strat2," << meanAlpha << "," << meanFitSS2 << "," << expAlpha  <<","<< seed<<","<<N<<","<<G<<","<<c1<<","<<c2<<","<<v1<<","<<v2<<","<<w1<<","<<w2<<","<<w3<<","<<w4<<","<<m<<","<<interactionPartners<<","<<mutRateAlpha<<","<<mutRateBeta<<","<<mutRateStrategySender<<","<<mutRateStrategyReceiver<<","<<mutStepAlpha<<","<<mutStepBeta<<","<<initializationType<<","<<initStrategySender<<","<<initStrategyReceiver<<","<<initAlpha<<","<<initBeta<<","<<replicates<<","<<alphaBetaMutation<<","<<cauchyDist;
				//		summaryStats << "\n" << rep << "," << g << ",Sender," << numSS3 << ",strat3," << meanAlpha << "," << meanFitSS3 << "," << expAlpha  <<","<< seed<<","<<N<<","<<G<<","<<c1<<","<<c2<<","<<v1<<","<<v2<<","<<w1<<","<<w2<<","<<w3<<","<<w4<<","<<m<<","<<interactionPartners<<","<<mutRateAlpha<<","<<mutRateBeta<<","<<mutRateStrategySender<<","<<mutRateStrategyReceiver<<","<<mutStepAlpha<<","<<mutStepBeta<<","<<initializationType<<","<<initStrategySender<<","<<initStrategyReceiver<<","<<initAlpha<<","<<initBeta<<","<<replicates<<","<<alphaBetaMutation<<","<<cauchyDist;
				//		summaryStats << "\n" << rep << "," << g << ",Receiver," << numRS1 << ",strat1," << meanBeta << "," << meanFitRS1 << "," << expBeta  <<","<< seed<<","<<N<<","<<G<<","<<c1<<","<<c2<<","<<v1<<","<<v2<<","<<w1<<","<<w2<<","<<w3<<","<<w4<<","<<m<<","<<interactionPartners<<","<<mutRateAlpha<<","<<mutRateBeta<<","<<mutRateStrategySender<<","<<mutRateStrategyReceiver<<","<<mutStepAlpha<<","<<mutStepBeta<<","<<initializationType<<","<<initStrategySender<<","<<initStrategyReceiver<<","<<initAlpha<<","<<initBeta<<","<<replicates<<","<<alphaBetaMutation<<","<<cauchyDist;
				//		summaryStats << "\n" << rep << "," << g << ",Receiver," << numRS2 << ",strat2," << meanBeta << "," << meanFitRS2 << "," << expBeta  <<","<< seed<<","<<N<<","<<G<<","<<c1<<","<<c2<<","<<v1<<","<<v2<<","<<w1<<","<<w2<<","<<w3<<","<<w4<<","<<m<<","<<interactionPartners<<","<<mutRateAlpha<<","<<mutRateBeta<<","<<mutRateStrategySender<<","<<mutRateStrategyReceiver<<","<<mutStepAlpha<<","<<mutStepBeta<<","<<initializationType<<","<<initStrategySender<<","<<initStrategyReceiver<<","<<initAlpha<<","<<initBeta<<","<<replicates<<","<<alphaBetaMutation<<","<<cauchyDist;
				//		summaryStats << "\n" << rep << "," << g << ",Receiver," << numRS3 << ",strat3," << meanBeta << "," << meanFitRS3 << "," << expBeta  <<","<< seed<<","<<N<<","<<G<<","<<c1<<","<<c2<<","<<v1<<","<<v2<<","<<w1<<","<<w2<<","<<w3<<","<<w4<<","<<m<<","<<interactionPartners<<","<<mutRateAlpha<<","<<mutRateBeta<<","<<mutRateStrategySender<<","<<mutRateStrategyReceiver<<","<<mutStepAlpha<<","<<mutStepBeta<<","<<initializationType<<","<<initStrategySender<<","<<initStrategyReceiver<<","<<initAlpha<<","<<initBeta<<","<<replicates<<","<<alphaBetaMutation<<","<<cauchyDist;

						//Write means to summaryStats log.
					}
				}
			 */

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
	//	params.close();
	//	summaryStats.close();





	std::cout << "\nDone";

	return 0;
}
