/*
 * signaller.cpp
 *
 *  Created on: Apr 26, 2023
 *      Author: owner
 */

#include "agents.h"

Sender::Sender(std::array<double, 39> annS_stats, double qualityInit) :		//Should I use pointers to the array to speed things up..?
	annS(annS_stats), quality(qualityInit)
{
}

double Sender::annS_output(double s, double q){

	double n1 = rlu(annS[1] + s);
	double n2 = rlu(annS[2] + q);

	double n3 = rlu(annS[3] + n1 * annS[11] + n2 * annS[15]);
	double n4 = rlu(annS[4] + n1 * annS[12] + n2 * annS[16]);
	double n5 = rlu(annS[5] + n1 * annS[13] + n2 * annS[17]);
	double n6 = rlu(annS[6] + n1 * annS[14] + n2 * annS[18]);

	double n7 = rlu(annS[7] + n3 * annS[19] + n4 * annS[23] + n5 * annS[27] + n6 * annS[31]);
	double n8 = rlu(annS[8] + n3 * annS[20] + n4 * annS[24] + n5 * annS[28] + n6 * annS[32]);
	double n9 = rlu(annS[9] + n3 * annS[21] + n4 * annS[25] + n5 * annS[29] + n6 * annS[33]);
	double n10 = rlu(annS[10] + n3 * annS[22] + n4 * annS[26] + n5 * annS[30] + n6 * annS[34]);

	double output = rlu(annS[0] + n7 * annS[35] + n8 * annS[36] + n9 * annS[37] + n10 *annS[38]);

	return output;
}

void Sender::ann_mutate(int var, double size){
	annS[var] += size;
	if (annS[var] > 100.0){	//Constrain network variables to be between -100 and 100
		annS[var] = 100.0;
	} else if (annS[var] < -100.0){
		annS[var] = -100.0;
	}
}

void Sender::set_quality(double qual){
	quality = qual;
}

double Sender::get_quality(){
	return quality;
}

double Sender::get_fitness(){
	return fitness;
}

void Sender::change_fitness(double change){
	fitness += change;
}

void Sender::reset_fitness(){
	fitness = 0;
}

double Sender::get_ann(int var){
	return annS[var];
}

Receiver::Receiver(std::array<double, 34> annR_stats) :
	annR(annR_stats)
{
}

double Receiver::annR_output(double s){

	double n1 = rlu(annR[1] + s);

	double n2 = rlu(annR[2] + n1 * annR[10]);
	double n3 = rlu(annR[3] + n1 * annR[11]);
	double n4 = rlu(annR[4] + n1 * annR[12]);
	double n5 = rlu(annR[5] + n1 * annR[13]);

	double n6 = rlu(annR[6] + n2 * annR[14] + n3 * annR[18] + n4 * annR[22] + n5 * annR[26]);
	double n7 = rlu(annR[7] + n2 * annR[15] + n3 + annR[19] + n4 * annR[23] + n5 * annR[27]);
	double n8 = rlu(annR[8] + n2 * annR[16] + n3 * annR[20] + n4 * annR[24] + n5 * annR[28]);
	double n9 = rlu(annR[9] + n2 * annR[17] + n3 * annR[21] + n4 * annR[25] + n5 * annR[29]);

	double output = rlu(annR[0] + n6 * annR[30] + n7 * annR[31] + n8 * annR[32] + n9 * annR[33]);

	return output;
}

void Receiver::ann_mutate(int var, double size){
	annR[var] += size;
	if (annR[var] > 100.0){	//Constrain network variables to be between -100 and 100
		annR[var] = 100.0;
	} else if (annR[var] < -100.0){
		annR[var] = -100.0;
	}
}

void Receiver::change_fitness(double change){
	fitness += change;
}

double Receiver::get_fitness(){
	return fitness;
}

void Receiver::reset_fitness(){
	fitness = 0;
}

double Receiver::get_ann(int var){
	return annR[var];
}


double rlu(double input){
	return std::min(std::max(0.0,input),1.0);
}