/*
 * signaller.h
 *
 *  Created on: Apr 26, 2023
 *      Author: owner
 */

//#ifndef SIGNALLER_H_
//#define SIGNALLER_H_
#pragma once
#include <array>

class Sender {
public:
	Sender(const std::array<double, 39> annS_stats, double qualityInit);
	double annS_output(double s, double q);
	void change_fitness(double change);
	void reset_fitness();
	void set_quality(double qual);
	double get_quality();
	double get_fitness();
	void ann_mutate(int var, double size);
	double get_ann(int var);

private:
	std::array<double, 39> annS; //weights and biases for annS
	double fitness = 0.0;
	double quality;

};

class Receiver {
public:
	Receiver(const std::array<double, 34> annR_stats);
	double annR_output(double s);
	void change_fitness(double change);
	void reset_fitness();
	double get_fitness();
	void ann_mutate(int var, double size);
	double get_ann(int var);

private:
	std::array<double, 34> annR;//weights and biases for annR
	double fitness = 0.0;

};

double rlu(double input);
bool annS_test(int resolution, std::array<double, 39> ann);
bool annR_test(int resolution, std::array<double, 34> ann);


//#endif /* SIGNALLER_H_ */
