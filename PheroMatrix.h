#pragma once
#include <vector>
#include <float.h>

typedef std::vector<std::vector<double> > MatrixArrayTypeDouble;


class PheroMatrix
{
public:

	double evaporation_rate_ , initial_pheromone_;
	MatrixArrayTypeDouble *matrix_;
	~PheroMatrix(void){
		delete matrix_;
	}
	PheroMatrix(int size, double evaporation_rate, double initial_pheromone){
		matrix_ = new MatrixArrayTypeDouble(size , size);
		evaporation_rate_ = evaporation_rate;
		initial_pheromone_ = initial_pheromone;
	}
	PheroMatrix(int size, double evaporation_rate){
		matrix_ = new MatrixArrayTypeDouble(size , size);
		evaporation_rate_ = evaporation_rate;
		initial_pheromone_ = 1/size;
		for (size_t i=0;i<matrix_->size();i++)	{
			for (size_t j=0;j<matrix_->size();j++)
				(*matrix_)[i][j] = initial_pheromone_;// 0.1 * rand() / (double)RAND_MAX;
		}
	}
	double get(unsigned int v, unsigned int w) {
		return (*matrix_)[v][w];
	}
	void add(unsigned int v, unsigned int w, double amount) {
		(*matrix_)[v][w] += amount;
		(*matrix_)[w][v] += amount;
	}
	void multipy(unsigned int v, unsigned int w, double amount) {
		(*matrix_)[v][w] *= amount;
		(*matrix_)[w][v] *= amount;
	}
	void subtract(unsigned int v, unsigned int w, double amount){
		(*matrix_)[v][w] -= amount;
		(*matrix_)[w][v] -= amount;
	}
	void evaporate(unsigned int v, unsigned int w) {
		(*matrix_)[v][w] *= 1 - evaporation_rate_;
		(*matrix_)[w][v] = (*matrix_)[v][w];
	}
	void set(unsigned int v, unsigned int w , double value) {
		(*matrix_)[v][w] = value;
		(*matrix_)[w][v] = value;
	}
	void evaporate_all() {
		for(unsigned int i=0;i<matrix_->size();i++) {
			for(unsigned int j=0;j<matrix_->size();j++) {
				evaporate(i,j);
			}
		}
	}
	double get_evaporation_rate() {
		return evaporation_rate_;
	}
	unsigned int  size() {
		return static_cast<unsigned int>(matrix_->size());
	}
	double lambda_branching_factor(unsigned int v, double lambda) {
		double min_pheromone = std::numeric_limits<double>::max();
		double max_pheromone = 0.0;
		for(unsigned int i=0;i<this->size();i++) {
			double pheromone = (*matrix_)[v][i];
			if(min_pheromone > pheromone) {
				min_pheromone = pheromone;
			}
			if(max_pheromone < pheromone) {
				max_pheromone = pheromone;
			}
		}
		double limit = min_pheromone + lambda * (max_pheromone - min_pheromone);
		unsigned int branching_factor = 0;
		for(unsigned int j=0;j<this->size();j++) {
			if((*matrix_)[v][j] >= limit) {
				branching_factor++;
			}
		}
		return branching_factor;
	}

	double average_lambda_branching_factor(double lambda) {
		double sum = 0.0;
		for(unsigned int i=0;i<this->size();i++) {
			sum += lambda_branching_factor(i, lambda);
		}
		return sum / (this->size() *2);
	}
};

