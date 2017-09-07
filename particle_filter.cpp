/*
 * particle_filter.cpp
 *
 *  Created on: September 06, 2017
 *      Author: Tanuj Kush
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h> 
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// Set the number of particles. Initialize all particles to first position (based on estimates of 
	// x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	
	// Set up a random generator
	default_random_engine gen;

	// Set up a normal distribution for x, y and theta
	normal_distribution<double> dist_x(x, std[0]);
	normal_distribution<double> dist_y(y, std[1]);
	normal_distribution<double> dist_theta(theta, std[2]);

	// Number of particles in the filter
	num_particles = 100;
	
	// Add particles, initialized to initial GPS coordinates plus some noise
	for (int i=0;i<num_particles;i++)
	{
		Particle p;
		p.x = dist_x(gen);
		p.y = dist_y(gen);
		//p.weight = 1.;
		p.theta = dist_theta(gen);
		weights.push_back(1.);

		//theta is measured positive CCW and is bound by [0,2*pi)
		//if theta is out of this range, it must be wrapped back
		wrap_theta(p.theta);

		particles.push_back(p);
	}

	//initialization is done
	is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// Predict the next location for each particle based on motion model and GPS sensor noise
	
	// Setup a random generator
	default_random_engine gen;

	double x, y, theta;
	double x_new, y_new, theta_new;

	// update the positions and heading based on motion model + sensor noise
	for (int i=0;i<num_particles;i++)
	{
		x = particles[i].x;
		y = particles[i].y;
		theta = particles[i].theta;

		// EOM for almost zero yaw rate
		if ( fabs(yaw_rate) < EPSILON )
		{
			x_new = x + ( velocity*delta_t )*cos( theta );
			y_new = y + ( velocity*delta_t )*sin( theta );
			theta_new = theta;	
		}
		// EOM for non zero yaw rate
		else
		{
			x_new = x + ( velocity/yaw_rate )*(  sin( theta + yaw_rate*delta_t ) - sin( theta ) );
			y_new = y + ( velocity/yaw_rate )*( -cos( theta + yaw_rate*delta_t ) + cos( theta ) );
			theta_new = theta + yaw_rate*delta_t;
			wrap_theta(theta_new);
		}
		
		// Generate Gaussian noise for GPS sensor noise
		normal_distribution<double> dist_x(x_new, std_pos[0]);
		normal_distribution<double> dist_y(y_new, std_pos[1]);
		normal_distribution<double> dist_theta(theta_new, std_pos[2]);
		
		// Set particle data
		particles[i].x = dist_x(gen);
		particles[i].y = dist_y(gen);
		particles[i].theta = dist_theta(gen);

		wrap_theta(particles[i].theta);
	}
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs>& predicted, Map map_lm) {
	// Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.

	// loop over all the transformed observations
	for (int i=0;i<predicted.size();i++)
	{
		double minDist = 1.e10;
		double dist = 1.e10;
		int minId = -1;
	
		// loop over all the landmark information
		for (int j=0;j<map_lm.landmark_list.size();j++)
		{
			dist = (map_lm.landmark_list[j].x_f - predicted[i].x)*(map_lm.landmark_list[j].x_f - predicted[i].x) \
				 + (map_lm.landmark_list[j].y_f - predicted[i].y)*(map_lm.landmark_list[j].y_f - predicted[i].y);
			
			// find minimum distance (no need for sqrt here, this speeds things up a bit!)
			if ( dist < minDist )
			{
				minDist = dist;
				minId = j;
			}
		}

		//associate transformed observation with landmark id minId
		predicted[i].id = map_lm.landmark_list[minId].id_i;		
	}
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		std::vector<LandmarkObs> observations, Map map_landmarks) {
	//   Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html

    double dist = 0.;
	double weight_sum = 0.;
	double exponent = 0.;
	double divisor = 2 * M_PI * std_landmark[0] * std_landmark[1];
			
	double x, y;
	double mu_x, mu_y;
	
	// loop over all particles
	for (int i=0;i<num_particles;i++)
	{
		double weight_prod = 1.;

		// create a vector of transformed obs
		std::vector<LandmarkObs> transformed_obs;
		for (int j=0;j<observations.size();j++)
		{
			LandmarkObs lm_pred;
			// transform car's observations to particle observations
			lm_pred.x = particles[i].x + ( observations[j].x * cos( particles[i].theta ) - observations[j].y * sin ( particles[i].theta ) );
			lm_pred.y = particles[i].y + ( observations[j].x * sin( particles[i].theta ) + observations[j].y * cos ( particles[i].theta ) );

			// check if this observation is within sensor range
			// if not, then ignore
			dist = lm_pred.x * lm_pred.x + lm_pred.y * lm_pred.y;
			if ( true ) //dist < sensor_range * sensor_range )
				transformed_obs.push_back(lm_pred);
		}
		
		// associate transformed observations with map landmarks
		dataAssociation(transformed_obs, map_landmarks);
		
		// update weights based on observations and associations
		for (int j=0;j<transformed_obs.size();j++)
		{
			x = transformed_obs[j].x;
			y = transformed_obs[j].y;
			mu_x = map_landmarks.landmark_list[transformed_obs[j].id-1].x_f;
			mu_y = map_landmarks.landmark_list[transformed_obs[j].id-1].y_f;
			exponent = -( ( x - mu_x )*( x - mu_x ) / ( 2*std_landmark[0]*std_landmark[0] ) + 
				          ( y - mu_y )*( y - mu_y ) / ( 2*std_landmark[1]*std_landmark[1] ) );
						
			// avoid division by zero
			if ( fabs( divisor ) > EPSILON )
				weight_prod *= exp( exponent ) / divisor;
			else
				//ignore zero variance 
				continue;
		}
		
		// update particle weight
		particles[i].weight = weight_prod;
		weights[i] = particles[i].weight;

		// calculate weight sum
		weight_sum += weights[i];	
	}

	// normalize weights for resampling
	// avoid division by zero
	if ( weight_sum < EPSILON*EPSILON*EPSILON )
		cout << "Sum of weights for resampling is almost zero! ("<<weight_sum<<")\n";
	else
	{
		for (int i=0;i<num_particles;i++)
			weights[i] /= weight_sum;
	}
}

void ParticleFilter::resample() {
	// Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
	
	// new vector for particles
	std::vector<Particle> resampled;

	// random generator
	std::default_random_engine gen;

	// discrete distribution based on weights vector
    std::discrete_distribution<> distribution( weights.begin(), weights.end() );
	int sample_id = -1;
	
	// pick particles based on weights with replacement
	for (int i=0;i<num_particles;i++)
	{
		sample_id = distribution(gen);
		resampled.push_back(particles[sample_id]);
	}
	
	// update particles vector with the sampled one
	particles.clear();
	particles = resampled;
}

Particle ParticleFilter::SetAssociations(Particle particle, std::vector<int> associations, std::vector<double> sense_x, std::vector<double> sense_y)
{
	//particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
	// associations: The landmark id that goes along with each listed association
	// sense_x: the associations x mapping already converted to world coordinates
	// sense_y: the associations y mapping already converted to world coordinates

	//Clear the previous associations
	particle.associations.clear();
	particle.sense_x.clear();
	particle.sense_y.clear();

	particle.associations= associations;
 	particle.sense_x = sense_x;
 	particle.sense_y = sense_y;

 	return particle;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}

void ParticleFilter::printDebugInfo(string method, Particle p)
{
	cout << "In " << method <<"\n";
	cout << "Particle " << ":\n";
	cout << "X:\t " << p.x << "\n";
	cout << "Y:\t " << p.y << "\n";
	cout << "Th:\t " << p.theta << "\n";
	cout << "Wt:\t " << p.weight << "\n";
	cout << "***\n";	
}
