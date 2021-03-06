/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
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
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
    
    num_particles = 10;
    
    default_random_engine gen;
    
    normal_distribution<double> N_x(x, std[0]);
    normal_distribution<double> N_y(y, std[1]);
    normal_distribution<double> N_theta(theta, std[2]);
    
    for(int i = 0; i < num_particles; i++){
        Particle particle;
        particle.id = i;
        particle.x = N_x(gen);
        particle.y = N_y(gen);
        particle.theta = N_theta(gen);
        particle.weight = 1;
        particles.push_back(particle);
        weights.push_back(particle.weight);
    }
    
    is_initialized = true;

}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

    default_random_engine gen;
    
    for(int i = 0; i < num_particles; i++){
        double x_prediction;
        double y_prediction;
        double theta_prediction;
        
        if(yaw_rate == 0){
            x_prediction = particles[i].x + velocity*delta_t*cos(particles[i].theta);
            y_prediction = particles[i].y + velocity*delta_t*sin(particles[i].theta);
            theta_prediction = particles[i].theta;
        }
        else{
            x_prediction = particles[i].x + velocity/yaw_rate*(sin(particles[i].theta+yaw_rate*delta_t)-sin(particles[i].theta));
            y_prediction = particles[i].y + velocity/yaw_rate*(cos(particles[i].theta) - cos(particles[i].theta+yaw_rate*delta_t));
            theta_prediction = particles[i].theta+yaw_rate*delta_t;
        }
        
        normal_distribution<double> N_x(x_prediction,std_pos[0]);
        normal_distribution<double> N_y(y_prediction,std_pos[1]);
        normal_distribution<double> N_theta(theta_prediction,std_pos[2]);
        
        particles[i].x = N_x(gen);
        particles[i].y = N_y(gen);
        particles[i].theta = N_theta(gen);
    }

}


void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

}




void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		std::vector<LandmarkObs> observations, Map map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html
    
    
    weights.clear();
    
    for(int p=0; p<num_particles ; p++){
        
        vector<int> landmarks_in_range;


        for(int i=0;i<map_landmarks.landmark_list.size();i++){

            double lm_distance = dist(particles[p].x, particles[p].y,map_landmarks.landmark_list[i].x_f,map_landmarks.landmark_list[i].y_f);
            
            if(lm_distance <= sensor_range){
                landmarks_in_range.push_back(i);
            }
        }

        if (landmarks_in_range.size() == 0){
            cout << "No landmarks in range" << endl;
        }
        
        
        long double total_weight = 1.0;
        
        for(int i = 0; i < observations.size(); i++){
            LandmarkObs trans_obs;
            
            trans_obs.x = particles[p].x + (observations[i].x*cos(particles[p].theta) - observations[i].y*sin(particles[p].theta));
            trans_obs.y = particles[p].y + (observations[i].x*sin(particles[p].theta) + observations[i].y*cos(particles[p].theta));
            
            int best_landmark = landmarks_in_range[0];
            double best_distance = dist(trans_obs.x,trans_obs.y,map_landmarks.landmark_list[best_landmark].x_f,map_landmarks.landmark_list[best_landmark].y_f);
                
            for( int l = 1 ; l < landmarks_in_range.size() ; l++){
                
                int lm = landmarks_in_range[l];
                    
                double this_distance = dist(trans_obs.x,trans_obs.y,map_landmarks.landmark_list[lm].x_f,map_landmarks.landmark_list[lm].y_f);
                
                if(this_distance < best_distance){
                    best_landmark = lm;
                    best_distance = this_distance;
                }
                    
            }
            
            trans_obs.id = best_landmark;
            
            
            double observation_weight;
            
            double meas_x = trans_obs.x;
            double meas_y = trans_obs.y;
            int landmark_id = trans_obs.id;
            double mu_x = map_landmarks.landmark_list[landmark_id].x_f;
            double mu_y = map_landmarks.landmark_list[landmark_id].y_f;
            
            double sig_x = std_landmark[0];
            double sig_y = std_landmark[1];
            
            observation_weight = 1/(2*M_PI*sig_x*sig_y)*exp(-(pow(meas_x-mu_x,2)/(2*sig_x*sig_x) + pow(meas_y-mu_y,2)/(2*sig_y*sig_y)));
            
            if (observation_weight > 0){
                total_weight = total_weight * observation_weight;
            }
            
        
        }
        
        particles[p].weight = total_weight;
        weights.push_back(particles[p].weight);
    }
    
    
}


void ParticleFilter::resample() {
    // TODO: Resample particles with replacement with probability proportional to their weight.
    // NOTE: You may find std::discrete_distribution helpful here.
    //   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
    
    default_random_engine gen;
    discrete_distribution<int> distribution(weights.begin(), weights.end());
    
    vector<Particle> resampled;
    vector<double> new_weights;
    
    for(int i=0;i<num_particles;i++){
        resampled.push_back(particles[distribution(gen)]);
    }
    
    particles = resampled;
    
    for(int i=0;i<particles.size();i++){
        new_weights.push_back(particles[i].weight);
    }
    
    weights = new_weights;
    
    
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
