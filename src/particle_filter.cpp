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

#define null 0.0001

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).

	// setting number of particles
	num_particles =50;
	
	
	// creating Gaussians
	normal_distribution<double> dist_x(x, std[0]);
	normal_distribution<double> dist_y(y, std[1]);
	normal_distribution<double> dist_theta(theta, std[2]);

	
	default_random_engine gen;
	
	// generating particle set
	for (int i = 0; i < num_particles; i++) {
	  Particle current_particle;
	  current_particle.id = i;
	  current_particle.x = dist_x(gen);
	  current_particle.y = dist_y(gen);
	  current_particle.theta = dist_theta(gen);
	  current_particle.weight = 1.0;
	  
	  particles.push_back(current_particle);
	  weights.push_back(current_particle.weight);
	}
	is_initialized = true;

}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
	

	default_random_engine gen;
	
	
	for (auto& i:particles) {
		double particle_x = i.x;
		double particle_y = i.y;
		double particle_theta = i.theta;
		 
		double pred_x;
		double pred_y;
		double pred_theta;
		 
		if (fabs(yaw_rate) < null) {
			pred_x = particle_x + velocity * cos(particle_theta) * delta_t;
			pred_y = particle_y + velocity * sin(particle_theta) * delta_t;
			pred_theta = particle_theta;
		} else {
			pred_x = particle_x + (velocity/yaw_rate) * (sin(particle_theta + (yaw_rate * delta_t)) - sin(particle_theta));
			pred_y = particle_y + (velocity/yaw_rate) * (cos(particle_theta) - cos(particle_theta + (yaw_rate * delta_t)));
			pred_theta = particle_theta + (yaw_rate * delta_t);
		}
		  
		normal_distribution<double> dist_x(pred_x, std_pos[0]);
		normal_distribution<double> dist_y(pred_y, std_pos[1]);
		normal_distribution<double> dist_theta(pred_theta, std_pos[2]);
		  
		i.x = dist_x(gen);
		i.y = dist_y(gen);
		i.theta = dist_theta(gen);
	}
	
	
		
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.
		
	
	for (auto& i: observations) {
	
		double closest_landmark_dist = INFINITY;
		int closest_landmark_id = -1;
		
		double obs_x = i.x;
		double obs_y = i.y;

		for (auto& j : predicted) {
		  double pred_x = j.x;
		  double pred_y = j.y;
		  int 	 pred_id = j.id;
		  double current_dist = dist(obs_x, obs_y, pred_x, pred_y);
		  

		  if (current_dist < closest_landmark_dist) {
		    closest_landmark_dist = current_dist;
		    closest_landmark_id = pred_id;
		  }
		}
		i.id = closest_landmark_id;
	}
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		std::vector<LandmarkObs> observations, Map map_landmarks) {
	// TODO: Update the weights of each particle using a multi-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html
	

  
  
  double weights_sum = 0.0;
  double sensor_range_doubled = pow(sensor_range,2);

  for (auto&  i : particles) {
    double particle_x = i.x;
    double particle_y = i.y;
    double particle_theta = i.theta;

    
    vector<LandmarkObs> transformed_observations;

    //transforming current observation from vehicle to map coordinates
    for (unsigned int j = 0; j < observations.size(); j++) {
      LandmarkObs current_transf_observation;
      current_transf_observation.id = j;
      current_transf_observation.x = particle_x + (cos(particle_theta) * observations[j].x) - (sin(particle_theta) * observations[j].y);
      current_transf_observation.y = particle_y + (sin(particle_theta) * observations[j].x) + (cos(particle_theta) * observations[j].y);
      transformed_observations.push_back(current_transf_observation);
    }

    // filtering the landmarks - saving only those in the sensor range
    vector<LandmarkObs> predicted_landmarks;
   
   for (auto&  k : map_landmarks.landmark_list) {

	Map::single_landmark_s current_landmark; 
	
	current_landmark = k;
	
	
	double diff_x = particle_x - current_landmark.x_f;
	double diff_y = particle_y - current_landmark.y_f;
	
    if ( diff_x*diff_x + diff_y*diff_y <= sensor_range_doubled ) {
        predicted_landmarks.push_back(LandmarkObs {current_landmark.id_i, current_landmark.x_f, current_landmark.y_f});
      }
    }

    
    //associating observations with predicted (filterered) landmarks
    dataAssociation(predicted_landmarks, transformed_observations);

    // calculating the particle weights 
    i.weight = 1.0;

    double sigma_xx = pow(std_landmark[0], 2);
    double sigma_yy = pow(std_landmark[1], 2);
    double normalizer = (1.0/(2.0 * M_PI * std_landmark[0] * std_landmark[1]));
    
    
    /*Calculate the weight of particle based on the multivariate Gaussian probability function*/
    for (auto& l : transformed_observations) {
      
	  double xy_p = 1.0;
	  
	  double transf_obs_x = l.x;
      double transf_obs_y = l.y;
      double transf_obs_id = l.id;
      

      for (auto& m : predicted_landmarks) {
        double pred_landmark_x = m.x;
        double pred_landmark_y = m.y;
        double pred_landmark_id = m.id;

        if (transf_obs_id == pred_landmark_id) {
          xy_p = normalizer * exp(-1.0 * ((pow((transf_obs_x - pred_landmark_x), 2)/(2.0 * sigma_xx)) + (pow((transf_obs_y - pred_landmark_y), 2)/(2.0 * sigma_yy))));
          i.weight *= xy_p;
        }
      }
    }
    weights_sum += i.weight;
  }

  // normalizing weights according to the sum to [0,1]
   for (unsigned int i = 0; i < particles.size(); i++) {
    particles[i].weight /= weights_sum;
    weights[i] = particles[i].weight;
  }
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
	
	vector<Particle> resampled_particles;

	
	default_random_engine gen;
	
	// wheel algorithm
	uniform_int_distribution<int> particle_index(0, num_particles - 1);
	// random starting particle
	int current_index = particle_index(gen);
	
	double beta = 0.0;
	
	double max_weight_doubled = 2.0 * *max_element(weights.begin(), weights.end());
	
	for (const auto& i : particles) {
		uniform_real_distribution<double> random_weight(0.0, max_weight_doubled);
		beta += random_weight(gen);

	  while (beta > weights[current_index]) {
	    beta -= weights[current_index];
	    current_index = (current_index + 1) % num_particles;
	  }
	  resampled_particles.push_back(particles[current_index]);
	}
	particles = resampled_particles;
}

Particle ParticleFilter::SetAssociations(Particle particle, std::vector<int> associations, std::vector<double> sense_x, std::vector<double> sense_y)
{
	//particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
	// associations: The landmark id that goes along with each listed association
	// sense_x: the associations x mapping already converted to world coordinates
	// sense_y: the associations y mapping already converted to world coordinates

	//clearing the previous associations
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
