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

std::vector<LandmarkObs> transformObs(Particle particle, std::vector<LandmarkObs> obs);

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of
	//   x, y, theta and their uncertainties from GPS) and all weights to 1.
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).

  num_particles = 50;
  particles.resize(num_particles);
  weights.resize(num_particles);

  default_random_engine rnd;
  normal_distribution<double> gen_x(x, std[0]);
  normal_distribution<double> gen_y(y, std[1]);
  normal_distribution<double> gen_theta(theta, std[2]);

  for (Particle& prtcl : particles)
  {
    prtcl.x = gen_x(rnd);
    prtcl.y = gen_y(rnd);
    prtcl.theta = gen_theta(rnd);
    prtcl.weight = 1;
  }

}

void ParticleFilter::prediction(double dt, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

  default_random_engine rnd;
  normal_distribution<double> noise_x(0, std_pos[0]);
  normal_distribution<double> noise_y(0, std_pos[1]);
  normal_distribution<double> noise_theta(0, std_pos[2]);

  for (Particle& prtcl : particles)
  {
    double yaw = prtcl.theta;

    // Update particle location/heading
    if (yaw_rate < 0.001) {
      prtcl.x += velocity * cos(yaw) * dt;
      prtcl.y += velocity * sin(yaw) * dt;
    }
    else {
      prtcl.x += (velocity / yaw_rate) * (sin(yaw + yaw_rate * dt) - sin(yaw));
      prtcl.y += (velocity / yaw_rate) * (cos(yaw) - cos(yaw + yaw_rate * dt));
    }

    prtcl.theta += yaw_rate * dt;

    // Add noise to particle location/heading
    prtcl.x += noise_x(rnd);
    prtcl.y += noise_y(rnd);
    prtcl.theta += noise_theta(rnd);
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


  int idx;
  double weight;
  double max_dist = sensor_range*sensor_range;
  double min_dist;

  for (unsigned int p=0; p<num_particles; p++) {
    // observations relative to current particle
    auto prtObsrvtns = transformObs(particles[p], observations); // vector of landmarks

    // create association and obervation lists
    particles[p].associations.clear();
  	particles[p].sense_x.clear();
  	particles[p].sense_y.clear();

    // initialize particle's new weight
    weight = 1;

    // Get list of landmarks within range of particle
    vector<LandmarkObs> lm_inR;
    for (unsigned int ii=0; ii<map_landmarks.landmark_list.size(); ii++) {
      double particle_landmark_dist = dist2(particles[p].x,particles[p].y, map_landmarks.landmark_list[ii].x_f, map_landmarks.landmark_list[ii].y_f);
      if (particle_landmark_dist < max_dist) {
        LandmarkObs lm;
        lm.x = map_landmarks.landmark_list[ii].x_f;
        lm.y = map_landmarks.landmark_list[ii].y_f;
        lm.id = map_landmarks.landmark_list[ii].id_i;
        lm_inR.push_back(lm);
      }
    }

    // Find the closest map landmark to the observations relative to the current
    // particle
    for (LandmarkObs plm : prtObsrvtns) {
      min_dist = max_dist;

      // iterate through each landmark in range of the particle
      for (unsigned int ii=0; ii<lm_inR.size(); ii++) {
        double tmp_dist = dist2(plm.x, plm.y, lm_inR[ii].x, lm_inR[ii].y);
        if (tmp_dist < min_dist) {
          min_dist = tmp_dist; // update the minimum distance
          idx = ii; // save the index of the current landmark
        }
      }

      if (min_dist < max_dist) {
        // update the particle's weight using the probability of the observation
        // being associated with the map landmark
        weight *= Nxy(plm.x, plm.y, lm_inR[idx].x, lm_inR[idx].y, std_landmark[0], std_landmark[1]);

        // save associations for plotting.
        particles[p].associations.push_back(lm_inR[idx].id);
        particles[p].sense_x.push_back(plm.x);
        particles[p].sense_y.push_back(plm.y);
      } else {
        weight = 0;
      }
    }

    // update current particles weight
    particles[p].weight = weight;
    weights[p] = weight;
  }
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight.
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

  // create a discrete distribution for re-sampling particles
  default_random_engine rnd;
  discrete_distribution<int> probableParticleIdx(weights.begin(), weights.end());

  // resample particles
  vector<Particle> newParticles;
  newParticles.resize(num_particles);

  // auto newParticles = particles;
  for (unsigned int p=0; p<num_particles; p++)
    newParticles[p] = particles[probableParticleIdx(rnd)];

  particles = newParticles;
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

std::vector<LandmarkObs> transformObs(Particle prtcl, std::vector<LandmarkObs> obs)
{
  std::vector<LandmarkObs> transObs;
  transObs.resize(obs.size());

  for (unsigned int idx=0; idx < obs.size(); idx++)
  {
    double cosTh = cos(prtcl.theta);
    double sinTh = sin(prtcl.theta);
    transObs[idx].x = prtcl.x + obs[idx].x * cosTh - obs[idx].y * sinTh;
    transObs[idx].y = prtcl.y + obs[idx].x * sinTh + obs[idx].y * cosTh;
  }

  return transObs;
}
