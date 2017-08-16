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

  num_particles = 100;
  particles.resize(num_particles);

  default_random_engine rnd;
  normal_distribution gen_x(x, std[0]);
  normal_distribution gen_y(y, std[1]);
  normal_distribution gen_theta(theta, std[2]);

  // for (unsigned int p=0; p < num_particles; p++)
  // {
  //   particles[p].x = gen_x(rnd);
  //   particles[p].y = gen_y(rnd);
  //   particles[p].theta = gen_theta(rnd);
  //   particles[p].weight = 1;
  // }
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
  normal_distribution noise_x(0, std_pos[0]);
  normal_distribution noise_y(0, std_pos[1]);
  normal_distribution noise_theta(0, std_pos[2]);

  for (Particle& prtcl : particles)
  {
    double yaw = prtcl.theta;
    double cosTh = cos(yaw);
    double sinTh = sin(yaw);

    // Update particle location/heading
    prtcl.x += (yaw_rate < 0.001) ? velocity * cosTh * dt :
      (velocity / yaw_rate) * (sin(yaw + yaw_rate * dt) - sinTh);

    prtcl.y += (yaw_rate < 0.001) ? velocity * sinTh * dt :
      (velocity / yaw_rate) * (cosTh - cos(yaw + yaw_rate * dt));

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

  for (Particle & prtcl : particles) {
    // observations relative to current particle
    auto prtclObsrvtns = transformObs(prtcl, observations); // vector of landmarks

    // make copy of landmark list
    auto landmark_list = map_landmarks.landmark_list;

    // create association and obervation lists
    vector<int> assoc;
    vector<double> sense_x;
    vector<double> sense_y;

    // initialize particle's new weight
    double weight = 1;

    // Find the closest map landmark to the observations relative to the current
    // particle
    for (auto & plm : prtclObsrvtns) {
      double min_dist = sensor_range;
      int idx;

      // iterate through each landmark in the map.
      for (unsigned int ii=0; ii<landmark_list.size(); ii++) {
        auto lm = landmark_list[ii];
        double tmp_dist = dist2(plm.x, plm.y, landmark_list[ii].x_f, landmark_list[ii].y_f);
        if (tmp_dist < min_dist) {
          // update the minimum distance
          min_dist = tmp_dist;
          // save the index of the current landmark
          idx = ii;
        }
      }

      if (min_dist < sensor_range) {
        // save map landmark ID of closest landmark
        assoc.push_back(landmark_idx[idx].id_i);

        // save xy location of observation closest to a map landmark
        sense_x.push_back(plm.x);
        sense_y.push_back(plm.y);

        // update the particle's weight using the probability of the observation
        // being associated with the map landmark
        weight *= Nxy(plm.x,
                      plm.y,
                      landmark_list[idx].x_f,
                      landmark_list[idx].y_f,
                      std_landmark[0],
                      std_landmark[1]);

        // remove the closest landmark from the list
        //  --> This assumes we can only have one observation for each landmark
        landmark_list.erase(landmark_list.begin() + idx);
      }
    }

    // update current particles weight
    prtcl.weight = weight;

    // set particle assotiations
    prtcl = SetAssociations(prtcl, assoc, sense_x, sense_y);
  }
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight.
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

  // Create array of particle weights
  double prtclWeight[num_particles];

  for (unsigned int p=0; p<num_particles; p++)
    prtclWeight[p] = particles[p].weight;

  // create a discrete distribution for re-sampling particles
  default_random_engine rnd;
  discrete_distribution<int> probableParticleIdx(prtclWeights);

  // resample particles
  auto newParticles = particles;
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
  transObs = obs;

  for (LandmarkObs & lm : transObs)
  {
    double cosTh = cos(-prtcl.theta);
    double sinTh = sin(-prtcl.theta);
    double xn = prtcl.x + lm.x * cosTh - lm.y * sinTh;
    double yn = prtcl.y + lm.x * sinTh + lm.y * cosTh;
    lm.x = xn;
    lm.y = yn;
  }

  return transObs;
}
