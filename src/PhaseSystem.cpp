/*************************************************************************
	> File Name: PhaseSystem.cpp
	> Author: Yi Zheng 
	> Mail: hczhengcq@gmail.com
	> Created Time: Wed 09 May 2018 02:26:58 PM CEST
 ************************************************************************/

#include <iostream>
#include <string>
#include <vector>
#include <set>
#include <map>
#include <algorithm>
#include <cmath>
#include <cassert>
#include <stdexcept>
#include <ProMP/PhaseSystem.hpp>
using namespace std;

namespace ProMP
{
    PhaseSystem::PhaseSystem(int num_basis, double width, const double goal_t, const double start_t) : num_basis_(num_basis), width_(width), z_(start_t) 
    {
        if (start_t < 0)
            throw std::invalid_argument("The start time point of the trajectory should >= 0!");
        if (goal_t <= start_t)
            throw std::invalid_argument("Goal must be chronologically after the start!");
        if (num_basis <= 1)
            throw std::invalid_argument("The number of basis function should > 1!");

        traj_duration_ = goal_t - start_t;
        z_dot_ = 1.0 / traj_duration_;
        execute_ = false;
    }
    
    /**Generate the vector for the phase system
     */
    void PhaseSystem::init()
    {
        double dis = 1 / (num_basis_ - 1);
        center_vec_(0) = 0.0;
        for (int i = 1; i <= num_basis_ - 1; i++)
        {
            center_vec_(i) = i * dis;  
        }
        center_vec_(num_basis_) = 1;

        num_basis_ += 1;
        execute_ = true;
    }

    void PhaseSystem::eval(Eigen::Ref<Eigen::ArrayXd> phase)
    {
        Eigen::ArrayXd Phi = ((z_ - center_vec_.array()).pow(2) / (-2 * width_)).exp();
        phase = Phi / Phi.sum();
    }

    void PhaseSystem::eval_derivative(Eigen::Ref<Eigen::ArrayXd> phase_derivative)
    {
        Eigen::ArrayXd Phi = ((z_ - center_vec_.array()).pow(2) / (-2 * width_)).exp();
        phase_derivative = 
    }

    
}



