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
    PhaseSystem::PhaseSystem(int num_basis, double width, int traj_timesteps) : num_basis_(num_basis), width_(width), traj_timesteps_(traj_timesteps), z_(0) 
    {
        z_dot_ = 1;
        execute_ = false;
        rollout_steps_ = std::floor(traj_timesteps_ / z_dot_);
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

        phase_prealloc_ = Eigen::ArrayXd::Zero(num_basis_);
        phase_dot_prealloc_ = Eigen::ArrayXd::Zero(num_basis_);
        phase_jerk_prealloc_ = Eigen::ArrayXd::Zero(num_basis_);
        phase_terms_ = Eigen::ArrayXXd::Zero(num_basis_, 3);
    }

    // evaluate the phase value 
    void PhaseSystem::eval(Eigen::Ref<Eigen::ArrayXd> phase)
    {
        phase = ((z_ - center_vec_.array()).pow(2) / (-2 * width_)).exp();
        phase /= phase.sum();
    }

    // evaluate the derivative of the phase 
    void PhaseSystem::eval_d(Eigen::Ref<Eigen::ArrayXd> phase_dot, const Eigen::Ref<const Eigen::ArrayXd> phase)
    {
        phase_dot = phase * (z_ - center_vec_.array()) / width_ * -1;
        phase_dot /=  phase_dot.sum();
    }

    // evaluate the third derivative (jerk) of the phase 
    void PhaseSystem::eval_ddd(Eigen::Ref<Eigen::ArrayXd> phase_jerk, const Eigen::Ref<const Eigen::ArrayXd> phase)
    {
        phase_jerk = phase * ((z_ - center_vec_.array() / width_).pow(2) - 1 / width_);
        phase_jerk /= phase_jerk.sum();
    }

    void PhaseSystem::step(Eigen::Ref<Eigen::ArrayXd> phase,
                           Eigen::Ref<Eigen::ArrayXd> phase_dot,
                           Eigen::Ref<Eigen::ArrayXd> phase_jerk)
    {
        if (!execute_)
            throw std::runtime_error("The phase system should not be execuing a step right now!");
        eval(phase);
        eval_d(phase_dot, phase);
        eval_ddd(phase_jerk, phase);

        z_ = z_ + z_dot_;
    }
    
    void PhaseSystem::reset()
    {
        z_dot_ = 1;
        execute_ = false;
        rollout_steps_ = std::floor(traj_timesteps_ / z_dot_);
        init();
    }

    void PhaseSystem::rollout()
    {
        if (execute_ == true)
        {
            while (z_ <= 1)
            {
                
                step(phase_prealloc_, phase_dot_prealloc_, phase_jerk_prealloc_);
                phase_terms_.col(0) = phase_prealloc_;
                phase_terms_.col(1) = phase_dot_prealloc_;
                phase_terms_.col(2) = phase_jerk_prealloc_;
                rollout_.push_back(Eigen::Map<Eigen::MatrixXd>(phase_terms_.data(), 3, num_basis_));
            }
        }
        else 
            throw std::runtime_error("The phase system should not be rolling out right now!"); 

        execute_ = false;
        reset();
    }    
}



