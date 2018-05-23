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
    PhaseSystem::PhaseSystem(double traj_timesteps, double num_basis, double width) : num_basis_(num_basis), width_(width), traj_timesteps_(traj_timesteps), z_(0) 
    {
        z_dot_ = 1;
        execute_ = false;
        rollout_steps_ = std::floor(traj_timesteps_ / z_dot_);
        center_vec_.resize(num_basis_ + 1);
        
    }
    
    /**Generate the vector for the phase system
     */
    void PhaseSystem::init()
    {
        double dis = 1 / num_basis_;
        center_vec_(0) = 0.0;
        for (int i = 1; i <= num_basis_ - 1; i++)
        {
            center_vec_(i) = i * dis;  
        }
        num_basis_ += 1;
        center_vec_(num_basis_ - 1) = 1;
     
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
        //cout << phase(0) << endl;
        //phase /= phase.sum();
    }

    // evaluate the derivative of the phase 
    void PhaseSystem::eval_d(Eigen::Ref<Eigen::ArrayXd> phase_dot, const Eigen::Ref<const Eigen::ArrayXd> phase)
    {
        phase_dot = phase * (z_ - center_vec_.array()) / width_ * -1;
        //phase_dot /=  phase_dot.sum();
        //phase_dot /= phase.sum();
    }

    // evaluate the third derivative (jerk) of the phase 
    void PhaseSystem::eval_ddd(Eigen::Ref<Eigen::ArrayXd> phase_jerk, const Eigen::Ref<const Eigen::ArrayXd> phase)
    {
        phase_jerk = phase * (((z_ - center_vec_.array()) / width_).pow(2) - 1 / width_);
        //phase_jerk /= phase_jerk.sum();
        //phase_jerk /= phase.sum();
    }

    void PhaseSystem::step(Eigen::Ref<Eigen::ArrayXd> phase,
                           Eigen::Ref<Eigen::ArrayXd> phase_dot,
                           Eigen::Ref<Eigen::ArrayXd> phase_jerk)
    {
        //cout << z_ << endl;
        if (z_ <= (1 + 1e-10))
        {
            z_ = z_ + 1.0 / rollout_steps_;
        }
        else 
            z_ = z_;


        eval(phase);
        eval_d(phase_dot, phase);
        eval_ddd(phase_jerk, phase);
        
       
        // Normalize everything
        double p = phase.sum();
        double pd = phase_dot.sum();
        double pdd = phase_jerk.sum();

        phase_jerk = phase_jerk / p - 2 * phase_dot * pd / std::pow(p, 2) - phase * pdd / std::pow(p, 2) + 2 * phase * std::pow(pd, 2) / std::pow(p, 3);  
        phase_dot = (phase_dot * phase.sum() - phase * phase_dot.sum()) / std::pow(p, 2);
        phase = phase / phase.sum();

        // Multiply with z_, z_dot_ 
        phase_dot = phase_dot * z_dot_;
        phase_jerk = phase_jerk * (pow(z_dot_, 2));
        z_vecs_.push_back(z_);    
    }
    
    void PhaseSystem::reset()
    {
        rollout_steps_ = std::floor(traj_timesteps_ / z_dot_);
 
        z_ = 0;
        //std::cout << rollout_steps_ << std::endl;
        //center_vec_.resize(num_basis_ + 1);
       
        z_vecs_.clear();
        rollout_.clear();

        phase_prealloc_ = Eigen::ArrayXd::Zero(num_basis_);
        phase_dot_prealloc_ = Eigen::ArrayXd::Zero(num_basis_);
        phase_jerk_prealloc_ = Eigen::ArrayXd::Zero(num_basis_);
        phase_terms_ = Eigen::ArrayXXd::Zero(num_basis_, 3);
        execute_ = true;
    }

    void PhaseSystem::rollout()
    {
        //std::cout << rollout_steps_ << std::endl;

        if (execute_ == true)
        {
            while (z_ <= (1.0 + 1e-10))
            {
                //cout << z_ << endl; 
                step(phase_prealloc_, phase_dot_prealloc_, phase_jerk_prealloc_);
                phase_terms_.col(0) = phase_prealloc_;
                phase_terms_.col(1) = phase_dot_prealloc_;
                phase_terms_.col(2) = phase_jerk_prealloc_; 
                rollout_.push_back(phase_terms_);
            }
            rollout_steps_ += 1;
        }
        else 
            throw std::runtime_error("The phase system should not be rolling out right now!"); 

        //std::cout << rollout_.size() << " " << rollout_[1].cols() << " " << rollout_[1].rows() << std::endl;

        execute_ = false;
        //reset();
    }   
}



