/*************************************************************************
	> File Name: JS_ProMP.cpp
	> Author: Yi Zheng 
	> Mail: hczhengcq@gmail.com
	> Created Time: Wed 09 May 2018 07:56:25 PM CEST
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
#include <sstream>
#include <ProMP/JS_ProMP.hpp>
using namespace std;

namespace ProMP
{
    JS_ProMP::JS_ProMP(double regular_coeff, int num_joints, double traj_timesteps, double traj_dt, int num_basis, double width) : phase_system_(traj_timesteps, num_basis, width), regular_coeff_(regular_coeff), num_joints_(num_joints), traj_timesteps_(traj_timesteps), traj_dt_(traj_dt)
    {
        MG_ = Eigen::EigenMultivariateNormal<double>(Eigen::VectorXd::Zero(2 * num_joints_), Eigen::MatrixXd::Identity(2 * num_joints_, 2 * num_joints_)); 
    }

    void JS_ProMP::init()
    {
        phase_system_.init();
        this->get_rollout_steps(rollout_steps_);
        this->get_num_basis(num_basis_);

        Mu_W_      = Eigen::VectorXd::Zero(num_basis_ * num_joints_); 
        Sigma_W_   = Eigen::MatrixXd::Zero(num_basis_ * num_joints_, num_basis_ * num_joints_);
        PHI_       = Eigen::MatrixXd::Zero(2 * num_joints_ * rollout_steps_, num_joints_ * num_basis_);
        PHI_t_     = Eigen::MatrixXd::Zero(2 * num_joints_, num_joints_ * num_basis_);
        Y_rollout_ = Eigen::MatrixXd::Zero(2 * num_joints_, rollout_steps_);
    }

    void JS_ProMP::AddDemo(Eigen::ArrayXXd demo) 
    {
        assert(demo.rows() == 2 * num_joints_); // need to provide the vlocity vector 
        if (demo.cols()!= traj_timesteps_)
        {
            cout << "[Warning]::The provided demo has longer timestamps than traj_timesteps_, will be cuted!" << std::endl;

            //TODO ADD BRACH INTERMS OF TRAINING ON POS AND VEL V.S ONLY ON POS
            demo = demo.block(0, 0, 2 * num_joints_, traj_timesteps_);
            demo = demo.transpose();//.resize(2 * traj_timesteps_, 1);
            demo.resize(2 * traj_timesteps_, 1);
        }
        
        Y_.push_back(demo);
    }


    void JS_ProMP::BuildDesignMatrix()
    {
        phase_system_.rollout();
        phase_system_.get_rollout(phase_rollout_);
        assert(phase_rollout_.size() == rollout_steps_);
        for (int i = 0; i < rollout_steps_; i++)
        {
            //for (int j = 0; j < num_joints_; j++)
            //{
            //    PHI_.block(2 * num_joints_ * i + 2 * j, num_basis_ * j, 2, num_basis_) = phase_rollout_[i].block(0, 0, num_basis_, 2).transpose();
            //    TAU_.block(num_joints_ * i + j, num_basis_ * j, 1, num_basis_) = phase_rollout_[i].col(2).transpose(); 
            //}
            for (int j = 0; j < num_joints_; j++)
            {
                PHI_.block(2 * num_joints_ * i + j, num_basis_ * j, 1, num_basis_) = phase_rollout_[i].col(0).transpose();
                PHI_.block(2 * num_joints_ * i + num_joints_ + j, num_basis_ * j, 1, num_basis_) = phase_rollout_[i].col(1).transpose();
                TAU_.block(num_joints_ * i + j, num_basis_ * j, 1, num_basis_) = phase_rollout_[i].col(2).transpose(); 
            }
        }
    }

    void JS_ProMP::L2Regression()
    {
        // Make sure the Y_ has more than one demos 
        if (Y_.size() <= 1)
        {
            std::stringstream ss;
            ss << Y_.size() << " demos provided! (need more than 1 demo)"<< std::endl;
            throw std::runtime_error(ss.str()) ;
        }

        Eigen::MatrixXd Mu_W_samples(num_basis_ * num_joints_, Y_.size());
       
        for (int i = 0; i < Y_.size(); i++)
        {
             // Make sure the demo dimensionality is the same as PHI_ 
            if (PHI_.rows() != Y_[i].size())
            {
                std::stringstream ss;
                ss << "time steps of PHI_(" << PHI_.rows() << ") is different with that of Y_(" << Y_.size() << ")" << std::endl;
                throw std::runtime_error(ss.str());
            }

            else 
            { 
                Mu_W_samples.col(i) = (PHI_.transpose() * PHI_ + regular_coeff_ * TAU_.transpose() * TAU_).inverse() * PHI_.transpose() * Y_[i].matrix();
               
            }
        }
        
        // compute mean of thw weights 
        Mu_W_ = Mu_W_samples.colwise().sum() / Y_.size();

        // compute covariance of the weights 
        for (int i = 0; i < Y_.size(); i++)
        {
            Sigma_W_ += (Mu_W_samples.col(i) - Mu_W_) * (Mu_W_samples.col(i) - Mu_W_).transpose(); 
        }
        // based on inverse-wishart ditribution for the prior of Sigma_W
        Sigma_W_ = (Y_.size() * Sigma_W_ + lambda_W_ * Eigen::MatrixXd::Identity(num_basis_ * num_joints_, num_basis_ * num_joints_)) / (Y_.size() + regular_coeff_);
    }

    

    void JS_ProMP::rollout()
    {
        // Rollout the learned ProMP 
        for (int i = 0; i < phase_rollout_.size(); i++)
        {
            this->step(i);
        }
        Y_rollout_vec_.push_back(Y_rollout_);
        
    }

    void JS_ProMP::step(int step_ind)
    {
        PHI_t_ = PHI_.block(2 * num_joints_ * step_ind, 0, 2 * num_joints_, num_basis_ * num_joints_); 
        MG_.setMean(PHI_t_ * Mu_W_);
        MG_.setCovar(PHI_t_ * Sigma_W_ * PHI_t_.transpose()); // + prior y;
        Y_rollout_.col(step_ind) = MG_.samples(1);
    }

    void JS_ProMP::reset()
    {
        phase_system_.reset();
        this->get_rollout_steps(rollout_steps_);
        this->get_num_basis(num_basis_);

        Mu_W_      = Eigen::VectorXd::Zero(num_basis_ * num_joints_); 
        Sigma_W_   = Eigen::MatrixXd::Zero(num_basis_ * num_joints_, num_basis_ * num_joints_);
        PHI_       = Eigen::MatrixXd::Zero(2 * num_joints_ * rollout_steps_, num_joints_ * num_basis_);
        PHI_t_     = Eigen::MatrixXd::Zero(2 * num_joints_, num_joints_ * num_basis_);
        Y_rollout_ = Eigen::MatrixXd::Zero(2 * num_joints_, rollout_steps_);

       // BuildDesignMatrix();
    }
 
}

