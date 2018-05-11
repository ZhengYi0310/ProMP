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
    JS_ProMP::JS_ProMP(int num_basis, double width, double regular_coeff, int num_joints, int traj_timesteps, double traj_dt) : phase_system_(num_basis, width, traj_timesteps), regular_coeff_(regular_coeff), num_joints_(num_joints), traj_timesteps_(traj_timesteps), traj_dt_(traj_dt)
    {
        phase_system_.init();
        this->get_rollout_steps(rollout_steps_);
        this->get_num_basis(num_basis_);
        W_prior_mean_           = Eigen::VectorXd::Zero(num_basis_ * num_joints_);
        W_prior_mean_samples_   = Eigen::MatrixXd::Zero(num_basis_ * num_joints, traj_timesteps_);
        W_prior_covar_          = Eigen::MatrixXd::Zero(num_basis_ * num_joints, num_basis_ * num_joints);
        PHI_                    = Eigen::MatrixXd::Zero(2 * num_joints_ * rollout_steps_, num_joints_ * num_basis_);
    }

    void JS_ProMP::AddDemo(Eigen::MatrixXd demo) 
    {
        assert(demo.rows() == 2 * num_joints_); // need to provide the vlocity vector 
        if (demo.cols()!= traj_timesteps_)
        {
            cout << "[Warning]::The provided demo has longer timestamps than traj_timesteps_, will be cuted!" << std::endl;
            demo = demo.block(0, 0, 2, traj_timesteps_);
        }
        
        demo.resize(demo.size(), 1);
        Y_.push_back(demo);
    }


    void JS_ProMP::BuildDesignMatrix()
    {
        phase_system_.rollout();
        phase_system_.get_rollout(phase_rollout_);
        assert(phase_rollout_.size() == rollout_steps_);
        for (int i = 0; i < rollout_steps_; i++)
        {
            for (int j = 0; j < num_joints_; j++)
            {
                PHI_.block(2 * num_joints_ * i + 2 * j, num_basis_ * j, 2, num_basis_) = phase_rollout_[i].block(0, 0, num_basis_, 2).transpose();
                TAU_.block(num_joints_ * i + j, num_basis_ * j, 1, num_basis_) = phase_rollout_[i].col(2).transpose(); 
            }
        }
    }

    void JS_ProMP::L2Regression()
    {
        // Make sure the demo dimensionality is the same as PHI_ 
        for (int i = 0; i < Y_.size(); i++)
        {
            if (PHI_.rows() != Y_[i].size())
            {
                std::stringstream ss;
                ss << "time steps of PHI_(" << PHI_.rows() << ") is different with that of Y_(" << Y_.size() << ")" << std::endl;
                throw std::runtime_error(ss.str());
            }

            else 
            { 
                W_prior_mean_samples_.col(i) += (PHI_.transpose() * PHI_ + regular_coeff_ * TAU_.transpose() * TAU_).inverse() * PHI_.transpose() * Y_[i];
                W_prior_covar_ += W_prior_mean_samples_ * W_prior_mean_samples_.transpose();
            }
        }

        W_prior_mean_ = W_prior_mean_samples_.colwise().sum() / Y_.size();
        W_prior_covar_ = (Y_.size() * W_prior_covar_ + lambda_W_ * Eigen::MatrixXd::Identity(num_basis_ * num_joints_, num_basis_ * num_joints_)) / (Y_.size() + regular_coeff_);
    }

    void JS_ProMP::AddViaPoints(double t, Eigen::VectorXd y, Eigen::MatrixXd y_covar)
    {
        Via_Points_.via_points_time_ind.push_back(t);
        int ind = Via_Points_.via_points_time_ind.size();
        assert(y.rows() == 2);
        assert(y.cols() == num_joints_);
        assert(y_covar.rows() == num_joints_);
        assert(y_covar.cols() == num_joints_);
        Via_Points_.obs_covar_map[ind] = y_covar;
        Via_Points_.obs_map[ind] = y;
    }

    void JS_ProMP::rollout(double randomness)
    {
        for (int i = 0; i < Via_Points_.via_points_time_ind.size(); i++)
        {
            // conditioning
            Eigen::MatrixXd phi_t =  PHI_.block(2 * num_joints_ * i/* timestamp */, 0, 2 * num_joints_, num_basis_ * num_joints_);

            Eigen::MatrixXd L =  W_prior_covar_conditioned_ * phi_t.transpose() * (Via_Points_.obs_covar_map[i] + phi_t * W_prior_covar_conditioned_ * phi_t.transpose()); 
            W_prior_mean_conditioned_ = W_prior_mean_conditioned_ + L * (Via_Points_.obs_map[i] - phi_t.transpose() * W_prior_mean_conditioned_);
            W_prior_covar_conditioned_ = W_prior_covar_conditioned_ - L * phi_t * W_prior_mean_conditioned_;
        }
    }


}

