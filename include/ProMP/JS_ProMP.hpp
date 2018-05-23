/*************************************************************************
	> File Name: JS_ProMP.hpp
	> Author: Yi Zheng 
	> Mail: hczhengcq@gmail.com
	> Created Time: Wed 09 May 2018 07:22:55 PM CEST
 ************************************************************************/

#ifndef _JS_PROMP_H
#define _JS_PROMP_H
#include <memory>
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/LU>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/StdVector>
#include <eigen3/Eigen/Geometry>
#include <ProMP/PhaseSystem.hpp>
#include <ProMP/eigenmvn.h>
#include <map>

namespace ProMP
{
    class JS_ProMP 
    {
        public:
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW

            JS_ProMP(double regular_coeff, int num_joints, double traj_timesteps, double traj_dt, int num_basis=100, double width = 0.05);

            ~JS_ProMP() {}

            inline void get_phase_centers(Eigen::VectorXd& center_vec)
            {
                phase_system_.get_centers(center_vec);
            }
        
            inline void set_phase_centers(Eigen::VectorXd center_vec)
            {
                phase_system_.set(centers_vec);
            }

            inline void get_phase_width_(double& width)
            {
                phase_system_.get_width(width);
            }
        
            inline void set_phase_width(double width)
            {
                phase_system_.set_width(width);
            }
        
            inline void get_num_basis(int& num_basis)
            {
                phase_system_.get_num_basis(num_basis);
            }
        
            inline void set_num_basis(int num_basis)
            {
                phase_system_.set_num_basis(num_basis);
            }

        
            inline void get_num_joints(int& num_joints)
            {
                num_joints = num_joints_;
            }
        
            inline void get_traj_timestamps(double& traj_timestamps)
            {
                traj_timestamps = traj_timestamps_;
            }

            inline void get_rollout_steps(double& rollout_steps)
            {
                phase_system_.get_rollout_steps(rollout_steps);
            }

            inline void get_num_demos(double& num_demos)
            {
                num_demos = Y_.size();
            }

            inline void get_mean(double time_ind, Eigen::Ref<Eigen::VectorXd> mean)
            {
                double step_ind = std::floor(time_ind / traj_dt_);
                get_mean_(step_ind, mean);
            }

            inline void get_mean_step(double step_ind, Eigen::Ref<Eigen::VectorXd> mean)
            {
                //double step_ind = std::floor(time_ind / traj_dt_);
                get_mean_(step_ind, mean);
            }
            
            inline void get_cov(double time_ind, Eigen::Ref<Eigen::MatrixXd> std)
            {
                double step_ind = std::floor(time_ind / traj_dt_);
                get_std_(step_ind, std);
            }

            inline void get_cov_step(double step_ind, Eigen::Ref<Eigen::MatrixXd> std)
            {
                //double step_ind = std::floor(time_ind / traj_dt_);
                get_std_(step_ind, std);
            }

            inline void get_bounds(double time_ind, Eigen::Ref<Eigen::VectorXd> lower_bounds, Eigen::Ref<Eigen::VectorXd> upper_bounds)
            {
                double step_ind = std::floor(time_ind / traj_dt_);
                get_bounds_(step_ind, lower_bounds, upper_bounds);
            }

            inline void get_bounds_step(double step_ind, Eigen::Ref<Eigen::VectorXd> lower_bounds, Eigen::Ref<Eigen::VectorXd> upper_bounds)
            {
                get_bounds_(step_ind, lower_bounds, upper_bounds);
            }

            inline void get_phase(Eigen::Ref<Eigen::MatrixXd> phi, Eigen::Ref<Eigen::MatrixXd> phi_d, Eigen::Ref<Eigen::MatrixXd> phi_dd, double time_ind)
            {
                double step_ind = std::floor(time_ind / traj_dt_);
                get_phase_(Eigen::Ref<Eigen::MatrixXd> phi, Eigen::Ref<Eigen::MatrixXd> phi_d, Eigen::Ref<Eigen::MatrixXd> phi_dd, double step_ind);
            }

            inline void get_phase_step(Eigen::Ref<Eigen::MatrixXd> phi, Eigen::Ref<Eigen::MatrixXd> phi_d, Eigen::Ref<Eigen::MatrixXd> phi_dd, double step_ind)
            {
                //double step_ind = std::floor(time_ind / traj_dt_);
                get_phase_(Eigen::Ref<Eigen::MatrixXd> phi, Eigen::Ref<Eigen::MatrixXd> phi_d, Eigen::Ref<Eigen::MatrixXd> phi_dd, double step_ind);
            }

            inline void AddViaPoints(double step_ind, Eigen::VectorXd y, Eigen::MatrixXd y_covar)
            {
                Via_Points_.push_back(ViaPoints(step_ind, y, y_covar)); 
            }

            inline void set_goal(Eigen::VectorXd goal, Eigen::MatrixXd goal_covar = 1e-6)
            {
                Via_Points_.push_back(ViaPoints(traj_timesteps_, goal, goal_covar));  
            }

            inline void set_start(Eigen::VectorXd start, Eigen::MatrixXd start_covar = 1e-6)
            {
                Via_Points_.push_back(ViaPoints(0, start, start_covar));   
            }


            inline void temporal_scaling(double scale)
            {
                assert(scale > 0);
                phase_system_.temporal_scaling(scale);
            }

            
            inline void clear_via_points()
            {
                Via_Points_.clear();
            }

            inline void clear_demos()
            {
                Y_.clear();
            }

            inline void clear_rollouts()
            {
                Y_rollout_vec_.clear();
            }

            /** Add a demonstration trajectory with timestamps T and joints number N to the system 
             * /param demo, a Matrix represents the trajectory (N * 2, T)
             */
            void AddDemo(Eigen::ArrayXXd demo);
            void L2Regression();

            /* Build the Phasis value design matrix PHI 
             */
            void BuildDesignMatrix();

            void AddViaPoints(double t, Eigen::VectorXd y, Eigen::MatrixXd y_covar);
            void rollout();
            void step(int step_ind);
            void reset()

            

        
        private:
            void init();
            PhaseSystem phase_system_;
            Eigen::MatrixXd PHI_;
            Eigen::MatrixXd PHI_t_;
            Eigen::MatrixXd TAU_; // used to store the jerk
            MatrixVector Y_; // stl container for training examples, each example should contain 2 * traj_timesteps_ points. 
            MatrixVector phase_rollout_;
           
            Eigen::VectorXd Mu_W_;   
            Eigen::MatrixXd Sigma_W_;
            Eigen::VectorXd Mu_W_cond_;
            Eigen::MatrixXd Sigma_W_cond_;

            Eigen::MatrixXd Y_rollout_;
            std::vector<Eigen::MatrixXd, Eigen::aligned_allocator<Eigen::MatrixXd> > Y_rollout_vec_;

            struct ViaPoints
            {
                /*
                std::vector<double> via_points_time_ind;
                std::map<int, Eigen::MatrixXd, std::less<int>, Eigen::aligned_allocator<std::pair<const int, Eigen::MatrixXd> > > obs_covar_map;
                 std::map<int, Eigen::VectorXd, std::less<int>, Eigen::aligned_allocator<std::pair<const int, Eigen::VectorXd> > > obs_map;
                */
                double step_ind_;

                Eigen::MatrixXd y_;
                Eigen::MatrixXd y_covar_;

                ViaPoints(double step_ind, Eigen::MatrixXd y, Eigen::MatrixXd y_covar) : step_ind_(step_ind), y_(y), y_covar_(y_covar)
            };
            
            std::vector<ViaPoints> Via_Points_;

            double regular_coeff_;
            double lambda_W_;
            int num_joints_;
            int num_basis_;
            double traj_dt_;
            double traj_timesteps_;
            double rollout_steps_;
            //ViaPoints Via_Points_;
            Eigen::EigenMultivariateNormal<double> MG_;

            /** Get the mean value for the specific time_steps
             *
             */
            inline void get_mean_(double step_ind, Eigen::Ref<Eigen::VectorXd> mean)
            {
                assert(0 <= step_ind && step_ind <= phase_rollout_.size());
                Eigen::MatrixXd phi = PHI_.block(2 * num_joints_ * time_ind, 0, 2 * num_joints_, num_basis_ * num_joints_); 
                mean =  phi * Mu_W_;
            }

            inline void get_cov_(double step_ind, Eigen::Ref<Eigen::MatrixXd> cov)
            {
                assert(0 <= step_ind && step_ind <= phase_rollout_.size());
                Eigen::MatrixXd phi = PHI_.block(2 * num_joints_ * time_ind, 0, 2 * num_joints_, num_basis_ * num_joints_);  
                cov = (phi * Sigma_W_ * phi.transpose()); //+ covairance of the demonstrations 
            }

            inline void get_phase_(Eigen::Ref<Eigen::MatrixXd> phi, Eigen::Ref<Eigen::MatrixXd> phi_d, Eigen::Ref<Eigen::MatrixXd> phi_dd, double step_ind)
            {
                assert(0 <= step_ind && step_ind <= phase_rollout_.size());
                phi = phase_rollout_[step_ind].col(0);
                phi = phase_rollout_[step_ind].col(1);
                phi = phase_rollout_[step_ind].col(2);
            }

            inline void get_bounds_(double step_ind, Eigen::Ref<Eigen::VectorXd> lower_bounds, Eigen::Ref<Eigen::VectorXd> upper_bounds)
            {
                Eigen::VectorXd mean;
                Eigen::MatrixXd cov;

                get_mean_(step_ind, mean);
                get_cov_(step_ind, cov);

                Eigen::ArrayXd std = cov.diagonal().array().sqrt();

                lower_bounds = mean - std.matrix() * 2;
                upper_bounds = mean + std.matrix() * 2;
            }

            //inline void get_bounds(int time_ind, )
    }
}
#endif
