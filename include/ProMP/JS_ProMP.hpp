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

            JS_ProMP(int num_basis=100, double width = 0.05, double regular_coeff, int num_joints, int traj_timesteps, double traj_dt);

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
                double step_ind = time_ind / traj_dt_;
                get_mean_(step_ind, mean);
            }
            
            inline void get_std(double time_ind, Eigen::Ref<Eigen::MatrixXd> std)
            {
                double step_ind = time_ind / traj_dt_;
                get_std_(step_ind, std);
            }
            /** Add a demonstration trajectory with timestamps T and joints number N to the system 
             * /param demo, a Matrix represents the trajectory (N * 2, T)
             */
            void AddDemo(Eigen::MatrixXd demo);
            void L2Regression();
            void BuildDesignMatrix();
            void AddViaPoints(double t, Eigen::VectorXd y, Eigen::MatrixXd y_covar);

        
        private:
            PhaseSystem phase_system_;
            Eigen::MatrixXd PHI_;
            Eigen::MatrixXd TAU_; // used to store the jerk
            MatrixVector Y_; // stl container for training examples, each example should contain 2 * traj_timesteps_ points. 
            MatrixVector phase_rollout_;
            Eigen::MatrixXd X_;
            Eigen::VectorXd W_prior_mean_;
            Eigen::MatrixXd W_prior_mean_samples_;
            Eigen::VectorXd W_prior_covar_;

            struct ViaPoints
            {
                std::vector<double> via_points_time_ind;
                std::map<int, Eigen::MatrixXd, std::less<int>, Eigen::aligned_allocator<std::pair<const int, Eigen::MatrixXd> > > obs_covar_map;
                 std::map<int, Eigen::VectorXd, std::less<int>, Eigen::aligned_allocator<std::pair<const int, Eigen::VectorXd> > > obs_map;
            };
            
           
            Eigen::MatrixXd via_points_;
            Eigen::MatrixXd via_points_covar_;
            double regular_coeff_;
            double lambda_W_;
            int num_joints_;
            int num_basis_;
            double traj_dt_;
            double traj_timesteps_;
            double rollout_steps_;
            ViaPoints Via_Points_;
            std::shared_pointer<Eigen::EigenMultivariateNormal> MG_;

            /** Get the mean value for the specific time_steps
             *
             */
            inline void get_mean_(double step_ind, Eigen::Ref<Eigen::VectorXd> mean)
            {
                assert(0 <= time_ind && time_ind <= phase_rollout_.size());
                Eigen::MatrixXd phi = PHI_.block(2 * num_joints_ * time_ind, 0, 2 * num_joints_, num_basis_ * num_joints_); 
                mean =  phi * W_prior_mean_;
            }

            inline void get_std_(double step_ind, Eigen::Ref<Eigen::MatrixXd> std)
            {
                assert(0 <= time_ind && time_ind <= phase_rollout_.size());
                Eigen::MatrixXd phi = PHI_.block(2 * num_joints_ * time_ind, 0, 2 * num_joints_, num_basis_ * num_joints_);  
                std = (phi * W_prior_covar_ * phi.transpose()).sqrt(); //+ covairance of the demonstrations 
            }

            //inline void get_bounds(int time_ind, )
    }
}
#endif
