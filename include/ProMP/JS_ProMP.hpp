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

namespace ProMP
{
    class JS_ProMP 
    {
        public:
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW

            JS_ProMP(int num_basis=100, double width = 0.05, double regular_coeff, int num_joints, int traj_timesteps);

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
                phase_system_.set_num_basis();
            }

        
            inline void get_num_joints(int& num_joints)
            {
                num_joints = num_joints_;
            }
        
            inline void get_traj_timestamps(double& traj_timestamps)
            {
                traj_timestamps = traj_timestamps;
            }

        
        private:
            PhaseSystem phase_system_;
            Eigen::MatrixXd PHI_;
            Eigen::MatrixXd Y_;
            Eigen::MatrixXd X_;
            Eigen::VectorXd W_prior_mean_;
            Eigen::VectorXd W_prior_covar_;
           
            Eigen::MatrixXd via_points_;
            Eigen::MatrixXd via_points_covar_;
            double regular_coeff_;
            int num_joints_;
            int traj_timesteps_;
            std::shared_pointer<Eigen::EigenMultivariateNormal> MG_;
    }
}
#endif
