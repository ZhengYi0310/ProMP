/*************************************************************************
	> File Name: PhaseSystem.hpp
	> Author: Yi Zheng 
	> Mail: hczhengcq@gmail.com
	> Created Time: Wed 09 May 2018 01:08:17 PM CEST
 ************************************************************************/

#ifndef _PHASESYSTEM_H
#define _PHASESYSTEM_H
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/LU>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/StdVector>
#include <eigen3/Eigen/Geometry>

typedef std::vector<Eigen::MatrixXd, Eigen::aligned_allocator<Eigen::MatrixXd> > MatrixVector;

namespace ProMP
{
    /**
     * the phase system for the probabilistic movement primitive
     */
    class PhaseSystem 
    {
        public:
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW

            PhaseSystem(int num_basis=100, double width = 0.05, const double goal_t, const double start_t); // : num_basis_(num_basis), overlap_(overlap), z_(0.0)

            ~PhaseSystem() {}

            void init()
            void eval(Eigen::Ref<Eigen::ArrayXd> phase);
            void eval_derivative(Eigen::Ref<Eigen::MatrixXd> phase_dot);
            void eval_jerk(Eigen::Ref<Eigen::MatrixXd> phase_jerk);

            void step(Eigen::Ref<Eigen::MatrixXd> phase, 
                      Eigen::Ref<Eigen::MatrixXd> phase_dot, 
                      Eigen::Ref<Eigen::MatrixXd> phase_jerk);

            void rollout()
            
            inline void get_centers(Eigen::VectorXd& centers_vec)
            {
                center_vec = center_vec_;
            }

            inline void get_width(double& width)
            {
                width = width_;
            }

            inline void set_centers(Eigen::VectorXd center_vec)
            {
                // Assert size of the two vecs
                assert(center_vec_.size() == center_vec.size());
                center_vec_ = center_vec;
            }

            inline void set_width(double width)
            {
                width_ = width;
            }

            inline void temoral_scaling(int scale)
            {
                assert(scale > 0);
                z_dot_ = z_dot_ * scale;
            }

        private:
            bool execute_;

            int num_basis_;
            double overlap_;

            double traj_duration_;
            double z_;
            double z_dot_;

            Eigen::VectorXd center_vec_;
            
            double width_;     

            MatrixVector rollout_;
    };
}
#endif
