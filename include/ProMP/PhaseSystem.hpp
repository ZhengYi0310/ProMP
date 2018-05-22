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

typedef std::vector<Eigen::ArrayXXd, Eigen::aligned_allocator<Eigen::ArrayXXd> > MatrixVector;

namespace ProMP
{
    /**
     * the phase system for the probabilistic movement primitive
     */
    class PhaseSystem 
    {
        public:
            EIGEN_MAKE_ALIGNED_OPERATOR_NEW

            PhaseSystem(double traj_timesteps, double num_basis=100, double width = 0.05); // : num_basis_(num_basis), overlap_(overlap), z_(0.0)

            ~PhaseSystem() {}

            void init();
            void eval(Eigen::Ref<Eigen::ArrayXd> phase);
            void eval_d(Eigen::Ref<Eigen::ArrayXd> phase_dot, const Eigen::Ref<const Eigen::ArrayXd> phase);
            void eval_ddd(Eigen::Ref<Eigen::ArrayXd> phase_jerk, const Eigen::Ref<const Eigen::ArrayXd> phase);

            void step(Eigen::Ref<Eigen::ArrayXd> phase, 
                      Eigen::Ref<Eigen::ArrayXd> phase_dot, 
                      Eigen::Ref<Eigen::ArrayXd> phase_jerk);

            void rollout();
            inline void get_rollout(MatrixVector& rollout) {rollout = rollout_;}
            void reset();
            
            
            inline void get_centers(Eigen::VectorXd& center_vec)
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
                assert(width > 0);
                width_ = width;
            }

            inline void set_num_basis(int num_basis)
            {
                assert(num_basis > 1);
                num_basis_ = num_basis;
            }

            inline void get_num_basis(int& num_basis)
            {
                num_basis = num_basis_;
            }

            inline void temporal_scaling(double scale)
            {
                assert(scale > 0);
                z_dot_ = z_dot_ * scale;
            }

            inline void get_rollout_steps(double& rollout_steps)
            {
                rollout_steps = rollout_steps_;
            }

            inline bool can_execute()
            {
                return execute_;
            }

            inline void get_phase_values(std::vector<double>& z_vecs)
            {
                z_vecs = z_vecs_;
            }

        private:
            bool execute_;
            double num_basis_;
            double traj_timesteps_;
            double rollout_steps_;
            double z_;
            double z_dot_;

            Eigen::VectorXd center_vec_;

            Eigen::ArrayXd phase_prealloc_;
            Eigen::ArrayXd phase_dot_prealloc_;
            Eigen::ArrayXd phase_jerk_prealloc_;
            Eigen::ArrayXXd phase_terms_;
            
            double width_;     

            MatrixVector rollout_;
            std::vector<double> z_vecs_;
    };
}
#endif
