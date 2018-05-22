/*************************************************************************
	> File Name: Phase_test.cpp
	> Author: Yi Zheng 
	> Mail: hczhengcq@gmail.com
	> Created Time: Mon 21 May 2018 05:21:34 PM CEST
 ************************************************************************/

#include <iostream>
#include <string>
#include <vector>
#include <set>
#include <map>
#include <algorithm>
#include <cmath>
#include <cassert>
#include <ProMP/PhaseSystem.hpp>
#include "catch/catch.hpp"
#include "matplotlib-cpp/matplotlibcpp.h"
using namespace std;
using namespace ProMP;
using namespace matplotlibcpp;

TEST_CASE("Phase system can be successfully initialized.", "[Phase_System]")
{
    PhaseSystem test_phase_sys1(100, 100, 0.02);
    SECTION("Phase system construction case 1")
    {
        
        REQUIRE(!test_phase_sys1.can_execute());
    }

    SECTION("Phase system initializetion case 2")
    {
        test_phase_sys1.init();
    
        Eigen::VectorXd centers_vec;
        double roll_out_steps;
    
        test_phase_sys1.get_centers(centers_vec);
        test_phase_sys1.get_rollout_steps(roll_out_steps);
    
        REQUIRE(roll_out_steps == 100);
        REQUIRE(centers_vec.rows() == 101);
        REQUIRE(centers_vec(0) == 0.0);
        REQUIRE(centers_vec(centers_vec.rows() - 1) == 1);
        REQUIRE(test_phase_sys1.can_execute());
    }

    PhaseSystem test_phase_sys2(300, 500, 0.02);
    SECTION("Phase system construction case 2")
    {
        
        REQUIRE(!test_phase_sys2.can_execute());
    }

    SECTION("Phase system initializetion case 2")
    {
        test_phase_sys2.init();
    
        Eigen::VectorXd centers_vec;
        double roll_out_steps;
    
        test_phase_sys2.get_centers(centers_vec);
        test_phase_sys2.get_rollout_steps(roll_out_steps);
    
        REQUIRE(roll_out_steps == 300);
        REQUIRE(centers_vec.rows() == 501);
        REQUIRE(centers_vec(0) == 0.0);
        REQUIRE(centers_vec(centers_vec.rows() - 1) == 1);
        REQUIRE(test_phase_sys2.can_execute());
    }
}

TEST_CASE("The centers of the phase system are placed correctly.", "[Phase_System]")
{
    SECTION("case 1")
    {
        double num_basis = 100;
        PhaseSystem test_phase_sys(100, num_basis, 0.02);
        test_phase_sys.init();


        double dis = 1 / num_basis;
        Eigen::VectorXd test_vec, centers_vec;
        test_vec.resize(num_basis + 1);

        test_vec(0) = 0.0;
        for (int i = 1; i <= num_basis - 1; i++)
        {
            test_vec(i) = i * dis;  
        }
        num_basis += 1;
        test_vec(num_basis - 1) = 1;
        
        test_phase_sys.get_centers(centers_vec);

        REQUIRE(test_vec.size() == centers_vec.size());
        
        
        for (int i = 0; i < test_vec.size(); i++)
        {
            REQUIRE(test_vec(i) == centers_vec(i));
            REQUIRE(test_vec(i) == i * dis);
        }
        
        
    }

    SECTION("case 1")
    {
        double num_basis = 30;
        PhaseSystem test_phase_sys(100, num_basis, 0.02);
        test_phase_sys.init();


        double dis = 1 / num_basis;
        Eigen::VectorXd test_vec, centers_vec;
        test_vec.resize(num_basis + 1);

        test_vec(0) = 0.0;
        for (int i = 1; i <= num_basis - 1; i++)
        {
            test_vec(i) = i * dis;  
        }
        num_basis += 1;
        test_vec(num_basis - 1) = 1;
        
        test_phase_sys.get_centers(centers_vec);

        REQUIRE(test_vec.size() == centers_vec.size());
        
        
        for (int i = 0; i < test_vec.size(); i++)
        {
            REQUIRE(test_vec(i) == centers_vec(i));
            REQUIRE(test_vec(i) == Approx(i * dis));
            //std::cout << test_vec(i) << std::endl;
        }
        
        
    }
}

TEST_CASE("Phase system rolls out.", "[Phase_System]")
{
    SECTION("case 1")
    {
        double num_basis = 100;
        PhaseSystem test_phase_sys(100, num_basis, 0.02);
        test_phase_sys.init();
        std::vector<double> z_vec;
        
        test_phase_sys.rollout();

        test_phase_sys.get_phase_values(z_vec);
        REQUIRE(z_vec.size() == 100);
        
        
        

        test_phase_sys.temporal_scaling(0.5);
        test_phase_sys.reset();
        test_phase_sys.rollout();
        test_phase_sys.get_phase_values(z_vec);
        REQUIRE(z_vec.size() == 200);

        test_phase_sys.temporal_scaling(8);
        test_phase_sys.reset();
        test_phase_sys.rollout();
        test_phase_sys.get_phase_values(z_vec);
        REQUIRE(z_vec.size() == 25);

        test_phase_sys.temporal_scaling(0.25);
        test_phase_sys.reset();
        test_phase_sys.rollout();
        test_phase_sys.get_phase_values(z_vec);
        REQUIRE(z_vec.size() == 100);
    }

    SECTION("case 2")
    {
        double num_basis = 50;
        PhaseSystem test_phase_sys(500, num_basis, 0.02);
        test_phase_sys.init();
        test_phase_sys.rollout();
        std::vector<double> z_vec;
        test_phase_sys.get_phase_values(z_vec);
        REQUIRE(z_vec.size() == 500);
        MatrixVector rollout;
        test_phase_sys.get_rollout(rollout);

        test_phase_sys.temporal_scaling(0.3);
        test_phase_sys.reset();
        test_phase_sys.rollout();
        test_phase_sys.get_phase_values(z_vec);
        REQUIRE(z_vec.size() == 1666);

        test_phase_sys.temporal_scaling(10);
        test_phase_sys.reset();
        test_phase_sys.rollout();
        test_phase_sys.get_phase_values(z_vec);
        REQUIRE(z_vec.size() == 167);

        test_phase_sys.temporal_scaling(1.0 / 3.0);
        test_phase_sys.reset();
        test_phase_sys.rollout();
        test_phase_sys.get_phase_values(z_vec);
        REQUIRE(z_vec.size() == 500);
    } 
}


TEST_CASE("results plotting", "[Phase_system]")
{
    PhaseSystem test_phase_sys(100, 30, 0.5);
    test_phase_sys.init();
    std::vector<double> z_vec;
        
    test_phase_sys.rollout();

    test_phase_sys.get_phase_values(z_vec);
    REQUIRE(z_vec.size() == 100);
    
    MatrixVector rollout;
    Eigen::VectorXd centers_vec;
    double rollout_steps; 
    int num_basis;
    test_phase_sys.get_rollout_steps(rollout_steps);
    test_phase_sys.get_num_basis(num_basis);
    test_phase_sys.get_rollout(rollout);
    test_phase_sys.get_centers(centers_vec);

    SECTION("rollout dimension is correct")
    {
        REQUIRE(rollout.size() == rollout_steps);
        REQUIRE(num_basis == 31);
        for (int i = 0; i < rollout.size(); i++)
        {
            REQUIRE(rollout[i].rows() == num_basis);
        }
    }

    SECTION("plot out the phase value")
    {
        Eigen::MatrixXd Basis_mat;
        std::vector<std::vector<double> > Basis_vec;
        Basis_mat.resize(num_basis, rollout_steps);
        /*
        for (int i = 0; i < centers_vec.rows(); i++)
        {
            cout << centers_vec.row(i) << std::endl;
        }
        */

        for (int i = 0; i < rollout.size(); i++)
        {
            Basis_mat.col(i) = rollout[i].col(2);  
        }

        for (int i = 0; i < num_basis; i++)
        {   
            std::vector<double> temp;
            for (int j = 0; j < Basis_mat.cols(); j++)
            {
                temp.push_back(Basis_mat(i, j));
            } 
            Basis_vec.push_back(temp);
            matplotlibcpp::plot(z_vec, Basis_vec[i]);
        }
        matplotlibcpp::show();
    }

    SECTION("after rescaling")
    {
        test_phase_sys.temporal_scaling(0.3);
        test_phase_sys.reset();
        test_phase_sys.rollout();
        test_phase_sys.get_rollout_steps(rollout_steps);
        test_phase_sys.get_num_basis(num_basis);
        test_phase_sys.get_rollout(rollout);
        test_phase_sys.get_centers(centers_vec);
        test_phase_sys.get_phase_values(z_vec);

        cout << z_vec.size() << endl;
        REQUIRE(z_vec.size() == rollout_steps);
        REQUIRE(rollout.size() == rollout_steps);
        REQUIRE(num_basis == 31);
        for (int i = 0; i < rollout.size(); i++)
        {
            REQUIRE(rollout[i].rows() == num_basis);
        }

        Eigen::MatrixXd Basis_mat;
        std::vector<std::vector<double> > Basis_vec;
        Basis_mat.resize(num_basis, rollout_steps);

        for (int i = 0; i < rollout.size(); i++)
        {
            Basis_mat.col(i) = rollout[i].col(1);  
        }

        for (int i = 0; i < num_basis; i++)
        {   
            std::vector<double> temp;
            for (int j = 0; j < Basis_mat.cols(); j++)
            {
                temp.push_back(Basis_mat(i, j));
            } 
            Basis_vec.push_back(temp);
            matplotlibcpp::plot(z_vec, Basis_vec[i]);
        }
        matplotlibcpp::show();
    }

    Eigen::MatrixXd test;
    cout << test.rows() << " " << test.cols() << endl;
}





