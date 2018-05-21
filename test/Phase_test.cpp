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
using namespace std;
using namespace ProMP;

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
        test_phase_sys.rollout();
    }

    SECTION("case 2")
    {
        double num_basis = 50;
        PhaseSystem test_phase_sys(300, num_basis, 0.02);
        test_phase_sys.init();
        test_phase_sys.rollout();
    }

    SECTION("case 2")
    {
        double num_basis = 30;
        PhaseSystem test_phase_sys(320, num_basis, 0.02);
        test_phase_sys.init();
        test_phase_sys.rollout();
    }

    SECTION("case 3")
    {
        double num_basis = 32;
        PhaseSystem test_phase_sys(320, num_basis, 0.02);
        test_phase_sys.init();
        test_phase_sys.rollout();
    }
}





