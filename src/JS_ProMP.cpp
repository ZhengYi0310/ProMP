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
#include <ProMP/JS_ProMP.hpp>
using namespace std;

namespace ProMP
{
    JS_ProMP::JS_ProMP(int num_basis, double width, double regular_coeff, int num_joints, int traj_timesteps) : phase_system_(num_basis, width, traj_timesteps), regular_coeff_(regular_coeff), num_joints_(num_joints)
    {
        
    }
}

