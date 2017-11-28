#pragma once
#ifndef _GLOBAL_VARS_H_
#define _GLOBAL_VARS_H_

#include <cstdlib>
#include <cstdio>
#include <deal.II/base/parameter_handler.h>
#include "nr.h"

using namespace dealii;
namespace Parameters{

struct Geometry{
    unsigned int xstep;
    unsigned int ystep;
    std::string wound_shape;
    unsigned int wound_radius;
    unsigned int wound_a;
    unsigned int wound_b;
    unsigned int wound_side_length;

    static void
    declare_parameters(ParameterHandler &prm);
 
    void
    parse_parameters(ParameterHandler &prm);
};

struct Time{
    double time_total;
    double time_step;

    static void
    declare_parameters(ParameterHandler &prm);

    void
    parse_parameters(ParameterHandler &prm);
};

struct Fibroblast{
    unsigned int FNinit;
    double rho1;
    double rho2;
    
    static void
    declare_parameters(ParameterHandler &prm);

    void
    parse_parameters(ParameterHandler &prm);
};

struct ECM{
    double pc;
    double dc;
    double df;
    double kappa;
    unsigned int ECM_gridsize;
 
    static void
    declare_parameters(ParameterHandler &prm);

    void
    parse_parameters(ParameterHandler &prm);
};

struct Endocyte{
    unsigned int ENinit;
    
    static void
    declare_parameters(ParameterHandler &prm);

    void
    parse_parameters(ParameterHandler &prm);
};

struct AllParameters : public Geometry,
                       public Time,
                       public Fibroblast,
                       public ECM,
                       public Endocyte{
    AllParameters(const std::string &input_file);

    static void
    declare_parameters(ParameterHandler &prm);

    void
    parse_parameters(ParameterHandler &prm);
};
}

#endif
