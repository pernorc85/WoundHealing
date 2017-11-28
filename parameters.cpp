#include "parameters.h"
#include <cstdlib>
#include <cstdio>

using namespace dealii;
Mat3D_DP *fibroblast_density_ptr = nullptr;
Mat3D_DP *speedfield_ptr = nullptr;
Mat3D_DP *thetax_ptr = nullptr;
Mat3D_DP *thetay_ptr = nullptr;
Mat3D_DP *thetaz_ptr = nullptr;
Mat3D_DP *tissue_displacement_x_ptr = nullptr;
Mat3D_DP *tissue_displacement_y_ptr = nullptr;
Mat3D_DP *tissue_displacement_z_ptr = nullptr;

namespace Parameters{

void Geometry::declare_parameters(ParameterHandler &prm){
    prm.enter_subsection("Geometry");
    {
        prm.declare_entry("xstep", "700",
                          Patterns::Integer(0),
                          "domain dimension in x-direction");
        prm.declare_entry("ystep", "700",
                          Patterns::Integer(0),
                          "domain dimension in y-direction");
        prm.declare_entry("wound_shape", "round",
                          Patterns::Selection("round|elliptic|disk"),
                          "Type of wound shape");
        prm.declare_entry("wound_radius", "250",
                          Patterns::Integer(0),
                          "wound radius, if wound is disk or half-sphere");
        prm.declare_entry("wound_a", "200",
                          Patterns::Integer(0),
                          "wound long radius, if wound is oval");
        prm.declare_entry("wound_b", "150",
                          Patterns::Integer(0),
                          "wound short radius, if wound is oval");
        prm.declare_entry("wound_side_length", "350",
                          Patterns::Integer(0),
                          "wound short radius, if wound is square");
    }
    prm.leave_subsection();
}

void Geometry::parse_parameters(ParameterHandler &prm){
    prm.enter_subsection("Geometry");
    {
        xstep = prm.get_integer("xstep");
        ystep = prm.get_integer("ystep");
        wound_shape = prm.get("wound_shape");
        wound_radius = prm.get_integer("wound_radius");
        wound_a = prm.get_integer("wound_a");
        wound_b = prm.get_integer("wound_b");
        wound_side_length = prm.get_integer("wound_side_length");
    }
    prm.leave_subsection();
}

void Time::declare_parameters(ParameterHandler &prm){
    prm.enter_subsection("Time");
    {
        prm.declare_entry("time_total", "400.0",
                          Patterns::Double(0.0),
                          "total simulated time in hr");
        prm.declare_entry("time_step", "0.15",
                          Patterns::Double(0.00),
                          "time step length in hr");
    }
    prm.leave_subsection();
}

void Time::parse_parameters(ParameterHandler &prm){
    prm.enter_subsection("Time");
    {
        time_total = prm.get_double("time_total");
        time_step = prm.get_double("time_step");
    }
    prm.leave_subsection();
}

void Fibroblast::declare_parameters(ParameterHandler &prm){
    prm.enter_subsection("Fibroblast");
    {
        prm.declare_entry("FNinit", "1200",
                          Patterns::Integer(0),
                          "Number of FNinit in healthy tissue per mm^3");
        prm.declare_entry("rho1", "0.3",
                          Patterns::Double(0.0),
                          "the influence of chemokine gradient on direction of fibroblast");
        prm.declare_entry("rho2", "0.3",
                          Patterns::Double(0.0),
                          "the influence of collagen orientation on direction of fibroblast");
    }
    prm.leave_subsection();
}

void Fibroblast::parse_parameters(ParameterHandler &prm){
    prm.enter_subsection("Fibroblast");
    {
        FNinit = prm.get_integer("FNinit");
        rho1 = prm.get_double("rho1");
        rho2 = prm.get_double("rho2");
    }
    prm.leave_subsection();
}

void ECM::declare_parameters(ParameterHandler &prm){
    prm.enter_subsection("ECM");
    {
        prm.declare_entry("pc", "0.44",
                          Patterns::Double(0.0),
                          "Rate of collagen production");
        prm.declare_entry("dc", "0.44",
                          Patterns::Double(0.0),
                          "Rate of collagen degradation");
        prm.declare_entry("df", "0.6",
                          Patterns::Double(0.0),
                          "Rate of fibronectin degradation");
        prm.declare_entry("kappa", "20",
                          Patterns::Double(0.0),
                          "Influence of fibroblast on collagen orientation");
        prm.declare_entry("ECM_gridsize", "10",
                          Patterns::Integer(0),
                          "spatial discretization of ECM density/orientation profile");
    }
    prm.leave_subsection();
}

void ECM::parse_parameters(ParameterHandler &prm){
    prm.enter_subsection("ECM");
    {
        pc = prm.get_double("pc");
        dc = prm.get_double("dc");
        df = prm.get_double("df");
        kappa = prm.get_double("kappa");
        ECM_gridsize = prm.get_double("ECM_gridsize");
    }
    prm.leave_subsection();
}

void Endocyte::declare_parameters(ParameterHandler &prm){
    prm.enter_subsection("Endocyte");
    {
        prm.declare_entry("ENinit", "12",
                          Patterns::Integer(0),
                          "Number of initial sprouts");
    }
    prm.leave_subsection();
}

void Endocyte::parse_parameters(ParameterHandler &prm){
    prm.enter_subsection("Endocyte");
    {
        ENinit = prm.get_integer("ENinit");
    }
    prm.leave_subsection();
}

AllParameters::AllParameters(const std::string &intput_file){
    ParameterHandler prm;
    declare_parameters(prm);
    prm.read_input(intput_file);
    parse_parameters(prm);
}

void AllParameters::declare_parameters(ParameterHandler &prm){
    Geometry::declare_parameters(prm);
    Time::declare_parameters(prm);
    Fibroblast::declare_parameters(prm);
    ECM::declare_parameters(prm);
    Endocyte::delcare_parameters(prm);
}

void AllParameters::parse_parameters(ParameterHandler &prm){
    Geometry::parse_parameters(prm);
    Time::parse_parameters(prm);
    Fibroblast::parse_parameters(prm);
    ECM::parse_parameters(prm);
    Endocyte::parse_parameters(prm);
}

}//end namespace Parameters
