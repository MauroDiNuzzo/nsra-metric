// nsra_py.cpp
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include "nsra.cpp"   // the OpenMP optimized core

/* How to compile:
$ cd cpp

$ c++ -O3 -Wall -shared -std=c++17 -fPIC -fopenmp \
    $(python3 -m pybind11 --includes) \
    nsra_py.cpp -o nsra_cpp$(python3-config --extension-suffix)

Or (if you use a different Python in venv, e.g., 3.12 here):

$ c++ -O3 -Wall -shared -std=c++17 -fPIC -fopenmp \
    $(python -m pybind11 --includes) \
    nsra_py.cpp -o nsra_cpp.cpython-312-x86_64-linux-gnu.so   
*/

namespace py = pybind11;

/*
Python usage:
import nsra_cpp
score = nsra_cpp.nsra_cpp(meas, pred, epsilon=0.1, tie_score=0.5)
*/
double nsra_wrapper(py::array_t<double> meas_np,
                    py::array_t<double> pred_np,
                    double epsilon = 0.0,
                    double tie_score = 0.5)
{
    // Check shapes
    if (meas_np.size() != pred_np.size())
        throw std::runtime_error("meas and pred must have the same number of elements");

    size_t n = meas_np.size();

    // Get raw pointers to the data (no copy)
    double* meas = meas_np.mutable_data();
    double* pred = pred_np.mutable_data();

    return nsra_core_omp(meas, pred, n, epsilon, tie_score);
}

PYBIND11_MODULE(nsra_cpp, m) {
    m.doc() = "NSRA: Null-stratified Rank Accuracy C++ OpenMP wrapper";
    m.def("nsra_cpp", &nsra_wrapper,
          py::arg("meas"), py::arg("pred"),
          py::arg("epsilon") = 0.0,
          py::arg("tie_score") = 0.5,
          "Compute NSRA score using OpenMP");
}
