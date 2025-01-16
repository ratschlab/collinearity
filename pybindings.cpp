//
// Created by Sayan Goswami on 16.01.2025.
//

#include "collinearity.h"
#include <pybind11/detail/descr.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

struct rt_index_t {
    const int batch_sz;
    index_t *index = nullptr;
    std::vector<std::string> refnames;
    explicit rt_index_t(int batch_sz): batch_sz(batch_sz) {}
};