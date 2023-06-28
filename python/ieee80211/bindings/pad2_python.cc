/*
 * Copyright 2023 Free Software Foundation, Inc.
 *
 * This file is part of GNU Radio
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 *
 */

/***********************************************************************************/
/* This file is automatically generated using bindtool and can be manually edited  */
/* The following lines can be configured to regenerate this file during cmake      */
/* If manual edits are made, the following tags should be modified accordingly.    */
/* BINDTOOL_GEN_AUTOMATIC(0)                                                       */
/* BINDTOOL_USE_PYGCCXML(0)                                                        */
/* BINDTOOL_HEADER_FILE(pad2.h)                                        */
/* BINDTOOL_HEADER_FILE_HASH(1638e273624edbf0aa4f390829baf053)                     */
/***********************************************************************************/

#include <pybind11/complex.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

#include <gnuradio/ieee80211/pad2.h>
// pydoc.h is automatically generated in the build directory
#include <pad2_pydoc.h>

void bind_pad2(py::module& m)
{

    using pad2    = gr::ieee80211::pad2;


    py::class_<pad2, gr::block, gr::basic_block,
        std::shared_ptr<pad2>>(m, "pad2", D(pad2))

        .def(py::init(&pad2::make),
           D(pad2,make)
        )
        



        ;




}








