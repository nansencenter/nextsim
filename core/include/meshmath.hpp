/* -*- mode: c++; coding: utf-8; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; show-trailing-whitespace: t -*- vim: set fenc=utf-8 ft=cpp et sw=4 ts=4 sts=4: */

/**
 * @file   meshmath.cpp
 * @author Abdoulaye Samake <abdoulaye.samake@nersc.no>
 * @author Sylvain Bouillon <sylvain.bouillon@nersc.no>
 * @author Einar Olason <einar.olason@nersc.no>
 * @date   Wed 17 Aug 2022 07:59:06 CEST
 */

#ifndef __MeshMath_HPP
#define __MeshMath_HPP 1


namespace Nextsim
{

    class MeshMath
    {
        public:

            //------------------------------------------------------------------------------------------------------
            //! Calculates the Jacobian Matrix Determinant:  measure of the normals of the element faces relative to each other.
            //! It is the determinant of the transformation from the reference triangle with vertices
            //! (0,0), (1,0) and (0,1) to an arbitrary triangle.
            //! This transformation is:
            //!   x=x0+(x1-x0)\xi + (x2-x1)\eta,
            //!   y=y0+(y1-y0)\xi + (y2-y1)\eta,
            //! with \xi,\eta between 0 and 1, so the determinant is:
            //!   (x1-x0)*(y2-y0) - (x2-x0)*(y1-y0).
            //! Called by the flip(), measure() and shapeCoeff() functions.
            //! \note
            //! * This is used to calculate the finite element shape coefficient.
            //! * The Jacobian is an indicator of the distortion of the current mesh
            //!   with respect to an undistorted mesh.
            static double
            jacobian(std::vector<std::vector<double>> const& vertices)
            {
                double jac = (vertices[1][0]-vertices[0][0])*(vertices[2][1]-vertices[0][1]);
                jac -= (vertices[2][0]-vertices[0][0])*(vertices[1][1]-vertices[0][1]);
                return jac;
            }; //jacobian
    };
}

#endif
