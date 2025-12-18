/* -------------------------------------------------------------------- */
/* Copyright (C) 1991 - 2026 - EDF - www.code-aster.org             */
/* This file is part of code_aster.                                     */
/*                                                                      */
/* code_aster is free software: you can redistribute it and/or modify   */
/* it under the terms of the GNU General Public License as published by */
/* the Free Software Foundation, either version 3 of the License, or    */
/* (at your option) any later version.                                  */
/*                                                                      */
/* code_aster is distributed in the hope that it will be useful,        */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of       */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        */
/* GNU General Public License for more details.                         */
/*                                                                      */
/* You should have received a copy of the GNU General Public License    */
/* along with code_aster.  If not, see <http://www.gnu.org/licenses/>.  */
/* -------------------------------------------------------------------- */

#include "astercxx.h"

#include "tria_mitc.h"

#include <cmath>
#include <iomanip> // Pour setprecision
#include <iostream>
#include <vector>

#include <math.h>
#include <stdalign.h>
#include <stdlib.h>
#include <string.h>

VectorReal B_p1_tr6( const VectorReal &w, const VectorReal &coordinate_dofs,
                     const VectorInt &entity_local_index, const VectorReal &c ) {
    VectorReal A( 42 * 42, 0.0 ); // Vecteur 1D pour stocker les valeurs (size = nbddl)
    // Quadrature rules
    static const double weights_48e[3] = { 0.1666666666666667, 0.1666666666666667,
                                           0.1666666666666667 };
    // Precomputed values of basis functions and precomputations
    // FE* dimensions: [permutation][entities][points][dofs]
    static const double FE4_C1_D01_Q48e[1][1][3][6] = {
        { { { -1.666666666666667, 0.0, -0.3333333333333333, 0.6666666666666667, 2.0,
              -0.6666666666666669 },
            { 0.3333333333333328, 0.0, 1.666666666666666, 0.6666666666666669, -2.0,
              -0.6666666666666662 },
            { 0.333333333333333, 0.0, -0.3333333333333335, 2.666666666666667, 0.0,
              -2.666666666666666 } } }
    };
    static const double FE4_C1_D10_Q48e[1][1][3][6] = {
        { { { -1.666666666666667, -0.333333333333333, 0.0, 0.6666666666666664, -0.6666666666666667,
              2.0 },
            { 0.3333333333333325, -0.3333333333333328, 0.0, 2.666666666666667, -2.666666666666667,
              0.0 },
            { 0.3333333333333331, 1.666666666666667, 0.0, 0.6666666666666666, -0.6666666666666663,
              -1.999999999999999 } } }
    };
    static const double FE4_C2_Q48e[1][1][3][3] = {
        { { { 0.6666666666666667, 0.1666666666666666, 0.1666666666666667 },
            { 0.1666666666666667, 0.1666666666666666, 0.6666666666666665 },
            { 0.1666666666666668, 0.6666666666666665, 0.1666666666666667 } } }
    };
    static const double FE4_C3_Q48e[1][1][3][3] = {
        { { { -0.1666666666666667, 0.1666666666666667, 0.8333333333333333 },
            { -0.6666666666666667, 0.6666666666666667, 0.3333333333333333 },
            { -0.1666666666666666, 0.1666666666666666, 0.8333333333333335 } } }
    };
    static const double FE4_C4_Q48e[1][1][3][3] = {
        { { { 0.1666666666666667, 0.8333333333333333, 0.1666666666666667 },
            { 0.1666666666666666, 0.8333333333333335, 0.1666666666666666 },
            { 0.6666666666666667, 0.3333333333333333, 0.6666666666666667 } } }
    };
    static const double FE6_C0_D10_Q48e[1][1][1][3] = { { { { -1.0, 1.0, 0.0 } } } };
    static const double FE6_C1_D01_Q48e[1][1][1][3] = { { { { -1.0, 0.0, 1.0 } } } };
    // ------------------------
    // Section: Jacobian
    // Inputs: FE6_C0_D10_Q48e, FE6_C1_D01_Q48e, coordinate_dofs
    // Outputs: J_c2, J_c3, J_c0, J_c1
    double J_c0 = 0.0;
    double J_c3 = 0.0;
    double J_c1 = 0.0;
    double J_c2 = 0.0;
    {
        for ( int ic = 0; ic < 3; ++ic ) {
            J_c0 += coordinate_dofs[(ic)*3] * FE6_C0_D10_Q48e[0][0][0][ic];
            J_c3 += coordinate_dofs[(ic)*3 + 1] * FE6_C1_D01_Q48e[0][0][0][ic];
            J_c1 += coordinate_dofs[(ic)*3] * FE6_C1_D01_Q48e[0][0][0][ic];
            J_c2 += coordinate_dofs[(ic)*3 + 1] * FE6_C0_D10_Q48e[0][0][0][ic];
        }
    }
    // ------------------------
    double sp_48e_0 = J_c0 * J_c3;
    double sp_48e_1 = J_c1 * J_c2;
    double sp_48e_2 = -sp_48e_1;
    double sp_48e_3 = sp_48e_0 + sp_48e_2;
    double sp_48e_4 = J_c0 / sp_48e_3;
    double sp_48e_5 = -J_c1;
    double sp_48e_6 = sp_48e_5 / sp_48e_3;
    double sp_48e_7 = J_c3 / sp_48e_3;
    double sp_48e_8 = -J_c2;
    double sp_48e_9 = sp_48e_8 / sp_48e_3;
    double sp_48e_10 = c[0] * c[2];
    double sp_48e_11 = c[3] * sp_48e_10;
    double sp_48e_12 = 1.0 + c[1];
    double sp_48e_13 = 4.0 * sp_48e_12;
    double sp_48e_14 = sp_48e_11 / sp_48e_13;
    double sp_48e_15 = sp_48e_4 + sp_48e_4;
    double sp_48e_16 = sp_48e_6 + sp_48e_6;
    double sp_48e_17 = sp_48e_15 / 2;
    double sp_48e_18 = sp_48e_16 / 2;
    double sp_48e_19 = sp_48e_9 / 2;
    double sp_48e_20 = sp_48e_7 / 2;
    double sp_48e_21 = sp_48e_4 / 2;
    double sp_48e_22 = sp_48e_6 / 2;
    double sp_48e_23 = sp_48e_9 + sp_48e_9;
    double sp_48e_24 = sp_48e_7 + sp_48e_7;
    double sp_48e_25 = sp_48e_23 / 2;
    double sp_48e_26 = sp_48e_24 / 2;
    double sp_48e_27 = -c[1];
    double sp_48e_28 = 1.0 + sp_48e_27;
    double sp_48e_29 = 2 * sp_48e_17;
    double sp_48e_30 = 2 * sp_48e_18;
    double sp_48e_31 = 2 * sp_48e_25;
    double sp_48e_32 = 2 * sp_48e_26;
    double sp_48e_33 = pow( c[3], 3 );
    double sp_48e_34 = c[0] * sp_48e_33;
    double sp_48e_35 = pow( c[1], 2 );
    double sp_48e_36 = -sp_48e_35;
    double sp_48e_37 = 1.0 + sp_48e_36;
    double sp_48e_38 = 24.0 * sp_48e_37;
    double sp_48e_39 = sp_48e_34 / sp_48e_38;
    double sp_48e_40 = fabs( sp_48e_3 );
    double sp_48e_41 = c[4] * sp_48e_33;
    double sp_48e_42 = -sp_48e_41;
    double sp_48e_43 = sp_48e_42 * sp_48e_40;
    for ( int iq = 0; iq < 3; ++iq ) {
        // ------------------------
        // Section: Function
        // Inputs: w, FE4_C4_Q48e
        // Outputs: w0_c4
        double w0_c4 = 0.0;
        {
            for ( int ic = 0; ic < 3; ++ic ) {
                w0_c4 += w[( ic ) + 15] * FE4_C4_Q48e[0][0][iq][ic];
            }
        }
        // ------------------------
        // ------------------------
        // Section: Function
        // Inputs: w, FE4_C3_Q48e
        // Outputs: w0_c3
        double w0_c3 = 0.0;
        {
            for ( int ic = 0; ic < 3; ++ic ) {
                w0_c3 += w[( ic ) + 15] * FE4_C3_Q48e[0][0][iq][ic];
            }
        }
        // ------------------------
        // ------------------------
        // Section: Function
        // Inputs: w, FE4_C1_D01_Q48e
        // Outputs: w0_d01_c1
        double w0_d01_c1 = 0.0;
        {
            for ( int ic = 0; ic < 6; ++ic ) {
                w0_d01_c1 += w[(ic)*2 + 1] * FE4_C1_D01_Q48e[0][0][iq][ic];
            }
        }
        // ------------------------
        // ------------------------
        // Section: Function
        // Inputs: w, FE4_C1_D10_Q48e
        // Outputs: w0_d10_c1
        double w0_d10_c1 = 0.0;
        {
            for ( int ic = 0; ic < 6; ++ic ) {
                w0_d10_c1 += w[(ic)*2 + 1] * FE4_C1_D10_Q48e[0][0][iq][ic];
            }
        }
        // ------------------------
        // ------------------------
        // Section: Function
        // Inputs: w, FE4_C1_D01_Q48e
        // Outputs: w0_d01_c0
        double w0_d01_c0 = 0.0;
        {
            for ( int ic = 0; ic < 6; ++ic ) {
                w0_d01_c0 += w[(ic)*2] * FE4_C1_D01_Q48e[0][0][iq][ic];
            }
        }
        // ------------------------
        // ------------------------
        // Section: Function
        // Inputs: w, FE4_C1_D10_Q48e
        // Outputs: w0_d10_c0
        double w0_d10_c0 = 0.0;
        {
            for ( int ic = 0; ic < 6; ++ic ) {
                w0_d10_c0 += w[(ic)*2] * FE4_C1_D10_Q48e[0][0][iq][ic];
            }
        }
        // ------------------------
        // ------------------------
        // Section: Intermediates
        // Inputs: w0_c4, w0_c3, w0_d01_c1, w0_d10_c1, w0_d01_c0, w0_d10_c0
        // Outputs: fw0, fw1, fw2, fw3, fw4, fw5, fw6
        double fw0 = 0;
        double fw1 = 0;
        double fw2 = 0;
        double fw3 = 0;
        double fw4 = 0;
        double fw5 = 0;
        double fw6 = 0;
        {
            double sv_48e_0 = w0_c4 * sp_48e_4;
            double sv_48e_1 = w0_c3 * sp_48e_6;
            double sv_48e_2 = sv_48e_0 + sv_48e_1;
            double sv_48e_3 = sv_48e_2 * sp_48e_4;
            double sv_48e_4 = sv_48e_2 * sp_48e_6;
            double sv_48e_5 = sv_48e_3 + sv_48e_3;
            double sv_48e_6 = sv_48e_4 + sv_48e_4;
            double sv_48e_7 = w0_c3 * sp_48e_7;
            double sv_48e_8 = w0_c4 * sp_48e_9;
            double sv_48e_9 = sv_48e_7 + sv_48e_8;
            double sv_48e_10 = sv_48e_9 * sp_48e_9;
            double sv_48e_11 = sv_48e_9 * sp_48e_7;
            double sv_48e_12 = sv_48e_10 + sv_48e_10;
            double sv_48e_13 = sv_48e_11 + sv_48e_11;
            double sv_48e_14 = sv_48e_5 + sv_48e_12;
            double sv_48e_15 = sv_48e_13 + sv_48e_6;
            double sv_48e_16 = sv_48e_14 * sp_48e_14;
            double sv_48e_17 = sv_48e_15 * sp_48e_14;
            double sv_48e_18 = w0_d01_c1 * sp_48e_4;
            double sv_48e_19 = w0_d10_c1 * sp_48e_6;
            double sv_48e_20 = sv_48e_18 + sv_48e_19;
            double sv_48e_21 = sv_48e_20 + sv_48e_20;
            double sv_48e_22 = sv_48e_21 / 2;
            double sv_48e_23 = sv_48e_22 * sp_48e_17;
            double sv_48e_24 = sv_48e_22 * sp_48e_18;
            double sv_48e_25 = sv_48e_23 + sv_48e_23;
            double sv_48e_26 = sv_48e_24 + sv_48e_24;
            double sv_48e_27 = w0_d01_c0 * sp_48e_4;
            double sv_48e_28 = w0_d10_c0 * sp_48e_6;
            double sv_48e_29 = sv_48e_27 + sv_48e_28;
            double sv_48e_30 = w0_d10_c1 * sp_48e_7;
            double sv_48e_31 = w0_d01_c1 * sp_48e_9;
            double sv_48e_32 = sv_48e_30 + sv_48e_31;
            double sv_48e_33 = sv_48e_29 + sv_48e_32;
            double sv_48e_34 = sv_48e_33 / 2;
            double sv_48e_35 = sv_48e_34 * sp_48e_19;
            double sv_48e_36 = sv_48e_34 * sp_48e_20;
            double sv_48e_37 = sv_48e_34 * sp_48e_21;
            double sv_48e_38 = sv_48e_34 * sp_48e_22;
            double sv_48e_39 = sv_48e_35 + sv_48e_35;
            double sv_48e_40 = sv_48e_36 + sv_48e_36;
            double sv_48e_41 = sv_48e_37 + sv_48e_37;
            double sv_48e_42 = sv_48e_38 + sv_48e_38;
            double sv_48e_43 = sv_48e_25 + sv_48e_39;
            double sv_48e_44 = sv_48e_26 + sv_48e_40;
            double sv_48e_45 = w0_d10_c0 * sp_48e_7;
            double sv_48e_46 = w0_d01_c0 * sp_48e_9;
            double sv_48e_47 = sv_48e_45 + sv_48e_46;
            double sv_48e_48 = sv_48e_47 + sv_48e_47;
            double sv_48e_49 = sv_48e_48 / 2;
            double sv_48e_50 = sv_48e_49 * sp_48e_25;
            double sv_48e_51 = sv_48e_49 * sp_48e_26;
            double sv_48e_52 = sv_48e_50 + sv_48e_50;
            double sv_48e_53 = sv_48e_51 + sv_48e_51;
            double sv_48e_54 = sv_48e_52 + sv_48e_41;
            double sv_48e_55 = sv_48e_53 + sv_48e_42;
            double sv_48e_56 = sv_48e_43 + sv_48e_39;
            double sv_48e_57 = sv_48e_44 + sv_48e_40;
            double sv_48e_58 = sv_48e_54 + sv_48e_41;
            double sv_48e_59 = sv_48e_55 + sv_48e_42;
            double sv_48e_60 = sv_48e_56 * sp_48e_28;
            double sv_48e_61 = sv_48e_57 * sp_48e_28;
            double sv_48e_62 = sv_48e_58 * sp_48e_28;
            double sv_48e_63 = sv_48e_59 * sp_48e_28;
            double sv_48e_64 = sv_48e_22 + sv_48e_49;
            double sv_48e_65 = sv_48e_64 * sp_48e_29;
            double sv_48e_66 = sv_48e_64 * sp_48e_30;
            double sv_48e_67 = sv_48e_64 * sp_48e_31;
            double sv_48e_68 = sv_48e_64 * sp_48e_32;
            double sv_48e_69 = c[1] * sv_48e_65;
            double sv_48e_70 = c[1] * sv_48e_66;
            double sv_48e_71 = c[1] * sv_48e_67;
            double sv_48e_72 = c[1] * sv_48e_68;
            double sv_48e_73 = sv_48e_60 + sv_48e_69;
            double sv_48e_74 = sv_48e_61 + sv_48e_70;
            double sv_48e_75 = sv_48e_62 + sv_48e_71;
            double sv_48e_76 = sv_48e_63 + sv_48e_72;
            double sv_48e_77 = sv_48e_73 * sp_48e_39;
            double sv_48e_78 = sv_48e_74 * sp_48e_39;
            double sv_48e_79 = sv_48e_75 * sp_48e_39;
            double sv_48e_80 = sv_48e_76 * sp_48e_39;
            double sv_48e_81 = sv_48e_16 * sp_48e_40;
            double sv_48e_82 = sv_48e_17 * sp_48e_40;
            double sv_48e_83 = sv_48e_77 * sp_48e_40;
            double sv_48e_84 = sv_48e_78 * sp_48e_40;
            double sv_48e_85 = sv_48e_79 * sp_48e_40;
            double sv_48e_86 = sv_48e_80 * sp_48e_40;
            fw0 = sv_48e_86 * weights_48e[iq];
            fw1 = sv_48e_85 * weights_48e[iq];
            fw2 = sv_48e_84 * weights_48e[iq];
            fw3 = sv_48e_83 * weights_48e[iq];
            fw4 = sp_48e_43 * weights_48e[iq];
            fw5 = sv_48e_82 * weights_48e[iq];
            fw6 = sv_48e_81 * weights_48e[iq];
        }
        // ------------------------
        // ------------------------
        // Section: Tensor Computation
        // Inputs: fw3, fw2, fw1, FE4_C1_D01_Q48e, fw0, FE4_C1_D10_Q48e
        // Outputs: A
        {
            for ( int i = 0; i < 6; ++i ) {
                A[2 * ( i )] += fw0 * FE4_C1_D10_Q48e[0][0][iq][i];
                A[2 * ( i )] += fw1 * FE4_C1_D01_Q48e[0][0][iq][i];
                A[( 2 * ( i ) + 1 )] += fw2 * FE4_C1_D10_Q48e[0][0][iq][i];
                A[( 2 * ( i ) + 1 )] += fw3 * FE4_C1_D01_Q48e[0][0][iq][i];
            }
        }
        // ------------------------
        // ------------------------
        // Section: Tensor Computation
        // Inputs: fw4, FE4_C4_Q48e, fw6, FE4_C3_Q48e, FE4_C2_Q48e, fw5
        // Outputs: A
        {
            for ( int i = 0; i < 3; ++i ) {
                A[( ( i ) + 12 )] += fw4 * FE4_C2_Q48e[0][0][iq][i];
                A[( ( i ) + 15 )] += fw5 * FE4_C3_Q48e[0][0][iq][i];
                A[( ( i ) + 15 )] += fw6 * FE4_C4_Q48e[0][0][iq][i];
            }
        }
        // ------------------------
    }

    return A;
}

VectorReal B_p2_tr6( const VectorReal &w, const VectorReal &coordinate_dofs,
                     const VectorInt &entity_local_index, const VectorReal &c ) {
    VectorReal A( 42 * 42, 0.0 ); // Vecteur 1D pour stocker les valeurs (size = nbddl)
    // Quadrature rules
    static const double weights_4a8[2] = { 0.5, 0.5 };
    // Precomputed values of basis functions and precomputations
    // FE* dimensions: [permutation][entities][points][dofs]
    static const double FE1_C1_D01_F_Q4a8[1][1][1][3] = { { { { -1.0, 0.0, 1.0 } } } };
    static const double FE5_C0_F_Q4a8[1][3][2][6] = {
        { { { 0.0, 0.4553418012614797, -0.1220084679281462, 0.6666666666666667, 0.0, 0.0 },
            { 0.0, -0.1220084679281461, 0.4553418012614795, 0.6666666666666667, 0.0, 0.0 } },
          { { 0.4553418012614797, 0.0, -0.1220084679281462, 0.0, 0.6666666666666669, 0.0 },
            { -0.1220084679281461, 0.0, 0.4553418012614795, 0.0, 0.6666666666666667, 0.0 } },
          { { 0.4553418012614794, -0.1220084679281462, 0.0, 0.0, 0.0, 0.6666666666666665 },
            { -0.1220084679281462, 0.4553418012614795, 0.0, 0.0, 0.0, 0.6666666666666665 } } }
    };
    static const double FE5_C2_D10_F_Q4a8[1][1][1][3] = { { { { -1.0, 1.0, 0.0 } } } };
    static const double FE5_C3_F_Q4a8[1][3][2][3] = {
        { { { -0.211324865405187, 0.211324865405187, 0.788675134594813 },
            { -0.7886751345948131, 0.7886751345948131, 0.211324865405187 } },
          { { -0.2113248654051873, 0.2113248654051873, 0.7886751345948128 },
            { -0.7886751345948131, 0.7886751345948131, 0.2113248654051869 } },
          { { 0.0, 0.0, 1.0 }, { 0.0, 0.0, 1.0 } } }
    };
    static const double FE5_C4_F_Q4a8[1][3][2][3] = {
        { { { 0.788675134594813, 0.211324865405187, 0.788675134594813 },
            { 0.211324865405187, 0.788675134594813, 0.211324865405187 } },
          { { 0.0, 1.0, 0.0 }, { 0.0, 1.0, 0.0 } },
          { { 0.2113248654051872, 0.7886751345948128, 0.2113248654051872 },
            { 0.788675134594813, 0.211324865405187, 0.788675134594813 } } }
    };
    static const double triangle_reference_facet_jacobian[3][2][1] = { { { -1.0 }, { 1.0 } },
                                                                       { { 0.0 }, { 1.0 } },
                                                                       { { 1.0 }, { 0.0 } } };
    static const double triangle_reference_facet_normals[3][2] = {
        { 0.7071067811865475, 0.7071067811865475 }, { -1.0, -0.0 }, { 0.0, -1.0 }
    };
    // ------------------------
    // Section: Function
    // Inputs: w, FE5_C2_D10_F_Q4a8
    // Outputs: w0_d10_c2
    double w0_d10_c2 = 0.0;
    {
        for ( int ic = 0; ic < 3; ++ic ) {
            w0_d10_c2 += w[( ic ) + 12] * FE5_C2_D10_F_Q4a8[0][0][0][ic];
        }
    }
    // ------------------------
    // ------------------------
    // Section: Jacobian
    // Inputs: coordinate_dofs, FE1_C1_D01_F_Q4a8, FE5_C2_D10_F_Q4a8
    // Outputs: J_c1, J_c2, J_c0, J_c3
    double J_c3 = 0.0;
    double J_c0 = 0.0;
    double J_c1 = 0.0;
    double J_c2 = 0.0;
    {
        for ( int ic = 0; ic < 3; ++ic ) {
            J_c3 += coordinate_dofs[(ic)*3 + 1] * FE1_C1_D01_F_Q4a8[0][0][0][ic];
            J_c0 += coordinate_dofs[(ic)*3] * FE5_C2_D10_F_Q4a8[0][0][0][ic];
            J_c1 += coordinate_dofs[(ic)*3] * FE1_C1_D01_F_Q4a8[0][0][0][ic];
            J_c2 += coordinate_dofs[(ic)*3 + 1] * FE5_C2_D10_F_Q4a8[0][0][0][ic];
        }
    }
    // ------------------------
    // ------------------------
    // Section: Function
    // Inputs: w, FE1_C1_D01_F_Q4a8
    // Outputs: w0_d01_c2
    double w0_d01_c2 = 0.0;
    {
        for ( int ic = 0; ic < 3; ++ic ) {
            w0_d01_c2 += w[( ic ) + 12] * FE1_C1_D01_F_Q4a8[0][0][0][ic];
        }
    }
    // ------------------------
    double sp_4a8_0 = J_c0 * J_c3;
    double sp_4a8_1 = J_c1 * J_c2;
    double sp_4a8_2 = -sp_4a8_1;
    double sp_4a8_3 = sp_4a8_0 + sp_4a8_2;
    double sp_4a8_4 = J_c3 / sp_4a8_3;
    double sp_4a8_5 = w0_d10_c2 * sp_4a8_4;
    double sp_4a8_6 = -J_c2;
    double sp_4a8_7 = sp_4a8_6 / sp_4a8_3;
    double sp_4a8_8 = w0_d01_c2 * sp_4a8_7;
    double sp_4a8_9 = sp_4a8_5 + sp_4a8_8;
    double sp_4a8_10 = J_c0 / sp_4a8_3;
    double sp_4a8_11 = triangle_reference_facet_normals[entity_local_index[0]][1] * sp_4a8_10;
    double sp_4a8_12 = -J_c1;
    double sp_4a8_13 = sp_4a8_12 / sp_4a8_3;
    double sp_4a8_14 = triangle_reference_facet_normals[entity_local_index[0]][0] * sp_4a8_13;
    double sp_4a8_15 = sp_4a8_11 + sp_4a8_14;
    double sp_4a8_16 = sp_4a8_15 * sp_4a8_15;
    double sp_4a8_17 = triangle_reference_facet_normals[entity_local_index[0]][0] * sp_4a8_4;
    double sp_4a8_18 = triangle_reference_facet_normals[entity_local_index[0]][1] * sp_4a8_7;
    double sp_4a8_19 = sp_4a8_17 + sp_4a8_18;
    double sp_4a8_20 = sp_4a8_19 * sp_4a8_19;
    double sp_4a8_21 = sp_4a8_16 + sp_4a8_20;
    double sp_4a8_22 = sqrt( sp_4a8_21 );
    double sp_4a8_23 = sp_4a8_15 / sp_4a8_22;
    double sp_4a8_24 = -sp_4a8_23;
    double sp_4a8_25 = w0_d01_c2 * sp_4a8_10;
    double sp_4a8_26 = w0_d10_c2 * sp_4a8_13;
    double sp_4a8_27 = sp_4a8_25 + sp_4a8_26;
    double sp_4a8_28 = sp_4a8_19 / sp_4a8_22;
    double sp_4a8_29 = sp_4a8_24 * sp_4a8_4;
    double sp_4a8_30 = sp_4a8_24 * sp_4a8_7;
    double sp_4a8_31 = sp_4a8_13 * sp_4a8_28;
    double sp_4a8_32 = sp_4a8_10 * sp_4a8_28;
    double sp_4a8_33 = sp_4a8_29 + sp_4a8_31;
    double sp_4a8_34 = sp_4a8_30 + sp_4a8_32;
    double sp_4a8_35 = -sp_4a8_4;
    double sp_4a8_36 = -sp_4a8_7;
    double sp_4a8_37 = -sp_4a8_24;
    double sp_4a8_38 = sp_4a8_35 * sp_4a8_24;
    double sp_4a8_39 = sp_4a8_36 * sp_4a8_24;
    double sp_4a8_40 = -sp_4a8_13;
    double sp_4a8_41 = -sp_4a8_10;
    double sp_4a8_42 = sp_4a8_40 * sp_4a8_28;
    double sp_4a8_43 = sp_4a8_41 * sp_4a8_28;
    double sp_4a8_44 = -sp_4a8_28;
    double sp_4a8_45 = sp_4a8_38 + sp_4a8_42;
    double sp_4a8_46 = sp_4a8_39 + sp_4a8_43;
    double sp_4a8_47 = J_c0 * triangle_reference_facet_jacobian[entity_local_index[0]][0][0];
    double sp_4a8_48 = J_c1 * triangle_reference_facet_jacobian[entity_local_index[0]][1][0];
    double sp_4a8_49 = sp_4a8_47 + sp_4a8_48;
    double sp_4a8_50 = sp_4a8_49 * sp_4a8_49;
    double sp_4a8_51 = triangle_reference_facet_jacobian[entity_local_index[0]][0][0] * J_c2;
    double sp_4a8_52 = triangle_reference_facet_jacobian[entity_local_index[0]][1][0] * J_c3;
    double sp_4a8_53 = sp_4a8_51 + sp_4a8_52;
    double sp_4a8_54 = sp_4a8_53 * sp_4a8_53;
    double sp_4a8_55 = sp_4a8_50 + sp_4a8_54;
    double sp_4a8_56 = sqrt( sp_4a8_55 );
    for ( int iq = 0; iq < 2; ++iq ) {
        // ------------------------
        // Section: Function
        // Inputs: w, FE5_C0_F_Q4a8
        // Outputs: w0_c0
        double w0_c0 = 0.0;
        {
            for ( int ic = 0; ic < 6; ++ic ) {
                w0_c0 += w[(ic)*2] * FE5_C0_F_Q4a8[0][entity_local_index[0]][iq][ic];
            }
        }
        // ------------------------
        // ------------------------
        // Section: Function
        // Inputs: w, FE5_C3_F_Q4a8
        // Outputs: w0_c3
        double w0_c3 = 0.0;
        {
            for ( int ic = 0; ic < 3; ++ic ) {
                w0_c3 += w[( ic ) + 15] * FE5_C3_F_Q4a8[0][entity_local_index[0]][iq][ic];
            }
        }
        // ------------------------
        // ------------------------
        // Section: Function
        // Inputs: w, FE5_C4_F_Q4a8
        // Outputs: w0_c4
        double w0_c4 = 0.0;
        {
            for ( int ic = 0; ic < 3; ++ic ) {
                w0_c4 += w[( ic ) + 15] * FE5_C4_F_Q4a8[0][entity_local_index[0]][iq][ic];
            }
        }
        // ------------------------
        // ------------------------
        // Section: Function
        // Inputs: w, FE5_C0_F_Q4a8
        // Outputs: w0_c1
        double w0_c1 = 0.0;
        {
            for ( int ic = 0; ic < 6; ++ic ) {
                w0_c1 += w[(ic)*2 + 1] * FE5_C0_F_Q4a8[0][entity_local_index[0]][iq][ic];
            }
        }
        // ------------------------
        // ------------------------
        // Section: Function
        // Inputs: w, FE5_C3_F_Q4a8
        // Outputs: w0_c5
        double w0_c5 = 0.0;
        {
            for ( int ic = 0; ic < 3; ++ic ) {
                w0_c5 += w[( ic ) + 18] * FE5_C3_F_Q4a8[0][entity_local_index[0]][iq][ic];
            }
        }
        // ------------------------
        // ------------------------
        // Section: Function
        // Inputs: w, FE5_C4_F_Q4a8
        // Outputs: w0_c6
        double w0_c6 = 0.0;
        {
            for ( int ic = 0; ic < 3; ++ic ) {
                w0_c6 += w[( ic ) + 18] * FE5_C4_F_Q4a8[0][entity_local_index[0]][iq][ic];
            }
        }
        // ------------------------
        // ------------------------
        // Section: Intermediates
        // Inputs: w0_c0, w0_c3, w0_c4, w0_c1, w0_c5, w0_c6
        // Outputs: fw0, fw1, fw2, fw3, fw4, fw5, fw6, fw7
        double fw0 = 0;
        double fw1 = 0;
        double fw2 = 0;
        double fw3 = 0;
        double fw4 = 0;
        double fw5 = 0;
        double fw6 = 0;
        double fw7 = 0;
        {
            double sv_4a8_0 = -w0_c0;
            double sv_4a8_1 = sp_4a8_9 + sv_4a8_0;
            double sv_4a8_2 = w0_c3 * sp_4a8_4;
            double sv_4a8_3 = w0_c4 * sp_4a8_7;
            double sv_4a8_4 = sv_4a8_2 + sv_4a8_3;
            double sv_4a8_5 = -sv_4a8_4;
            double sv_4a8_6 = sv_4a8_1 + sv_4a8_5;
            double sv_4a8_7 = sv_4a8_6 * sp_4a8_24;
            double sv_4a8_8 = -w0_c1;
            double sv_4a8_9 = sp_4a8_27 + sv_4a8_8;
            double sv_4a8_10 = w0_c4 * sp_4a8_10;
            double sv_4a8_11 = w0_c3 * sp_4a8_13;
            double sv_4a8_12 = sv_4a8_10 + sv_4a8_11;
            double sv_4a8_13 = -sv_4a8_12;
            double sv_4a8_14 = sv_4a8_9 + sv_4a8_13;
            double sv_4a8_15 = sv_4a8_14 * sp_4a8_28;
            double sv_4a8_16 = sv_4a8_7 + sv_4a8_15;
            double sv_4a8_17 = sv_4a8_16 * sp_4a8_33;
            double sv_4a8_18 = sv_4a8_16 * sp_4a8_34;
            double sv_4a8_19 = w0_c5 * sp_4a8_4;
            double sv_4a8_20 = w0_c6 * sp_4a8_7;
            double sv_4a8_21 = sv_4a8_19 + sv_4a8_20;
            double sv_4a8_22 = sv_4a8_21 * sp_4a8_24;
            double sv_4a8_23 = w0_c6 * sp_4a8_10;
            double sv_4a8_24 = w0_c5 * sp_4a8_13;
            double sv_4a8_25 = sv_4a8_23 + sv_4a8_24;
            double sv_4a8_26 = sv_4a8_25 * sp_4a8_28;
            double sv_4a8_27 = sv_4a8_22 + sv_4a8_26;
            double sv_4a8_28 = sv_4a8_27 * sp_4a8_33;
            double sv_4a8_29 = sv_4a8_27 * sp_4a8_34;
            double sv_4a8_30 = sv_4a8_27 * sp_4a8_37;
            double sv_4a8_31 = sv_4a8_27 * sp_4a8_45;
            double sv_4a8_32 = sv_4a8_27 * sp_4a8_46;
            double sv_4a8_33 = sv_4a8_27 * sp_4a8_44;
            double sv_4a8_34 = sv_4a8_17 * sp_4a8_56;
            double sv_4a8_35 = sv_4a8_18 * sp_4a8_56;
            double sv_4a8_36 = sv_4a8_28 * sp_4a8_56;
            double sv_4a8_37 = sv_4a8_29 * sp_4a8_56;
            double sv_4a8_38 = sv_4a8_30 * sp_4a8_56;
            double sv_4a8_39 = sv_4a8_31 * sp_4a8_56;
            double sv_4a8_40 = sv_4a8_32 * sp_4a8_56;
            double sv_4a8_41 = sv_4a8_33 * sp_4a8_56;
            fw0 = sv_4a8_38 * weights_4a8[iq];
            fw1 = sv_4a8_41 * weights_4a8[iq];
            fw2 = sv_4a8_36 * weights_4a8[iq];
            fw3 = sv_4a8_37 * weights_4a8[iq];
            fw4 = sv_4a8_39 * weights_4a8[iq];
            fw5 = sv_4a8_40 * weights_4a8[iq];
            fw6 = sv_4a8_34 * weights_4a8[iq];
            fw7 = sv_4a8_35 * weights_4a8[iq];
        }
        // ------------------------
        // ------------------------
        // Section: Tensor Computation
        // Inputs: fw1, fw0, FE5_C0_F_Q4a8
        // Outputs: A
        {
            for ( int i = 0; i < 6; ++i ) {
                A[2 * ( i )] += fw0 * FE5_C0_F_Q4a8[0][entity_local_index[0]][iq][i];
                A[( 2 * ( i ) + 1 )] += fw1 * FE5_C0_F_Q4a8[0][entity_local_index[0]][iq][i];
            }
        }
        // ------------------------
        // ------------------------
        // Section: Tensor Computation
        // Inputs: fw2, fw7, fw5, FE5_C3_F_Q4a8, FE5_C4_F_Q4a8, fw4, fw6, fw3, FE5_C2_D10_F_Q4a8,
        // FE1_C1_D01_F_Q4a8 Outputs: A
        {
            for ( int i = 0; i < 3; ++i ) {
                A[( ( i ) + 12 )] += fw2 * FE5_C2_D10_F_Q4a8[0][0][0][i];
                A[( ( i ) + 12 )] += fw3 * FE1_C1_D01_F_Q4a8[0][0][0][i];
                A[( ( i ) + 15 )] += fw4 * FE5_C3_F_Q4a8[0][entity_local_index[0]][iq][i];
                A[( ( i ) + 15 )] += fw5 * FE5_C4_F_Q4a8[0][entity_local_index[0]][iq][i];
                A[( ( i ) + 18 )] += fw6 * FE5_C3_F_Q4a8[0][entity_local_index[0]][iq][i];
                A[( ( i ) + 18 )] += fw7 * FE5_C4_F_Q4a8[0][entity_local_index[0]][iq][i];
            }
        }
        // ------------------------
    }

    return A;
}

VectorReal B_p4_tr6( const VectorReal &w, const VectorReal &coordinate_dofs,
                     const VectorInt &entity_local_index, const VectorReal &c ) {
    VectorReal A( 42 * 42, 0.0 ); // Vecteur 1D pour stocker les valeurs (size = nbddl)
    // Quadrature rules
    static const double weights_48e[3] = { 0.1666666666666667, 0.1666666666666667,
                                           0.1666666666666667 };
    // Precomputed values of basis functions and precomputations
    // FE* dimensions: [permutation][entities][points][dofs]
    static const double FE4_C1_D01_Q48e[1][1][3][6] = {
        { { { -1.666666666666667, 0.0, -0.3333333333333333, 0.6666666666666667, 2.0,
              -0.6666666666666669 },
            { 0.3333333333333328, 0.0, 1.666666666666666, 0.6666666666666669, -2.0,
              -0.6666666666666662 },
            { 0.333333333333333, 0.0, -0.3333333333333335, 2.666666666666667, 0.0,
              -2.666666666666666 } } }
    };
    static const double FE4_C1_D10_Q48e[1][1][3][6] = {
        { { { -1.666666666666667, -0.333333333333333, 0.0, 0.6666666666666664, -0.6666666666666667,
              2.0 },
            { 0.3333333333333325, -0.3333333333333328, 0.0, 2.666666666666667, -2.666666666666667,
              0.0 },
            { 0.3333333333333331, 1.666666666666667, 0.0, 0.6666666666666666, -0.6666666666666663,
              -1.999999999999999 } } }
    };
    static const double FE4_C3_Q48e[1][1][3][3] = {
        { { { -0.1666666666666667, 0.1666666666666667, 0.8333333333333333 },
            { -0.6666666666666667, 0.6666666666666667, 0.3333333333333333 },
            { -0.1666666666666666, 0.1666666666666666, 0.8333333333333335 } } }
    };
    static const double FE4_C4_Q48e[1][1][3][3] = {
        { { { 0.1666666666666667, 0.8333333333333333, 0.1666666666666667 },
            { 0.1666666666666666, 0.8333333333333335, 0.1666666666666666 },
            { 0.6666666666666667, 0.3333333333333333, 0.6666666666666667 } } }
    };
    static const double FE6_C0_D10_Q48e[1][1][1][3] = { { { { -1.0, 1.0, 0.0 } } } };
    static const double FE6_C1_D01_Q48e[1][1][1][3] = { { { { -1.0, 0.0, 1.0 } } } };
    // ------------------------
    // Section: Jacobian
    // Inputs: FE6_C0_D10_Q48e, FE6_C1_D01_Q48e, coordinate_dofs
    // Outputs: J_c2, J_c3, J_c0, J_c1
    double J_c0 = 0.0;
    double J_c3 = 0.0;
    double J_c1 = 0.0;
    double J_c2 = 0.0;
    {
        for ( int ic = 0; ic < 3; ++ic ) {
            J_c0 += coordinate_dofs[(ic)*3] * FE6_C0_D10_Q48e[0][0][0][ic];
            J_c3 += coordinate_dofs[(ic)*3 + 1] * FE6_C1_D01_Q48e[0][0][0][ic];
            J_c1 += coordinate_dofs[(ic)*3] * FE6_C1_D01_Q48e[0][0][0][ic];
            J_c2 += coordinate_dofs[(ic)*3 + 1] * FE6_C0_D10_Q48e[0][0][0][ic];
        }
    }
    // ------------------------
    double sp_48e_0 = J_c0 * J_c3;
    double sp_48e_1 = J_c1 * J_c2;
    double sp_48e_2 = -sp_48e_1;
    double sp_48e_3 = sp_48e_0 + sp_48e_2;
    double sp_48e_4 = J_c0 / sp_48e_3;
    double sp_48e_5 = -J_c1;
    double sp_48e_6 = sp_48e_5 / sp_48e_3;
    double sp_48e_7 = sp_48e_4 * sp_48e_4;
    double sp_48e_8 = sp_48e_4 * sp_48e_6;
    double sp_48e_9 = sp_48e_6 * sp_48e_6;
    double sp_48e_10 = sp_48e_7 + sp_48e_7;
    double sp_48e_11 = sp_48e_8 + sp_48e_8;
    double sp_48e_12 = sp_48e_9 + sp_48e_9;
    double sp_48e_13 = J_c3 / sp_48e_3;
    double sp_48e_14 = -J_c2;
    double sp_48e_15 = sp_48e_14 / sp_48e_3;
    double sp_48e_16 = sp_48e_15 * sp_48e_15;
    double sp_48e_17 = sp_48e_13 * sp_48e_15;
    double sp_48e_18 = sp_48e_13 * sp_48e_13;
    double sp_48e_19 = sp_48e_16 + sp_48e_16;
    double sp_48e_20 = sp_48e_17 + sp_48e_17;
    double sp_48e_21 = sp_48e_18 + sp_48e_18;
    double sp_48e_22 = sp_48e_10 + sp_48e_19;
    double sp_48e_23 = sp_48e_11 + sp_48e_20;
    double sp_48e_24 = sp_48e_21 + sp_48e_12;
    double sp_48e_25 = c[0] * c[2];
    double sp_48e_26 = c[3] * sp_48e_25;
    double sp_48e_27 = 1.0 + c[1];
    double sp_48e_28 = 4.0 * sp_48e_27;
    double sp_48e_29 = sp_48e_26 / sp_48e_28;
    double sp_48e_30 = sp_48e_22 * sp_48e_29;
    double sp_48e_31 = sp_48e_23 * sp_48e_29;
    double sp_48e_32 = sp_48e_24 * sp_48e_29;
    double sp_48e_33 = sp_48e_4 + sp_48e_4;
    double sp_48e_34 = sp_48e_6 + sp_48e_6;
    double sp_48e_35 = sp_48e_33 / 2;
    double sp_48e_36 = sp_48e_34 / 2;
    double sp_48e_37 = sp_48e_35 * sp_48e_35;
    double sp_48e_38 = sp_48e_35 * sp_48e_36;
    double sp_48e_39 = sp_48e_36 * sp_48e_36;
    double sp_48e_40 = sp_48e_37 + sp_48e_37;
    double sp_48e_41 = sp_48e_38 + sp_48e_38;
    double sp_48e_42 = sp_48e_39 + sp_48e_39;
    double sp_48e_43 = sp_48e_15 / 2;
    double sp_48e_44 = sp_48e_13 / 2;
    double sp_48e_45 = sp_48e_4 / 2;
    double sp_48e_46 = sp_48e_6 / 2;
    double sp_48e_47 = sp_48e_43 * sp_48e_43;
    double sp_48e_48 = sp_48e_44 * sp_48e_43;
    double sp_48e_49 = sp_48e_45 * sp_48e_43;
    double sp_48e_50 = sp_48e_46 * sp_48e_43;
    double sp_48e_51 = sp_48e_44 * sp_48e_44;
    double sp_48e_52 = sp_48e_45 * sp_48e_44;
    double sp_48e_53 = sp_48e_44 * sp_48e_46;
    double sp_48e_54 = sp_48e_45 * sp_48e_45;
    double sp_48e_55 = sp_48e_45 * sp_48e_46;
    double sp_48e_56 = sp_48e_46 * sp_48e_46;
    double sp_48e_57 = sp_48e_47 + sp_48e_47;
    double sp_48e_58 = sp_48e_48 + sp_48e_48;
    double sp_48e_59 = sp_48e_49 + sp_48e_49;
    double sp_48e_60 = sp_48e_50 + sp_48e_50;
    double sp_48e_61 = sp_48e_51 + sp_48e_51;
    double sp_48e_62 = sp_48e_52 + sp_48e_52;
    double sp_48e_63 = sp_48e_53 + sp_48e_53;
    double sp_48e_64 = sp_48e_54 + sp_48e_54;
    double sp_48e_65 = sp_48e_55 + sp_48e_55;
    double sp_48e_66 = sp_48e_56 + sp_48e_56;
    double sp_48e_67 = sp_48e_40 + sp_48e_57;
    double sp_48e_68 = sp_48e_41 + sp_48e_58;
    double sp_48e_69 = sp_48e_42 + sp_48e_61;
    double sp_48e_70 = sp_48e_15 + sp_48e_15;
    double sp_48e_71 = sp_48e_13 + sp_48e_13;
    double sp_48e_72 = sp_48e_70 / 2;
    double sp_48e_73 = sp_48e_71 / 2;
    double sp_48e_74 = sp_48e_72 * sp_48e_72;
    double sp_48e_75 = sp_48e_73 * sp_48e_72;
    double sp_48e_76 = sp_48e_73 * sp_48e_73;
    double sp_48e_77 = sp_48e_74 + sp_48e_74;
    double sp_48e_78 = sp_48e_75 + sp_48e_75;
    double sp_48e_79 = sp_48e_76 + sp_48e_76;
    double sp_48e_80 = sp_48e_77 + sp_48e_64;
    double sp_48e_81 = sp_48e_78 + sp_48e_65;
    double sp_48e_82 = sp_48e_79 + sp_48e_66;
    double sp_48e_83 = sp_48e_67 + sp_48e_57;
    double sp_48e_84 = sp_48e_68 + sp_48e_58;
    double sp_48e_85 = sp_48e_59 + sp_48e_59;
    double sp_48e_86 = sp_48e_60 + sp_48e_60;
    double sp_48e_87 = sp_48e_69 + sp_48e_61;
    double sp_48e_88 = sp_48e_62 + sp_48e_62;
    double sp_48e_89 = sp_48e_63 + sp_48e_63;
    double sp_48e_90 = sp_48e_80 + sp_48e_64;
    double sp_48e_91 = sp_48e_81 + sp_48e_65;
    double sp_48e_92 = sp_48e_82 + sp_48e_66;
    double sp_48e_93 = -c[1];
    double sp_48e_94 = 1.0 + sp_48e_93;
    double sp_48e_95 = sp_48e_83 * sp_48e_94;
    double sp_48e_96 = sp_48e_84 * sp_48e_94;
    double sp_48e_97 = sp_48e_85 * sp_48e_94;
    double sp_48e_98 = sp_48e_86 * sp_48e_94;
    double sp_48e_99 = sp_48e_87 * sp_48e_94;
    double sp_48e_100 = sp_48e_88 * sp_48e_94;
    double sp_48e_101 = sp_48e_89 * sp_48e_94;
    double sp_48e_102 = sp_48e_90 * sp_48e_94;
    double sp_48e_103 = sp_48e_91 * sp_48e_94;
    double sp_48e_104 = sp_48e_92 * sp_48e_94;
    double sp_48e_105 = 2 * sp_48e_35;
    double sp_48e_106 = 2 * sp_48e_36;
    double sp_48e_107 = 2 * sp_48e_72;
    double sp_48e_108 = 2 * sp_48e_73;
    double sp_48e_109 = sp_48e_105 * sp_48e_35;
    double sp_48e_110 = sp_48e_106 * sp_48e_35;
    double sp_48e_111 = sp_48e_107 * sp_48e_35;
    double sp_48e_112 = sp_48e_108 * sp_48e_35;
    double sp_48e_113 = sp_48e_105 * sp_48e_36;
    double sp_48e_114 = sp_48e_106 * sp_48e_36;
    double sp_48e_115 = sp_48e_107 * sp_48e_36;
    double sp_48e_116 = sp_48e_108 * sp_48e_36;
    double sp_48e_117 = sp_48e_105 * sp_48e_72;
    double sp_48e_118 = sp_48e_106 * sp_48e_72;
    double sp_48e_119 = sp_48e_107 * sp_48e_72;
    double sp_48e_120 = sp_48e_108 * sp_48e_72;
    double sp_48e_121 = sp_48e_105 * sp_48e_73;
    double sp_48e_122 = sp_48e_106 * sp_48e_73;
    double sp_48e_123 = sp_48e_107 * sp_48e_73;
    double sp_48e_124 = sp_48e_108 * sp_48e_73;
    double sp_48e_125 = c[1] * sp_48e_109;
    double sp_48e_126 = c[1] * sp_48e_113;
    double sp_48e_127 = c[1] * sp_48e_117;
    double sp_48e_128 = c[1] * sp_48e_121;
    double sp_48e_129 = c[1] * sp_48e_110;
    double sp_48e_130 = c[1] * sp_48e_114;
    double sp_48e_131 = c[1] * sp_48e_118;
    double sp_48e_132 = c[1] * sp_48e_122;
    double sp_48e_133 = c[1] * sp_48e_111;
    double sp_48e_134 = c[1] * sp_48e_112;
    double sp_48e_135 = c[1] * sp_48e_115;
    double sp_48e_136 = c[1] * sp_48e_116;
    double sp_48e_137 = c[1] * sp_48e_119;
    double sp_48e_138 = c[1] * sp_48e_123;
    double sp_48e_139 = c[1] * sp_48e_120;
    double sp_48e_140 = c[1] * sp_48e_124;
    double sp_48e_141 = sp_48e_95 + sp_48e_125;
    double sp_48e_142 = sp_48e_96 + sp_48e_126;
    double sp_48e_143 = sp_48e_97 + sp_48e_127;
    double sp_48e_144 = sp_48e_98 + sp_48e_128;
    double sp_48e_145 = sp_48e_96 + sp_48e_129;
    double sp_48e_146 = sp_48e_99 + sp_48e_130;
    double sp_48e_147 = sp_48e_100 + sp_48e_131;
    double sp_48e_148 = sp_48e_101 + sp_48e_132;
    double sp_48e_149 = sp_48e_97 + sp_48e_133;
    double sp_48e_150 = sp_48e_98 + sp_48e_134;
    double sp_48e_151 = sp_48e_100 + sp_48e_135;
    double sp_48e_152 = sp_48e_101 + sp_48e_136;
    double sp_48e_153 = sp_48e_102 + sp_48e_137;
    double sp_48e_154 = sp_48e_103 + sp_48e_138;
    double sp_48e_155 = sp_48e_103 + sp_48e_139;
    double sp_48e_156 = sp_48e_104 + sp_48e_140;
    double sp_48e_157 = pow( c[3], 3 );
    double sp_48e_158 = c[0] * sp_48e_157;
    double sp_48e_159 = pow( c[1], 2 );
    double sp_48e_160 = -sp_48e_159;
    double sp_48e_161 = 1.0 + sp_48e_160;
    double sp_48e_162 = 24.0 * sp_48e_161;
    double sp_48e_163 = sp_48e_158 / sp_48e_162;
    double sp_48e_164 = sp_48e_141 * sp_48e_163;
    double sp_48e_165 = sp_48e_142 * sp_48e_163;
    double sp_48e_166 = sp_48e_143 * sp_48e_163;
    double sp_48e_167 = sp_48e_144 * sp_48e_163;
    double sp_48e_168 = sp_48e_145 * sp_48e_163;
    double sp_48e_169 = sp_48e_146 * sp_48e_163;
    double sp_48e_170 = sp_48e_147 * sp_48e_163;
    double sp_48e_171 = sp_48e_148 * sp_48e_163;
    double sp_48e_172 = sp_48e_149 * sp_48e_163;
    double sp_48e_173 = sp_48e_150 * sp_48e_163;
    double sp_48e_174 = sp_48e_151 * sp_48e_163;
    double sp_48e_175 = sp_48e_152 * sp_48e_163;
    double sp_48e_176 = sp_48e_153 * sp_48e_163;
    double sp_48e_177 = sp_48e_154 * sp_48e_163;
    double sp_48e_178 = sp_48e_155 * sp_48e_163;
    double sp_48e_179 = sp_48e_156 * sp_48e_163;
    double sp_48e_180 = fabs( sp_48e_3 );
    double sp_48e_181 = sp_48e_30 * sp_48e_180;
    double sp_48e_182 = sp_48e_31 * sp_48e_180;
    double sp_48e_183 = sp_48e_32 * sp_48e_180;
    double sp_48e_184 = sp_48e_164 * sp_48e_180;
    double sp_48e_185 = sp_48e_165 * sp_48e_180;
    double sp_48e_186 = sp_48e_166 * sp_48e_180;
    double sp_48e_187 = sp_48e_167 * sp_48e_180;
    double sp_48e_188 = sp_48e_168 * sp_48e_180;
    double sp_48e_189 = sp_48e_169 * sp_48e_180;
    double sp_48e_190 = sp_48e_170 * sp_48e_180;
    double sp_48e_191 = sp_48e_171 * sp_48e_180;
    double sp_48e_192 = sp_48e_172 * sp_48e_180;
    double sp_48e_193 = sp_48e_173 * sp_48e_180;
    double sp_48e_194 = sp_48e_174 * sp_48e_180;
    double sp_48e_195 = sp_48e_175 * sp_48e_180;
    double sp_48e_196 = sp_48e_176 * sp_48e_180;
    double sp_48e_197 = sp_48e_177 * sp_48e_180;
    double sp_48e_198 = sp_48e_178 * sp_48e_180;
    double sp_48e_199 = sp_48e_179 * sp_48e_180;
    for ( int iq = 0; iq < 3; ++iq ) {
        // ------------------------
        // Section: Intermediates
        // Inputs:
        // Outputs: fw0, fw1, fw2, fw3, fw4, fw5, fw6, fw7, fw8, fw9, fw10, fw11, fw12, fw13, fw14,
        // fw15, fw16, fw17, fw18
        double fw0 = 0;
        double fw1 = 0;
        double fw2 = 0;
        double fw3 = 0;
        double fw4 = 0;
        double fw5 = 0;
        double fw6 = 0;
        double fw7 = 0;
        double fw8 = 0;
        double fw9 = 0;
        double fw10 = 0;
        double fw11 = 0;
        double fw12 = 0;
        double fw13 = 0;
        double fw14 = 0;
        double fw15 = 0;
        double fw16 = 0;
        double fw17 = 0;
        double fw18 = 0;
        {
            fw0 = sp_48e_199 * weights_48e[iq];
            fw1 = sp_48e_198 * weights_48e[iq];
            fw2 = sp_48e_197 * weights_48e[iq];
            fw3 = sp_48e_196 * weights_48e[iq];
            fw4 = sp_48e_195 * weights_48e[iq];
            fw5 = sp_48e_193 * weights_48e[iq];
            fw6 = sp_48e_194 * weights_48e[iq];
            fw7 = sp_48e_192 * weights_48e[iq];
            fw8 = sp_48e_191 * weights_48e[iq];
            fw9 = sp_48e_190 * weights_48e[iq];
            fw10 = sp_48e_187 * weights_48e[iq];
            fw11 = sp_48e_186 * weights_48e[iq];
            fw12 = sp_48e_189 * weights_48e[iq];
            fw13 = sp_48e_188 * weights_48e[iq];
            fw14 = sp_48e_185 * weights_48e[iq];
            fw15 = sp_48e_184 * weights_48e[iq];
            fw16 = sp_48e_183 * weights_48e[iq];
            fw17 = sp_48e_182 * weights_48e[iq];
            fw18 = sp_48e_181 * weights_48e[iq];
        }
        // ------------------------
        // ------------------------
        // Section: Tensor Computation
        // Inputs: fw3, fw12, fw2, fw9, fw4, fw14, fw13, fw6, fw8, fw1, fw7, fw11, fw15,
        // FE4_C1_D01_Q48e, fw5, fw10, fw0, FE4_C1_D10_Q48e Outputs: A
        {
            double temp_0[6] = { 0 };
            for ( int j = 0; j < 6; ++j ) {
                temp_0[j] = fw0 * FE4_C1_D10_Q48e[0][0][iq][j];
            }
            double temp_1[6] = { 0 };
            for ( int j = 0; j < 6; ++j ) {
                temp_1[j] = fw1 * FE4_C1_D01_Q48e[0][0][iq][j];
            }
            double temp_2[6] = { 0 };
            for ( int j = 0; j < 6; ++j ) {
                temp_2[j] = fw2 * FE4_C1_D10_Q48e[0][0][iq][j];
            }
            double temp_3[6] = { 0 };
            for ( int j = 0; j < 6; ++j ) {
                temp_3[j] = fw3 * FE4_C1_D01_Q48e[0][0][iq][j];
            }
            double temp_4[6] = { 0 };
            for ( int j = 0; j < 6; ++j ) {
                temp_4[j] = fw4 * FE4_C1_D10_Q48e[0][0][iq][j];
            }
            double temp_5[6] = { 0 };
            for ( int j = 0; j < 6; ++j ) {
                temp_5[j] = fw5 * FE4_C1_D01_Q48e[0][0][iq][j];
            }
            double temp_6[6] = { 0 };
            for ( int j = 0; j < 6; ++j ) {
                temp_6[j] = fw6 * FE4_C1_D10_Q48e[0][0][iq][j];
            }
            double temp_7[6] = { 0 };
            for ( int j = 0; j < 6; ++j ) {
                temp_7[j] = fw7 * FE4_C1_D01_Q48e[0][0][iq][j];
            }
            double temp_8[6] = { 0 };
            for ( int j = 0; j < 6; ++j ) {
                temp_8[j] = fw8 * FE4_C1_D10_Q48e[0][0][iq][j];
            }
            double temp_9[6] = { 0 };
            for ( int j = 0; j < 6; ++j ) {
                temp_9[j] = fw9 * FE4_C1_D01_Q48e[0][0][iq][j];
            }
            double temp_10[6] = { 0 };
            for ( int j = 0; j < 6; ++j ) {
                temp_10[j] = fw10 * FE4_C1_D10_Q48e[0][0][iq][j];
            }
            double temp_11[6] = { 0 };
            for ( int j = 0; j < 6; ++j ) {
                temp_11[j] = fw11 * FE4_C1_D01_Q48e[0][0][iq][j];
            }
            double temp_12[6] = { 0 };
            for ( int j = 0; j < 6; ++j ) {
                temp_12[j] = fw12 * FE4_C1_D10_Q48e[0][0][iq][j];
            }
            double temp_13[6] = { 0 };
            for ( int j = 0; j < 6; ++j ) {
                temp_13[j] = fw13 * FE4_C1_D01_Q48e[0][0][iq][j];
            }
            double temp_14[6] = { 0 };
            for ( int j = 0; j < 6; ++j ) {
                temp_14[j] = fw14 * FE4_C1_D10_Q48e[0][0][iq][j];
            }
            double temp_15[6] = { 0 };
            for ( int j = 0; j < 6; ++j ) {
                temp_15[j] = fw15 * FE4_C1_D01_Q48e[0][0][iq][j];
            }
            for ( int j = 0; j < 6; ++j ) {
                for ( int i = 0; i < 6; ++i ) {
                    A[21 * ( 2 * ( i ) ) + 2 * ( j )] += FE4_C1_D10_Q48e[0][0][iq][i] * temp_0[j];
                    A[21 * ( 2 * ( i ) ) + 2 * ( j )] += FE4_C1_D10_Q48e[0][0][iq][i] * temp_1[j];
                    A[21 * ( 2 * ( i ) ) + 2 * ( j )] += FE4_C1_D01_Q48e[0][0][iq][i] * temp_2[j];
                    A[21 * ( 2 * ( i ) ) + 2 * ( j )] += FE4_C1_D01_Q48e[0][0][iq][i] * temp_3[j];
                    A[21 * ( 2 * ( i ) ) + ( 2 * ( j ) + 1 )] +=
                        FE4_C1_D10_Q48e[0][0][iq][i] * temp_4[j];
                    A[21 * ( 2 * ( i ) ) + ( 2 * ( j ) + 1 )] +=
                        FE4_C1_D10_Q48e[0][0][iq][i] * temp_5[j];
                    A[21 * ( 2 * ( i ) ) + ( 2 * ( j ) + 1 )] +=
                        FE4_C1_D01_Q48e[0][0][iq][i] * temp_6[j];
                    A[21 * ( 2 * ( i ) ) + ( 2 * ( j ) + 1 )] +=
                        FE4_C1_D01_Q48e[0][0][iq][i] * temp_7[j];
                    A[21 * ( 2 * ( i ) + 1 ) + 2 * ( j )] +=
                        FE4_C1_D10_Q48e[0][0][iq][i] * temp_8[j];
                    A[21 * ( 2 * ( i ) + 1 ) + 2 * ( j )] +=
                        FE4_C1_D10_Q48e[0][0][iq][i] * temp_9[j];
                    A[21 * ( 2 * ( i ) + 1 ) + 2 * ( j )] +=
                        FE4_C1_D01_Q48e[0][0][iq][i] * temp_10[j];
                    A[21 * ( 2 * ( i ) + 1 ) + 2 * ( j )] +=
                        FE4_C1_D01_Q48e[0][0][iq][i] * temp_11[j];
                    A[21 * ( 2 * ( i ) + 1 ) + ( 2 * ( j ) + 1 )] +=
                        FE4_C1_D10_Q48e[0][0][iq][i] * temp_12[j];
                    A[21 * ( 2 * ( i ) + 1 ) + ( 2 * ( j ) + 1 )] +=
                        FE4_C1_D10_Q48e[0][0][iq][i] * temp_13[j];
                    A[21 * ( 2 * ( i ) + 1 ) + ( 2 * ( j ) + 1 )] +=
                        FE4_C1_D01_Q48e[0][0][iq][i] * temp_14[j];
                    A[21 * ( 2 * ( i ) + 1 ) + ( 2 * ( j ) + 1 )] +=
                        FE4_C1_D01_Q48e[0][0][iq][i] * temp_15[j];
                }
            }
        }
        // ------------------------
        // ------------------------
        // Section: Tensor Computation
        // Inputs: FE4_C4_Q48e, fw16, FE4_C3_Q48e, fw17, fw18
        // Outputs: A
        {
            double temp_0[3] = { 0 };
            for ( int j = 0; j < 3; ++j ) {
                temp_0[j] = fw16 * FE4_C3_Q48e[0][0][iq][j];
            }
            double temp_1[3] = { 0 };
            for ( int j = 0; j < 3; ++j ) {
                temp_1[j] = fw17 * FE4_C4_Q48e[0][0][iq][j];
            }
            double temp_2[3] = { 0 };
            for ( int j = 0; j < 3; ++j ) {
                temp_2[j] = fw17 * FE4_C3_Q48e[0][0][iq][j];
            }
            double temp_3[3] = { 0 };
            for ( int j = 0; j < 3; ++j ) {
                temp_3[j] = fw18 * FE4_C4_Q48e[0][0][iq][j];
            }
            for ( int j = 0; j < 3; ++j ) {
                for ( int i = 0; i < 3; ++i ) {
                    A[21 * ( ( i ) + 15 ) + ( ( j ) + 15 )] += FE4_C3_Q48e[0][0][iq][i] * temp_0[j];
                    A[21 * ( ( i ) + 15 ) + ( ( j ) + 15 )] += FE4_C3_Q48e[0][0][iq][i] * temp_1[j];
                    A[21 * ( ( i ) + 15 ) + ( ( j ) + 15 )] += FE4_C4_Q48e[0][0][iq][i] * temp_2[j];
                    A[21 * ( ( i ) + 15 ) + ( ( j ) + 15 )] += FE4_C4_Q48e[0][0][iq][i] * temp_3[j];
                }
            }
        }
        // ------------------------
    }
    return A;
}

VectorReal B_p5_tr6( const VectorReal &w, const VectorReal &coordinate_dofs,
                     const VectorInt &entity_local_index, const VectorReal &c ) {
    VectorReal A( 42 * 42, 0.0 ); // Vecteur 1D pour stocker les valeurs (size = nbddl)

    // Quadrature rules
    static const double weights_4a8[2] = { 0.5, 0.5 };
    // Precomputed values of basis functions and precomputations
    // FE* dimensions: [permutation][entities][points][dofs]
    static const double FE4_C0_F_Q4a8[1][3][2][6] = {
        { { { 0.0, 0.4553418012614797, -0.1220084679281462, 0.6666666666666667, 0.0, 0.0 },
            { 0.0, -0.1220084679281461, 0.4553418012614795, 0.6666666666666667, 0.0, 0.0 } },
          { { 0.4553418012614797, 0.0, -0.1220084679281462, 0.0, 0.6666666666666669, 0.0 },
            { -0.1220084679281461, 0.0, 0.4553418012614795, 0.0, 0.6666666666666667, 0.0 } },
          { { 0.4553418012614794, -0.1220084679281462, 0.0, 0.0, 0.0, 0.6666666666666665 },
            { -0.1220084679281462, 0.4553418012614795, 0.0, 0.0, 0.0, 0.6666666666666665 } } }
    };
    static const double FE4_C2_D10_F_Q4a8[1][1][1][3] = { { { { -1.0, 1.0, 0.0 } } } };
    static const double FE4_C3_F_Q4a8[1][3][2][3] = {
        { { { -0.211324865405187, 0.211324865405187, 0.788675134594813 },
            { -0.7886751345948131, 0.7886751345948131, 0.211324865405187 } },
          { { -0.2113248654051873, 0.2113248654051873, 0.7886751345948128 },
            { -0.7886751345948131, 0.7886751345948131, 0.2113248654051869 } },
          { { 0.0, 0.0, 1.0 }, { 0.0, 0.0, 1.0 } } }
    };
    static const double FE4_C4_F_Q4a8[1][3][2][3] = {
        { { { 0.788675134594813, 0.211324865405187, 0.788675134594813 },
            { 0.211324865405187, 0.788675134594813, 0.211324865405187 } },
          { { 0.0, 1.0, 0.0 }, { 0.0, 1.0, 0.0 } },
          { { 0.2113248654051872, 0.7886751345948128, 0.2113248654051872 },
            { 0.788675134594813, 0.211324865405187, 0.788675134594813 } } }
    };
    static const double FE6_C1_D01_F_Q4a8[1][1][1][3] = { { { { -1.0, 0.0, 1.0 } } } };
    static const double triangle_reference_facet_jacobian[3][2][1] = { { { -1.0 }, { 1.0 } },
                                                                       { { 0.0 }, { 1.0 } },
                                                                       { { 1.0 }, { 0.0 } } };
    static const double triangle_reference_facet_normals[3][2] = {
        { 0.7071067811865475, 0.7071067811865475 }, { -1.0, -0.0 }, { 0.0, -1.0 }
    };
    // ------------------------
    // Section: Jacobian
    // Inputs: FE6_C1_D01_F_Q4a8, FE4_C2_D10_F_Q4a8, coordinate_dofs
    // Outputs: J_c2, J_c3, J_c0, J_c1
    double J_c3 = 0.0;
    double J_c0 = 0.0;
    double J_c1 = 0.0;
    double J_c2 = 0.0;
    {
        for ( int ic = 0; ic < 3; ++ic ) {
            J_c3 += coordinate_dofs[(ic)*3 + 1] * FE6_C1_D01_F_Q4a8[0][0][0][ic];
            J_c0 += coordinate_dofs[(ic)*3] * FE4_C2_D10_F_Q4a8[0][0][0][ic];
            J_c1 += coordinate_dofs[(ic)*3] * FE6_C1_D01_F_Q4a8[0][0][0][ic];
            J_c2 += coordinate_dofs[(ic)*3 + 1] * FE4_C2_D10_F_Q4a8[0][0][0][ic];
        }
    }
    // ------------------------
    double sp_4a8_0 = J_c0 * J_c3;
    double sp_4a8_1 = J_c1 * J_c2;
    double sp_4a8_2 = -sp_4a8_1;
    double sp_4a8_3 = sp_4a8_0 + sp_4a8_2;
    double sp_4a8_4 = J_c3 / sp_4a8_3;
    double sp_4a8_5 = -J_c2;
    double sp_4a8_6 = sp_4a8_5 / sp_4a8_3;
    double sp_4a8_7 = -sp_4a8_4;
    double sp_4a8_8 = -sp_4a8_6;
    double sp_4a8_9 = J_c0 / sp_4a8_3;
    double sp_4a8_10 = triangle_reference_facet_normals[entity_local_index[0]][1] * sp_4a8_9;
    double sp_4a8_11 = -J_c1;
    double sp_4a8_12 = sp_4a8_11 / sp_4a8_3;
    double sp_4a8_13 = triangle_reference_facet_normals[entity_local_index[0]][0] * sp_4a8_12;
    double sp_4a8_14 = sp_4a8_10 + sp_4a8_13;
    double sp_4a8_15 = sp_4a8_14 * sp_4a8_14;
    double sp_4a8_16 = triangle_reference_facet_normals[entity_local_index[0]][0] * sp_4a8_4;
    double sp_4a8_17 = triangle_reference_facet_normals[entity_local_index[0]][1] * sp_4a8_6;
    double sp_4a8_18 = sp_4a8_16 + sp_4a8_17;
    double sp_4a8_19 = sp_4a8_18 * sp_4a8_18;
    double sp_4a8_20 = sp_4a8_15 + sp_4a8_19;
    double sp_4a8_21 = sqrt( sp_4a8_20 );
    double sp_4a8_22 = sp_4a8_14 / sp_4a8_21;
    double sp_4a8_23 = -sp_4a8_22;
    double sp_4a8_24 = sp_4a8_23 * sp_4a8_4;
    double sp_4a8_25 = sp_4a8_23 * sp_4a8_6;
    double sp_4a8_26 = -sp_4a8_23;
    double sp_4a8_27 = sp_4a8_7 * sp_4a8_23;
    double sp_4a8_28 = sp_4a8_8 * sp_4a8_23;
    double sp_4a8_29 = -sp_4a8_12;
    double sp_4a8_30 = -sp_4a8_9;
    double sp_4a8_31 = sp_4a8_18 / sp_4a8_21;
    double sp_4a8_32 = sp_4a8_12 * sp_4a8_31;
    double sp_4a8_33 = sp_4a8_9 * sp_4a8_31;
    double sp_4a8_34 = sp_4a8_29 * sp_4a8_31;
    double sp_4a8_35 = sp_4a8_30 * sp_4a8_31;
    double sp_4a8_36 = -sp_4a8_31;
    double sp_4a8_37 = sp_4a8_24 + sp_4a8_32;
    double sp_4a8_38 = sp_4a8_25 + sp_4a8_33;
    double sp_4a8_39 = sp_4a8_27 + sp_4a8_34;
    double sp_4a8_40 = sp_4a8_28 + sp_4a8_35;
    double sp_4a8_41 = sp_4a8_37 * sp_4a8_37;
    double sp_4a8_42 = sp_4a8_38 * sp_4a8_37;
    double sp_4a8_43 = sp_4a8_38 * sp_4a8_38;
    double sp_4a8_44 = sp_4a8_37 * sp_4a8_26;
    double sp_4a8_45 = sp_4a8_38 * sp_4a8_26;
    double sp_4a8_46 = sp_4a8_39 * sp_4a8_37;
    double sp_4a8_47 = sp_4a8_39 * sp_4a8_38;
    double sp_4a8_48 = sp_4a8_40 * sp_4a8_37;
    double sp_4a8_49 = sp_4a8_40 * sp_4a8_38;
    double sp_4a8_50 = sp_4a8_37 * sp_4a8_36;
    double sp_4a8_51 = sp_4a8_38 * sp_4a8_36;
    double sp_4a8_52 = J_c0 * triangle_reference_facet_jacobian[entity_local_index[0]][0][0];
    double sp_4a8_53 = J_c1 * triangle_reference_facet_jacobian[entity_local_index[0]][1][0];
    double sp_4a8_54 = sp_4a8_52 + sp_4a8_53;
    double sp_4a8_55 = sp_4a8_54 * sp_4a8_54;
    double sp_4a8_56 = triangle_reference_facet_jacobian[entity_local_index[0]][0][0] * J_c2;
    double sp_4a8_57 = triangle_reference_facet_jacobian[entity_local_index[0]][1][0] * J_c3;
    double sp_4a8_58 = sp_4a8_56 + sp_4a8_57;
    double sp_4a8_59 = sp_4a8_58 * sp_4a8_58;
    double sp_4a8_60 = sp_4a8_55 + sp_4a8_59;
    double sp_4a8_61 = sqrt( sp_4a8_60 );
    double sp_4a8_62 = sp_4a8_41 * sp_4a8_61;
    double sp_4a8_63 = sp_4a8_42 * sp_4a8_61;
    double sp_4a8_64 = sp_4a8_43 * sp_4a8_61;
    double sp_4a8_65 = sp_4a8_44 * sp_4a8_61;
    double sp_4a8_66 = sp_4a8_45 * sp_4a8_61;
    double sp_4a8_67 = sp_4a8_46 * sp_4a8_61;
    double sp_4a8_68 = sp_4a8_47 * sp_4a8_61;
    double sp_4a8_69 = sp_4a8_48 * sp_4a8_61;
    double sp_4a8_70 = sp_4a8_49 * sp_4a8_61;
    double sp_4a8_71 = sp_4a8_50 * sp_4a8_61;
    double sp_4a8_72 = sp_4a8_51 * sp_4a8_61;
    for ( int iq = 0; iq < 2; ++iq ) {
        // ------------------------
        // Section: Intermediates
        // Inputs:
        // Outputs: fw0, fw1, fw2, fw3, fw4, fw5, fw6, fw7, fw8, fw9, fw10
        double fw0 = 0;
        double fw1 = 0;
        double fw2 = 0;
        double fw3 = 0;
        double fw4 = 0;
        double fw5 = 0;
        double fw6 = 0;
        double fw7 = 0;
        double fw8 = 0;
        double fw9 = 0;
        double fw10 = 0;
        {
            fw0 = sp_4a8_65 * weights_4a8[iq];
            fw1 = sp_4a8_66 * weights_4a8[iq];
            fw2 = sp_4a8_71 * weights_4a8[iq];
            fw3 = sp_4a8_72 * weights_4a8[iq];
            fw4 = sp_4a8_62 * weights_4a8[iq];
            fw5 = sp_4a8_63 * weights_4a8[iq];
            fw6 = sp_4a8_64 * weights_4a8[iq];
            fw7 = sp_4a8_67 * weights_4a8[iq];
            fw8 = sp_4a8_68 * weights_4a8[iq];
            fw9 = sp_4a8_69 * weights_4a8[iq];
            fw10 = sp_4a8_70 * weights_4a8[iq];
        }
        // ------------------------
        // ------------------------
        // Section: Tensor Computation
        // Inputs: fw3, fw2, FE4_C0_F_Q4a8, FE4_C4_F_Q4a8, fw1, FE4_C3_F_Q4a8, fw0
        // Outputs: A
        {
            double temp_0[3] = { 0 };
            for ( int j = 0; j < 3; ++j ) {
                temp_0[j] = fw0 * FE4_C3_F_Q4a8[0][entity_local_index[0]][iq][j];
            }
            double temp_1[3] = { 0 };
            for ( int j = 0; j < 3; ++j ) {
                temp_1[j] = fw1 * FE4_C4_F_Q4a8[0][entity_local_index[0]][iq][j];
            }
            double temp_2[3] = { 0 };
            for ( int j = 0; j < 3; ++j ) {
                temp_2[j] = fw2 * FE4_C3_F_Q4a8[0][entity_local_index[0]][iq][j];
            }
            double temp_3[3] = { 0 };
            for ( int j = 0; j < 3; ++j ) {
                temp_3[j] = fw3 * FE4_C4_F_Q4a8[0][entity_local_index[0]][iq][j];
            }
            for ( int j = 0; j < 3; ++j ) {
                for ( int i = 0; i < 6; ++i ) {
                    A[21 * ( 2 * ( i ) ) + ( ( j ) + 18 )] +=
                        FE4_C0_F_Q4a8[0][entity_local_index[0]][iq][i] * temp_0[j];
                    A[21 * ( 2 * ( i ) ) + ( ( j ) + 18 )] +=
                        FE4_C0_F_Q4a8[0][entity_local_index[0]][iq][i] * temp_1[j];
                    A[21 * ( 2 * ( i ) + 1 ) + ( ( j ) + 18 )] +=
                        FE4_C0_F_Q4a8[0][entity_local_index[0]][iq][i] * temp_2[j];
                    A[21 * ( 2 * ( i ) + 1 ) + ( ( j ) + 18 )] +=
                        FE4_C0_F_Q4a8[0][entity_local_index[0]][iq][i] * temp_3[j];
                }
            }
        }
        // ------------------------
        // ------------------------
        // Section: Tensor Computation
        // Inputs: fw9, fw4, FE4_C4_F_Q4a8, fw6, fw8, FE4_C3_F_Q4a8, FE6_C1_D01_F_Q4a8, fw7, fw5,
        // fw10, FE4_C2_D10_F_Q4a8 Outputs: A
        {
            double temp_0[3] = { 0 };
            for ( int j = 0; j < 3; ++j ) {
                temp_0[j] = fw4 * FE4_C3_F_Q4a8[0][entity_local_index[0]][iq][j];
            }
            double temp_1[3] = { 0 };
            for ( int j = 0; j < 3; ++j ) {
                temp_1[j] = fw5 * FE4_C4_F_Q4a8[0][entity_local_index[0]][iq][j];
            }
            double temp_2[3] = { 0 };
            for ( int j = 0; j < 3; ++j ) {
                temp_2[j] = fw5 * FE4_C3_F_Q4a8[0][entity_local_index[0]][iq][j];
            }
            double temp_3[3] = { 0 };
            for ( int j = 0; j < 3; ++j ) {
                temp_3[j] = fw6 * FE4_C4_F_Q4a8[0][entity_local_index[0]][iq][j];
            }
            double temp_4[3] = { 0 };
            for ( int j = 0; j < 3; ++j ) {
                temp_4[j] = fw7 * FE4_C3_F_Q4a8[0][entity_local_index[0]][iq][j];
            }
            double temp_5[3] = { 0 };
            for ( int j = 0; j < 3; ++j ) {
                temp_5[j] = fw8 * FE4_C4_F_Q4a8[0][entity_local_index[0]][iq][j];
            }
            double temp_6[3] = { 0 };
            for ( int j = 0; j < 3; ++j ) {
                temp_6[j] = fw9 * FE4_C3_F_Q4a8[0][entity_local_index[0]][iq][j];
            }
            double temp_7[3] = { 0 };
            for ( int j = 0; j < 3; ++j ) {
                temp_7[j] = fw10 * FE4_C4_F_Q4a8[0][entity_local_index[0]][iq][j];
            }
            double temp_8[3] = { 0 };
            for ( int j = 0; j < 3; ++j ) {
                temp_8[j] = fw4 * FE4_C2_D10_F_Q4a8[0][0][0][j];
            }
            double temp_9[3] = { 0 };
            for ( int j = 0; j < 3; ++j ) {
                temp_9[j] = fw5 * FE6_C1_D01_F_Q4a8[0][0][0][j];
            }
            double temp_10[3] = { 0 };
            for ( int j = 0; j < 3; ++j ) {
                temp_10[j] = fw5 * FE4_C2_D10_F_Q4a8[0][0][0][j];
            }
            double temp_11[3] = { 0 };
            for ( int j = 0; j < 3; ++j ) {
                temp_11[j] = fw6 * FE6_C1_D01_F_Q4a8[0][0][0][j];
            }
            double temp_12[3] = { 0 };
            for ( int j = 0; j < 3; ++j ) {
                temp_12[j] = fw7 * FE4_C3_F_Q4a8[0][entity_local_index[0]][iq][j];
            }
            double temp_13[3] = { 0 };
            for ( int j = 0; j < 3; ++j ) {
                temp_13[j] = fw9 * FE4_C4_F_Q4a8[0][entity_local_index[0]][iq][j];
            }
            double temp_14[3] = { 0 };
            for ( int j = 0; j < 3; ++j ) {
                temp_14[j] = fw8 * FE4_C3_F_Q4a8[0][entity_local_index[0]][iq][j];
            }
            double temp_15[3] = { 0 };
            for ( int j = 0; j < 3; ++j ) {
                temp_15[j] = fw10 * FE4_C4_F_Q4a8[0][entity_local_index[0]][iq][j];
            }
            for ( int j = 0; j < 3; ++j ) {
                for ( int i = 0; i < 3; ++i ) {
                    A[21 * ( ( i ) + 12 ) + ( ( j ) + 18 )] +=
                        FE4_C2_D10_F_Q4a8[0][0][0][i] * temp_0[j];
                    A[21 * ( ( i ) + 12 ) + ( ( j ) + 18 )] +=
                        FE4_C2_D10_F_Q4a8[0][0][0][i] * temp_1[j];
                    A[21 * ( ( i ) + 12 ) + ( ( j ) + 18 )] +=
                        FE6_C1_D01_F_Q4a8[0][0][0][i] * temp_2[j];
                    A[21 * ( ( i ) + 12 ) + ( ( j ) + 18 )] +=
                        FE6_C1_D01_F_Q4a8[0][0][0][i] * temp_3[j];
                    A[21 * ( ( i ) + 15 ) + ( ( j ) + 18 )] +=
                        FE4_C3_F_Q4a8[0][entity_local_index[0]][iq][i] * temp_4[j];
                    A[21 * ( ( i ) + 15 ) + ( ( j ) + 18 )] +=
                        FE4_C3_F_Q4a8[0][entity_local_index[0]][iq][i] * temp_5[j];
                    A[21 * ( ( i ) + 15 ) + ( ( j ) + 18 )] +=
                        FE4_C4_F_Q4a8[0][entity_local_index[0]][iq][i] * temp_6[j];
                    A[21 * ( ( i ) + 15 ) + ( ( j ) + 18 )] +=
                        FE4_C4_F_Q4a8[0][entity_local_index[0]][iq][i] * temp_7[j];
                    A[21 * ( ( i ) + 18 ) + ( ( j ) + 12 )] +=
                        FE4_C3_F_Q4a8[0][entity_local_index[0]][iq][i] * temp_8[j];
                    A[21 * ( ( i ) + 18 ) + ( ( j ) + 12 )] +=
                        FE4_C3_F_Q4a8[0][entity_local_index[0]][iq][i] * temp_9[j];
                    A[21 * ( ( i ) + 18 ) + ( ( j ) + 12 )] +=
                        FE4_C4_F_Q4a8[0][entity_local_index[0]][iq][i] * temp_10[j];
                    A[21 * ( ( i ) + 18 ) + ( ( j ) + 12 )] +=
                        FE4_C4_F_Q4a8[0][entity_local_index[0]][iq][i] * temp_11[j];
                    A[21 * ( ( i ) + 18 ) + ( ( j ) + 15 )] +=
                        FE4_C3_F_Q4a8[0][entity_local_index[0]][iq][i] * temp_12[j];
                    A[21 * ( ( i ) + 18 ) + ( ( j ) + 15 )] +=
                        FE4_C3_F_Q4a8[0][entity_local_index[0]][iq][i] * temp_13[j];
                    A[21 * ( ( i ) + 18 ) + ( ( j ) + 15 )] +=
                        FE4_C4_F_Q4a8[0][entity_local_index[0]][iq][i] * temp_14[j];
                    A[21 * ( ( i ) + 18 ) + ( ( j ) + 15 )] +=
                        FE4_C4_F_Q4a8[0][entity_local_index[0]][iq][i] * temp_15[j];
                }
            }
        }
        // ------------------------
        // ------------------------
        // Section: Tensor Computation
        // Inputs: fw3, fw2, FE4_C0_F_Q4a8, FE4_C4_F_Q4a8, fw1, FE4_C3_F_Q4a8, fw0
        // Outputs: A
        {
            double temp_0[6] = { 0 };
            for ( int j = 0; j < 6; ++j ) {
                temp_0[j] = fw0 * FE4_C0_F_Q4a8[0][entity_local_index[0]][iq][j];
            }
            double temp_1[6] = { 0 };
            for ( int j = 0; j < 6; ++j ) {
                temp_1[j] = fw1 * FE4_C0_F_Q4a8[0][entity_local_index[0]][iq][j];
            }
            double temp_2[6] = { 0 };
            for ( int j = 0; j < 6; ++j ) {
                temp_2[j] = fw2 * FE4_C0_F_Q4a8[0][entity_local_index[0]][iq][j];
            }
            double temp_3[6] = { 0 };
            for ( int j = 0; j < 6; ++j ) {
                temp_3[j] = fw3 * FE4_C0_F_Q4a8[0][entity_local_index[0]][iq][j];
            }
            for ( int j = 0; j < 6; ++j ) {
                for ( int i = 0; i < 3; ++i ) {
                    A[21 * ( ( i ) + 18 ) + 2 * ( j )] +=
                        FE4_C3_F_Q4a8[0][entity_local_index[0]][iq][i] * temp_0[j];
                    A[21 * ( ( i ) + 18 ) + 2 * ( j )] +=
                        FE4_C4_F_Q4a8[0][entity_local_index[0]][iq][i] * temp_1[j];
                    A[21 * ( ( i ) + 18 ) + ( 2 * ( j ) + 1 )] +=
                        FE4_C3_F_Q4a8[0][entity_local_index[0]][iq][i] * temp_2[j];
                    A[21 * ( ( i ) + 18 ) + ( 2 * ( j ) + 1 )] +=
                        FE4_C4_F_Q4a8[0][entity_local_index[0]][iq][i] * temp_3[j];
                }
            }
        }
        // ------------------------
    }

    return A;
}

extern "C" {
void BP1_tr6_Fortran( const double *w, const int nw, const double *coordinate_dofs, const int ncd,
                      const int *entities, const int ne, const double *cst, const int ncst,
                      double *A ) {
    if ( nw <= 0 || ncd <= 0 || ncst <= 0 || A == nullptr ) {
        std::cerr << "Erreur : tailles invalides dans K_elem_fortran." << std::endl;
        return;
    }

    // Création des vecteurs à partir des pointeurs C
    VectorReal w_vec( w, w + nw );
    VectorInt entities_vec( entities, entities + ne );
    VectorReal coordinate_dofs_vec( coordinate_dofs, coordinate_dofs + ncd );
    VectorReal cst_vec( cst, cst + ncst );

    // Appel de la fonction C++ de calcul
    VectorReal A_vec = B_p1_tr6( w_vec, coordinate_dofs_vec, entities_vec, cst_vec );

    // Copie du résultat dans le tableau C fourni par Fortran
    for ( size_t i = 0; i < A_vec.size(); ++i ) {
        A[i] = A_vec[i];
    }
}
void BP2_tr6_Fortran( const double *w, const int nw, const double *coordinate_dofs, const int ncd,
                      const int *entities, const int ne, const double *cst, const int ncst,
                      double *A ) {
    if ( nw <= 0 || ncd <= 0 || ncst <= 0 || A == nullptr ) {
        std::cerr << "Erreur : tailles invalides dans K_elem_fortran." << std::endl;
        return;
    }

    // Création des vecteurs à partir des pointeurs C
    VectorReal w_vec( w, w + nw );
    VectorInt entities_vec( entities, entities + ne );
    VectorReal coordinate_dofs_vec( coordinate_dofs, coordinate_dofs + ncd );
    VectorReal cst_vec( cst, cst + ncst );

    // Appel de la fonction C++ de calcul
    VectorReal A_vec = B_p2_tr6( w_vec, coordinate_dofs_vec, entities_vec, cst_vec );

    // Copie du résultat dans le tableau C fourni par Fortran
    for ( size_t i = 0; i < A_vec.size(); ++i ) {
        A[i] = A_vec[i];
    }
}
void BP4_tr6_Fortran( const double *w, const int nw, const double *coordinate_dofs, const int ncd,
                      const int *entities, const int ne, const double *cst, const int ncst,
                      double *A ) {
    if ( nw <= 0 || ncd <= 0 || ncst <= 0 || A == nullptr ) {
        std::cerr << "Erreur : tailles invalides dans K_elem_fortran." << std::endl;
        return;
    }

    // Création des vecteurs à partir des pointeurs C
    VectorReal w_vec( w, w + nw );
    VectorInt entities_vec( entities, entities + ne );
    VectorReal coordinate_dofs_vec( coordinate_dofs, coordinate_dofs + ncd );
    VectorReal cst_vec( cst, cst + ncst );

    // Appel de la fonction C++ de calcul
    VectorReal A_vec = B_p4_tr6( w_vec, coordinate_dofs_vec, entities_vec, cst_vec );

    // Copie du résultat dans le tableau C fourni par Fortran
    for ( size_t i = 0; i < A_vec.size(); ++i ) {
        A[i] = A_vec[i];
    }
}
void BP5_tr6_Fortran( const double *w, const int nw, const double *coordinate_dofs, const int ncd,
                      const int *entities, const int ne, const double *cst, const int ncst,
                      double *A ) {
    if ( nw <= 0 || ncd <= 0 || ncst <= 0 || A == nullptr ) {
        std::cerr << "Erreur : tailles invalides dans K_elem_fortran." << std::endl;
        return;
    }

    // Création des vecteurs à partir des pointeurs C
    VectorReal w_vec( w, w + nw );
    VectorInt entities_vec( entities, entities + ne );
    VectorReal coordinate_dofs_vec( coordinate_dofs, coordinate_dofs + ncd );
    VectorReal cst_vec( cst, cst + ncst );

    // Appel de la fonction C++ de calcul
    VectorReal A_vec = B_p5_tr6( w_vec, coordinate_dofs_vec, entities_vec, cst_vec );

    // Copie du résultat dans le tableau C fourni par Fortran
    for ( size_t i = 0; i < A_vec.size(); ++i ) {
        A[i] = A_vec[i];
    }
}
}
