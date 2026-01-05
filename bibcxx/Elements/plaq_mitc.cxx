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

#include "plaq_mitc.h"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

#include <math.h>
#include <stdalign.h>
#include <stdlib.h>
#include <string.h>

VectorReal B_p1_qu9( const VectorReal &w, const VectorReal &coordinate_dofs,
                     const VectorInt &entity_local_index, const VectorReal &c ) {
    VectorReal A( 42 * 42, 0.0 ); // Vecteur 1D pour stocker les valeurs (size = nbddl)
    // Quadrature rules
    static const double weights_0df[4] = { 0.25, 0.25, 0.25, 0.25 };
    static const double weights_e1c[9] = { 0.07716049382716043, 0.1234567901234567,
                                           0.07716049382716043, 0.1234567901234567,
                                           0.1975308641975309,  0.1234567901234567,
                                           0.07716049382716043, 0.1234567901234567,
                                           0.07716049382716043 };
    // Precomputed values of basis functions and precomputations
    // FE* dimensions: [permutation][entities][points][dofs]
    static const double FE4_C1_D01_Q0df[1][1][4][9] = {
        { { { -0.9811252243246881, 0.2628917115316043, -0.07044162180172914, 0.01887477567531192,
              -1.436467025586168, 1.051566846126417, -0.2817664872069161, -0.1031336922528346,
              1.539600717839002 },
            { 0.07044162180172903, -0.01887477567531187, 0.9811252243246882, -0.2628917115316043,
              0.1031336922528346, -1.051566846126417, 0.2817664872069161, 1.436467025586168,
              -1.539600717839002 },
            { 0.2628917115316045, -0.9811252243246883, 0.01887477567531164, -0.07044162180172886,
              -1.436467025586168, -0.2817664872069161, 1.051566846126417, -0.1031336922528346,
              1.539600717839002 },
            { -0.01887477567531198, 0.07044162180172914, -0.2628917115316041, 0.981125224324688,
              0.1031336922528346, 0.2817664872069162, -1.051566846126417, 1.436467025586168,
              -1.539600717839002 } } }
    };
    static const double FE4_C1_D10_Q0df[1][1][4][9] = {
        { { { -0.9811252243246882, -0.07044162180172919, 0.2628917115316043, 0.01887477567531193,
              1.051566846126417, -1.436467025586168, -0.1031336922528345, -0.2817664872069161,
              1.539600717839002 },
            { 0.2628917115316044, 0.01887477567531159, -0.9811252243246882, -0.07044162180172903,
              -0.2817664872069162, -1.436467025586168, -0.1031336922528345, 1.051566846126417,
              1.539600717839002 },
            { 0.07044162180172914, 0.9811252243246882, -0.01887477567531187, -0.2628917115316043,
              -1.051566846126417, 0.1031336922528346, 1.436467025586168, 0.2817664872069161,
              -1.539600717839002 },
            { -0.01887477567531198, -0.2628917115316041, 0.07044162180172925, 0.981125224324688,
              0.2817664872069162, 0.1031336922528345, 1.436467025586168, -1.051566846126417,
              -1.539600717839002 } } }
    };
    static const double FE4_C2_Qe1c[1][1][9][4] = {
        { { { 0.7872983346207416, 0.09999999999999995, 0.09999999999999992, 0.01270166537925832 },
            { 0.4436491673103709, 0.05635083268962912, 0.4436491673103709, 0.05635083268962912 },
            { 0.1000000000000001, 0.01270166537925832, 0.7872983346207417, 0.09999999999999989 },
            { 0.4436491673103709, 0.4436491673103709, 0.05635083268962912, 0.05635083268962912 },
            { 0.25, 0.25, 0.25, 0.25 },
            { 0.05635083268962918, 0.05635083268962918, 0.4436491673103709, 0.4436491673103709 },
            { 0.1000000000000001, 0.7872983346207417, 0.01270166537925832, 0.09999999999999992 },
            { 0.05635083268962918, 0.4436491673103709, 0.05635083268962918, 0.4436491673103709 },
            { 0.01270166537925829, 0.1, 0.09999999999999998, 0.7872983346207418 } } }
    };
    static const double FE4_C3_Q0df[1][1][4][4] = {
        { { { 0.7886751345948129, 0.0, 0.0, 0.2113248654051871 },
            { 0.2113248654051872, 0.0, 0.0, 0.7886751345948129 },
            { 0.7886751345948129, 0.0, 0.0, 0.2113248654051871 },
            { 0.2113248654051872, 0.0, 0.0, 0.7886751345948129 } } }
    };
    static const double FE4_C4_Q0df[1][1][4][4] = {
        { { { 0.0, 0.7886751345948129, 0.2113248654051871, 0.0 },
            { 0.0, 0.7886751345948129, 0.2113248654051871, 0.0 },
            { 0.0, 0.2113248654051872, 0.7886751345948129, 0.0 },
            { 0.0, 0.2113248654051872, 0.7886751345948129, 0.0 } } }
    };
    static const double FE5_C0_D10_Q0df[1][1][4][4] = {
        { { { -0.7886751345948126, 0.7886751345948129, -0.2113248654051871, 0.2113248654051871 },
            { -0.2113248654051872, 0.2113248654051871, -0.7886751345948129, 0.7886751345948129 },
            { -0.7886751345948126, 0.7886751345948129, -0.2113248654051871, 0.2113248654051871 },
            { -0.2113248654051872, 0.2113248654051871, -0.7886751345948129, 0.7886751345948129 } } }
    };
    static const double FE5_C0_D10_Qe1c[1][1][9][4] = {
        { { { -0.8872983346207415, 0.8872983346207417, -0.1127016653792582, 0.1127016653792582 },
            { -0.4999999999999999, 0.5, -0.4999999999999999, 0.5 },
            { -0.1127016653792585, 0.1127016653792583, -0.8872983346207417, 0.8872983346207418 },
            { -0.8872983346207415, 0.8872983346207417, -0.1127016653792582, 0.1127016653792582 },
            { -0.4999999999999999, 0.5, -0.4999999999999999, 0.5 },
            { -0.1127016653792585, 0.1127016653792583, -0.8872983346207417, 0.8872983346207418 },
            { -0.8872983346207415, 0.8872983346207417, -0.1127016653792582, 0.1127016653792582 },
            { -0.4999999999999999, 0.5, -0.4999999999999999, 0.5 },
            { -0.1127016653792585, 0.1127016653792583, -0.8872983346207417, 0.8872983346207418 } } }
    };
    static const double FE5_C1_D01_Q0df[1][1][4][4] = {
        { { { -0.7886751345948126, -0.2113248654051871, 0.7886751345948129, 0.2113248654051871 },
            { -0.7886751345948126, -0.2113248654051871, 0.7886751345948129, 0.2113248654051871 },
            { -0.2113248654051872, -0.7886751345948129, 0.2113248654051871, 0.7886751345948129 },
            { -0.2113248654051872, -0.7886751345948129, 0.2113248654051871, 0.7886751345948129 } } }
    };
    static const double FE5_C1_D01_Qe1c[1][1][9][4] = {
        { { { -0.8872983346207415, -0.1127016653792582, 0.8872983346207417, 0.1127016653792582 },
            { -0.8872983346207415, -0.1127016653792582, 0.8872983346207417, 0.1127016653792582 },
            { -0.8872983346207415, -0.1127016653792582, 0.8872983346207417, 0.1127016653792582 },
            { -0.4999999999999999, -0.4999999999999999, 0.5, 0.5 },
            { -0.4999999999999999, -0.4999999999999999, 0.5, 0.5 },
            { -0.4999999999999999, -0.4999999999999999, 0.5, 0.5 },
            { -0.1127016653792585, -0.8872983346207417, 0.1127016653792583, 0.8872983346207418 },
            { -0.1127016653792585, -0.8872983346207417, 0.1127016653792583, 0.8872983346207418 },
            { -0.1127016653792585, -0.8872983346207417, 0.1127016653792583, 0.8872983346207418 } } }
    };
    double sp_0df_0 = c[0] * c[2];
    double sp_0df_1 = c[3] * sp_0df_0;
    double sp_0df_2 = 1.0 + c[1];
    double sp_0df_3 = 4.0 * sp_0df_2;
    double sp_0df_4 = sp_0df_1 / sp_0df_3;
    double sp_0df_5 = -c[1];
    double sp_0df_6 = 1.0 + sp_0df_5;
    double sp_0df_7 = pow( c[3], 3 );
    double sp_0df_8 = c[0] * sp_0df_7;
    double sp_0df_9 = pow( c[1], 2 );
    double sp_0df_10 = -sp_0df_9;
    double sp_0df_11 = 1.0 + sp_0df_10;
    double sp_0df_12 = 24.0 * sp_0df_11;
    double sp_0df_13 = sp_0df_8 / sp_0df_12;
    double sp_e1c_0 = c[4] * sp_0df_7;
    double sp_e1c_1 = -sp_e1c_0;
    for ( int iq = 0; iq < 4; ++iq ) {
        // ------------------------
        // Section: Jacobian
        // Inputs: FE5_C1_D01_Q0df, FE5_C0_D10_Q0df, coordinate_dofs
        // Outputs: J_c3, J_c0, J_c2, J_c1
        double J_c0 = 0.0;
        double J_c3 = 0.0;
        double J_c1 = 0.0;
        double J_c2 = 0.0;
        {
            for ( int ic = 0; ic < 4; ++ic ) {
                J_c0 += coordinate_dofs[(ic)*3] * FE5_C0_D10_Q0df[0][0][iq][ic];
                J_c3 += coordinate_dofs[(ic)*3 + 1] * FE5_C1_D01_Q0df[0][0][iq][ic];
                J_c1 += coordinate_dofs[(ic)*3] * FE5_C1_D01_Q0df[0][0][iq][ic];
                J_c2 += coordinate_dofs[(ic)*3 + 1] * FE5_C0_D10_Q0df[0][0][iq][ic];
            }
        }
        // ------------------------
        // ------------------------
        // Section: Function
        // Inputs: w, FE4_C4_Q0df
        // Outputs: w0_c4
        double w0_c4 = 0.0;
        {
            for ( int ic = 0; ic < 4; ++ic ) {
                w0_c4 += w[( ic ) + 22] * FE4_C4_Q0df[0][0][iq][ic];
            }
        }
        // ------------------------
        // ------------------------
        // Section: Function
        // Inputs: w, FE4_C3_Q0df
        // Outputs: w0_c3
        double w0_c3 = 0.0;
        {
            for ( int ic = 0; ic < 4; ++ic ) {
                w0_c3 += w[( ic ) + 22] * FE4_C3_Q0df[0][0][iq][ic];
            }
        }
        // ------------------------
        // ------------------------
        // Section: Function
        // Inputs: w, FE4_C1_D01_Q0df
        // Outputs: w0_d01_c1
        double w0_d01_c1 = 0.0;
        {
            for ( int ic = 0; ic < 9; ++ic ) {
                w0_d01_c1 += w[(ic)*2 + 1] * FE4_C1_D01_Q0df[0][0][iq][ic];
            }
        }
        // ------------------------
        // ------------------------
        // Section: Function
        // Inputs: w, FE4_C1_D10_Q0df
        // Outputs: w0_d10_c1
        double w0_d10_c1 = 0.0;
        {
            for ( int ic = 0; ic < 9; ++ic ) {
                w0_d10_c1 += w[(ic)*2 + 1] * FE4_C1_D10_Q0df[0][0][iq][ic];
            }
        }
        // ------------------------
        // ------------------------
        // Section: Function
        // Inputs: w, FE4_C1_D01_Q0df
        // Outputs: w0_d01_c0
        double w0_d01_c0 = 0.0;
        {
            for ( int ic = 0; ic < 9; ++ic ) {
                w0_d01_c0 += w[(ic)*2] * FE4_C1_D01_Q0df[0][0][iq][ic];
            }
        }
        // ------------------------
        // ------------------------
        // Section: Function
        // Inputs: w, FE4_C1_D10_Q0df
        // Outputs: w0_d10_c0
        double w0_d10_c0 = 0.0;
        {
            for ( int ic = 0; ic < 9; ++ic ) {
                w0_d10_c0 += w[(ic)*2] * FE4_C1_D10_Q0df[0][0][iq][ic];
            }
        }
        // ------------------------
        // ------------------------
        // Section: Intermediates
        // Inputs: J_c3, J_c0, J_c2, J_c1, w0_c4, w0_c3, w0_d01_c1, w0_d10_c1, w0_d01_c0, w0_d10_c0
        // Outputs: fw0, fw1, fw2, fw3, fw4, fw5
        double fw0 = 0;
        double fw1 = 0;
        double fw2 = 0;
        double fw3 = 0;
        double fw4 = 0;
        double fw5 = 0;
        {
            double sv_0df_0 = J_c0 * J_c3;
            double sv_0df_1 = J_c1 * J_c2;
            double sv_0df_2 = -sv_0df_1;
            double sv_0df_3 = sv_0df_0 + sv_0df_2;
            double sv_0df_4 = J_c0 / sv_0df_3;
            double sv_0df_5 = -J_c1;
            double sv_0df_6 = sv_0df_5 / sv_0df_3;
            double sv_0df_7 = w0_c4 * sv_0df_4;
            double sv_0df_8 = w0_c3 * sv_0df_6;
            double sv_0df_9 = sv_0df_7 + sv_0df_8;
            double sv_0df_10 = sv_0df_9 * sv_0df_4;
            double sv_0df_11 = sv_0df_9 * sv_0df_6;
            double sv_0df_12 = sv_0df_10 + sv_0df_10;
            double sv_0df_13 = sv_0df_11 + sv_0df_11;
            double sv_0df_14 = J_c3 / sv_0df_3;
            double sv_0df_15 = -J_c2;
            double sv_0df_16 = sv_0df_15 / sv_0df_3;
            double sv_0df_17 = w0_c3 * sv_0df_14;
            double sv_0df_18 = w0_c4 * sv_0df_16;
            double sv_0df_19 = sv_0df_17 + sv_0df_18;
            double sv_0df_20 = sv_0df_19 * sv_0df_16;
            double sv_0df_21 = sv_0df_19 * sv_0df_14;
            double sv_0df_22 = sv_0df_20 + sv_0df_20;
            double sv_0df_23 = sv_0df_21 + sv_0df_21;
            double sv_0df_24 = sv_0df_12 + sv_0df_22;
            double sv_0df_25 = sv_0df_23 + sv_0df_13;
            double sv_0df_26 = sv_0df_24 * sp_0df_4;
            double sv_0df_27 = sv_0df_25 * sp_0df_4;
            double sv_0df_28 = sv_0df_4 + sv_0df_4;
            double sv_0df_29 = sv_0df_6 + sv_0df_6;
            double sv_0df_30 = sv_0df_28 / 2;
            double sv_0df_31 = sv_0df_29 / 2;
            double sv_0df_32 = w0_d01_c1 * sv_0df_4;
            double sv_0df_33 = w0_d10_c1 * sv_0df_6;
            double sv_0df_34 = sv_0df_32 + sv_0df_33;
            double sv_0df_35 = sv_0df_34 + sv_0df_34;
            double sv_0df_36 = sv_0df_35 / 2;
            double sv_0df_37 = sv_0df_36 * sv_0df_30;
            double sv_0df_38 = sv_0df_36 * sv_0df_31;
            double sv_0df_39 = sv_0df_37 + sv_0df_37;
            double sv_0df_40 = sv_0df_38 + sv_0df_38;
            double sv_0df_41 = sv_0df_16 / 2;
            double sv_0df_42 = sv_0df_14 / 2;
            double sv_0df_43 = sv_0df_4 / 2;
            double sv_0df_44 = sv_0df_6 / 2;
            double sv_0df_45 = w0_d01_c0 * sv_0df_4;
            double sv_0df_46 = w0_d10_c0 * sv_0df_6;
            double sv_0df_47 = sv_0df_45 + sv_0df_46;
            double sv_0df_48 = w0_d10_c1 * sv_0df_14;
            double sv_0df_49 = w0_d01_c1 * sv_0df_16;
            double sv_0df_50 = sv_0df_48 + sv_0df_49;
            double sv_0df_51 = sv_0df_47 + sv_0df_50;
            double sv_0df_52 = sv_0df_51 / 2;
            double sv_0df_53 = sv_0df_52 * sv_0df_41;
            double sv_0df_54 = sv_0df_52 * sv_0df_42;
            double sv_0df_55 = sv_0df_52 * sv_0df_43;
            double sv_0df_56 = sv_0df_52 * sv_0df_44;
            double sv_0df_57 = sv_0df_53 + sv_0df_53;
            double sv_0df_58 = sv_0df_54 + sv_0df_54;
            double sv_0df_59 = sv_0df_55 + sv_0df_55;
            double sv_0df_60 = sv_0df_56 + sv_0df_56;
            double sv_0df_61 = sv_0df_39 + sv_0df_57;
            double sv_0df_62 = sv_0df_40 + sv_0df_58;
            double sv_0df_63 = sv_0df_16 + sv_0df_16;
            double sv_0df_64 = sv_0df_14 + sv_0df_14;
            double sv_0df_65 = sv_0df_63 / 2;
            double sv_0df_66 = sv_0df_64 / 2;
            double sv_0df_67 = w0_d10_c0 * sv_0df_14;
            double sv_0df_68 = w0_d01_c0 * sv_0df_16;
            double sv_0df_69 = sv_0df_67 + sv_0df_68;
            double sv_0df_70 = sv_0df_69 + sv_0df_69;
            double sv_0df_71 = sv_0df_70 / 2;
            double sv_0df_72 = sv_0df_71 * sv_0df_65;
            double sv_0df_73 = sv_0df_71 * sv_0df_66;
            double sv_0df_74 = sv_0df_72 + sv_0df_72;
            double sv_0df_75 = sv_0df_73 + sv_0df_73;
            double sv_0df_76 = sv_0df_74 + sv_0df_59;
            double sv_0df_77 = sv_0df_75 + sv_0df_60;
            double sv_0df_78 = sv_0df_61 + sv_0df_57;
            double sv_0df_79 = sv_0df_62 + sv_0df_58;
            double sv_0df_80 = sv_0df_76 + sv_0df_59;
            double sv_0df_81 = sv_0df_77 + sv_0df_60;
            double sv_0df_82 = sv_0df_78 * sp_0df_6;
            double sv_0df_83 = sv_0df_79 * sp_0df_6;
            double sv_0df_84 = sv_0df_80 * sp_0df_6;
            double sv_0df_85 = sv_0df_81 * sp_0df_6;
            double sv_0df_86 = sv_0df_36 + sv_0df_71;
            double sv_0df_87 = 2 * sv_0df_30;
            double sv_0df_88 = 2 * sv_0df_31;
            double sv_0df_89 = 2 * sv_0df_65;
            double sv_0df_90 = 2 * sv_0df_66;
            double sv_0df_91 = sv_0df_86 * sv_0df_87;
            double sv_0df_92 = sv_0df_86 * sv_0df_88;
            double sv_0df_93 = sv_0df_86 * sv_0df_89;
            double sv_0df_94 = sv_0df_86 * sv_0df_90;
            double sv_0df_95 = c[1] * sv_0df_91;
            double sv_0df_96 = c[1] * sv_0df_92;
            double sv_0df_97 = c[1] * sv_0df_93;
            double sv_0df_98 = c[1] * sv_0df_94;
            double sv_0df_99 = sv_0df_82 + sv_0df_95;
            double sv_0df_100 = sv_0df_83 + sv_0df_96;
            double sv_0df_101 = sv_0df_84 + sv_0df_97;
            double sv_0df_102 = sv_0df_85 + sv_0df_98;
            double sv_0df_103 = sv_0df_99 * sp_0df_13;
            double sv_0df_104 = sv_0df_100 * sp_0df_13;
            double sv_0df_105 = sv_0df_101 * sp_0df_13;
            double sv_0df_106 = sv_0df_102 * sp_0df_13;
            double sv_0df_107 = fabs( sv_0df_3 );
            double sv_0df_108 = sv_0df_26 * sv_0df_107;
            double sv_0df_109 = sv_0df_27 * sv_0df_107;
            double sv_0df_110 = sv_0df_103 * sv_0df_107;
            double sv_0df_111 = sv_0df_104 * sv_0df_107;
            double sv_0df_112 = sv_0df_105 * sv_0df_107;
            double sv_0df_113 = sv_0df_106 * sv_0df_107;
            fw0 = sv_0df_113 * weights_0df[iq];
            fw1 = sv_0df_112 * weights_0df[iq];
            fw2 = sv_0df_111 * weights_0df[iq];
            fw3 = sv_0df_110 * weights_0df[iq];
            fw4 = sv_0df_109 * weights_0df[iq];
            fw5 = sv_0df_108 * weights_0df[iq];
        }
        // ------------------------
        // ------------------------
        // Section: Tensor Computation
        // Inputs: fw3, FE4_C1_D01_Q0df, fw1, FE4_C1_D10_Q0df, fw0, fw2
        // Outputs: A
        {
            for ( int i = 0; i < 9; ++i ) {
                A[2 * ( i )] += fw0 * FE4_C1_D10_Q0df[0][0][iq][i];
                A[2 * ( i )] += fw1 * FE4_C1_D01_Q0df[0][0][iq][i];
                A[( 2 * ( i ) + 1 )] += fw2 * FE4_C1_D10_Q0df[0][0][iq][i];
                A[( 2 * ( i ) + 1 )] += fw3 * FE4_C1_D01_Q0df[0][0][iq][i];
            }
        }
        // ------------------------
        // ------------------------
        // Section: Tensor Computation
        // Inputs: fw4, fw5, FE4_C4_Q0df, FE4_C3_Q0df
        // Outputs: A
        {
            for ( int i = 0; i < 4; ++i ) {
                A[( ( i ) + 22 )] += fw4 * FE4_C3_Q0df[0][0][iq][i];
                A[( ( i ) + 22 )] += fw5 * FE4_C4_Q0df[0][0][iq][i];
            }
        }
        // ------------------------
    }
    for ( int iq = 0; iq < 9; ++iq ) {
        // ------------------------
        // Section: Jacobian
        // Inputs: FE5_C1_D01_Qe1c, FE5_C0_D10_Qe1c, coordinate_dofs
        // Outputs: J_c3, J_c0, J_c2, J_c1
        double J_c0 = 0.0;
        double J_c3 = 0.0;
        double J_c1 = 0.0;
        double J_c2 = 0.0;
        {
            for ( int ic = 0; ic < 4; ++ic ) {
                J_c0 += coordinate_dofs[(ic)*3] * FE5_C0_D10_Qe1c[0][0][iq][ic];
                J_c3 += coordinate_dofs[(ic)*3 + 1] * FE5_C1_D01_Qe1c[0][0][iq][ic];
                J_c1 += coordinate_dofs[(ic)*3] * FE5_C1_D01_Qe1c[0][0][iq][ic];
                J_c2 += coordinate_dofs[(ic)*3 + 1] * FE5_C0_D10_Qe1c[0][0][iq][ic];
            }
        }
        // ------------------------
        // ------------------------
        // Section: Intermediates
        // Inputs: J_c3, J_c0, J_c2, J_c1
        // Outputs: fw6
        double fw6 = 0;
        {
            double sv_e1c_0 = J_c0 * J_c3;
            double sv_e1c_1 = J_c1 * J_c2;
            double sv_e1c_2 = -sv_e1c_1;
            double sv_e1c_3 = sv_e1c_0 + sv_e1c_2;
            double sv_e1c_4 = fabs( sv_e1c_3 );
            double sv_e1c_5 = sp_e1c_1 * sv_e1c_4;
            fw6 = sv_e1c_5 * weights_e1c[iq];
        }
        // ------------------------
        // ------------------------
        // Section: Tensor Computation
        // Inputs: fw6, FE4_C2_Qe1c
        // Outputs: A
        {
            for ( int i = 0; i < 4; ++i ) {
                A[( ( i ) + 18 )] += fw6 * FE4_C2_Qe1c[0][0][iq][i];
            }
        }
        // ------------------------
    }
    return A;
}

VectorReal B_p2_qu9( const VectorReal &w, const VectorReal &coordinate_dofs,
                     const VectorInt &entity_local_index, const VectorReal &c ) {
    VectorReal A( 42 * 42, 0.0 ); // Vecteur 1D pour stocker les valeurs (size = nbddl)
    // Quadrature rules
    static const double weights_4a8[2] = { 0.5, 0.5 };
    // Precomputed values of basis functions and precomputations
    // FE* dimensions: [permutation][entities][points][dofs]
    static const double FE4_C0_F_Q4a8[1][4][2][9] = {
        { { { 0.4553418012614795, -0.1220084679281462, 0.0, 0.0, 0.6666666666666667, 0.0, 0.0, 0.0,
              0.0 },
            { -0.1220084679281463, 0.4553418012614796, 0.0, 0.0, 0.6666666666666667, 0.0, 0.0, 0.0,
              0.0 } },
          { { 0.4553418012614795, 0.0, -0.1220084679281462, 0.0, 0.0, 0.6666666666666667, 0.0, 0.0,
              0.0 },
            { -0.1220084679281462, 0.0, 0.4553418012614796, 0.0, 0.0, 0.6666666666666667, 0.0, 0.0,
              0.0 } },
          { { 0.0, 0.4553418012614796, 0.0, -0.1220084679281462, 0.0, 0.0, 0.6666666666666667, 0.0,
              0.0 },
            { 0.0, -0.1220084679281461, 0.0, 0.4553418012614795, 0.0, 0.0, 0.6666666666666667, 0.0,
              0.0 } },
          { { 0.0, 0.0, 0.4553418012614796, -0.1220084679281462, 0.0, 0.0, 0.0, 0.6666666666666667,
              0.0 },
            { 0.0, 0.0, -0.1220084679281461, 0.4553418012614794, 0.0, 0.0, 0.0, 0.6666666666666667,
              0.0 } } }
    };
    static const double FE4_C2_D10_F_Q4a8[1][4][2][4] = {
        { { { -1.0, 1.0, 0.0, 0.0 }, { -1.0, 1.0, 0.0, 0.0 } },
          { { -0.7886751345948126, 0.7886751345948129, -0.2113248654051871, 0.2113248654051871 },
            { -0.2113248654051872, 0.2113248654051871, -0.7886751345948129, 0.7886751345948129 } },
          { { -0.7886751345948126, 0.7886751345948129, -0.2113248654051871, 0.2113248654051871 },
            { -0.2113248654051872, 0.2113248654051871, -0.7886751345948129, 0.7886751345948129 } },
          { { 0.0, 0.0, -1.0, 1.0 }, { 0.0, 0.0, -1.0, 1.0 } } }
    };
    static const double FE4_C3_F_Q4a8[1][4][2][4] = {
        { { { 1.0, 0.0, 0.0, 0.0 }, { 1.0, 0.0, 0.0, 0.0 } },
          { { 0.7886751345948129, 0.0, 0.0, 0.2113248654051871 },
            { 0.2113248654051872, 0.0, 0.0, 0.7886751345948129 } },
          { { 0.7886751345948129, 0.0, 0.0, 0.2113248654051871 },
            { 0.2113248654051872, 0.0, 0.0, 0.7886751345948129 } },
          { { 0.0, 0.0, 0.0, 1.0 }, { 0.0, 0.0, 0.0, 1.0 } } }
    };
    static const double FE4_C4_F_Q4a8[1][4][2][4] = {
        { { { 0.0, 0.7886751345948129, 0.2113248654051871, 0.0 },
            { 0.0, 0.2113248654051872, 0.7886751345948129, 0.0 } },
          { { 0.0, 1.0, 0.0, 0.0 }, { 0.0, 1.0, 0.0, 0.0 } },
          { { 0.0, 0.0, 1.0, 0.0 }, { 0.0, 0.0, 1.0, 0.0 } },
          { { 0.0, 0.7886751345948129, 0.2113248654051871, 0.0 },
            { 0.0, 0.2113248654051872, 0.7886751345948129, 0.0 } } }
    };
    static const double FE5_C1_D01_F_Q4a8[1][4][2][4] = {
        { { { -0.7886751345948126, -0.2113248654051871, 0.7886751345948129, 0.2113248654051871 },
            { -0.2113248654051872, -0.7886751345948129, 0.2113248654051871, 0.7886751345948129 } },
          { { -1.0, 0.0, 1.0, 0.0 }, { -1.0, 0.0, 1.0, 0.0 } },
          { { 0.0, -1.0, 0.0, 1.0 }, { 0.0, -1.0, 0.0, 1.0 } },
          { { -0.7886751345948126, -0.2113248654051871, 0.7886751345948129, 0.2113248654051871 },
            { -0.2113248654051872, -0.7886751345948129, 0.2113248654051871, 0.7886751345948129 } } }
    };
    static const double quadrilateral_reference_facet_jacobian[4][2][1] = {
        { { 1.0 }, { 0.0 } }, { { 0.0 }, { 1.0 } }, { { 0.0 }, { 1.0 } }, { { 1.0 }, { 0.0 } }
    };
    static const double quadrilateral_reference_facet_normals[4][2] = {
        { 0.0, -1.0 }, { -1.0, -0.0 }, { 1.0, 0.0 }, { -0.0, 1.0 }
    };
    for ( int iq = 0; iq < 2; ++iq ) {
        // ------------------------
        // Section: Function
        // Inputs: w, FE4_C2_D10_F_Q4a8
        // Outputs: w0_d10_c2
        double w0_d10_c2 = 0.0;
        {
            for ( int ic = 0; ic < 4; ++ic ) {
                w0_d10_c2 += w[( ic ) + 18] * FE4_C2_D10_F_Q4a8[0][entity_local_index[0]][iq][ic];
            }
        }
        // ------------------------
        // ------------------------
        // Section: Jacobian
        // Inputs: FE5_C1_D01_F_Q4a8, FE4_C2_D10_F_Q4a8, coordinate_dofs
        // Outputs: J_c3, J_c0, J_c2, J_c1
        double J_c3 = 0.0;
        double J_c0 = 0.0;
        double J_c1 = 0.0;
        double J_c2 = 0.0;
        {
            for ( int ic = 0; ic < 4; ++ic ) {
                J_c3 += coordinate_dofs[(ic)*3 + 1] *
                        FE5_C1_D01_F_Q4a8[0][entity_local_index[0]][iq][ic];
                J_c0 +=
                    coordinate_dofs[(ic)*3] * FE4_C2_D10_F_Q4a8[0][entity_local_index[0]][iq][ic];
                J_c1 +=
                    coordinate_dofs[(ic)*3] * FE5_C1_D01_F_Q4a8[0][entity_local_index[0]][iq][ic];
                J_c2 += coordinate_dofs[(ic)*3 + 1] *
                        FE4_C2_D10_F_Q4a8[0][entity_local_index[0]][iq][ic];
            }
        }
        // ------------------------
        // ------------------------
        // Section: Function
        // Inputs: w, FE5_C1_D01_F_Q4a8
        // Outputs: w0_d01_c2
        double w0_d01_c2 = 0.0;
        {
            for ( int ic = 0; ic < 4; ++ic ) {
                w0_d01_c2 += w[( ic ) + 18] * FE5_C1_D01_F_Q4a8[0][entity_local_index[0]][iq][ic];
            }
        }
        // ------------------------
        // ------------------------
        // Section: Function
        // Inputs: w, FE4_C0_F_Q4a8
        // Outputs: w0_c0
        double w0_c0 = 0.0;
        {
            for ( int ic = 0; ic < 9; ++ic ) {
                w0_c0 += w[(ic)*2] * FE4_C0_F_Q4a8[0][entity_local_index[0]][iq][ic];
            }
        }
        // ------------------------
        // ------------------------
        // Section: Function
        // Inputs: w, FE4_C3_F_Q4a8
        // Outputs: w0_c3
        double w0_c3 = 0.0;
        {
            for ( int ic = 0; ic < 4; ++ic ) {
                w0_c3 += w[( ic ) + 22] * FE4_C3_F_Q4a8[0][entity_local_index[0]][iq][ic];
            }
        }
        // ------------------------
        // ------------------------
        // Section: Function
        // Inputs: w, FE4_C4_F_Q4a8
        // Outputs: w0_c4
        double w0_c4 = 0.0;
        {
            for ( int ic = 0; ic < 4; ++ic ) {
                w0_c4 += w[( ic ) + 22] * FE4_C4_F_Q4a8[0][entity_local_index[0]][iq][ic];
            }
        }
        // ------------------------
        // ------------------------
        // Section: Function
        // Inputs: w, FE4_C0_F_Q4a8
        // Outputs: w0_c1
        double w0_c1 = 0.0;
        {
            for ( int ic = 0; ic < 9; ++ic ) {
                w0_c1 += w[(ic)*2 + 1] * FE4_C0_F_Q4a8[0][entity_local_index[0]][iq][ic];
            }
        }
        // ------------------------
        // ------------------------
        // Section: Function
        // Inputs: w, FE4_C3_F_Q4a8
        // Outputs: w0_c5
        double w0_c5 = 0.0;
        {
            for ( int ic = 0; ic < 4; ++ic ) {
                w0_c5 += w[( ic ) + 26] * FE4_C3_F_Q4a8[0][entity_local_index[0]][iq][ic];
            }
        }
        // ------------------------
        // ------------------------
        // Section: Function
        // Inputs: w, FE4_C4_F_Q4a8
        // Outputs: w0_c6
        double w0_c6 = 0.0;
        {
            for ( int ic = 0; ic < 4; ++ic ) {
                w0_c6 += w[( ic ) + 26] * FE4_C4_F_Q4a8[0][entity_local_index[0]][iq][ic];
            }
        }
        // ------------------------
        // ------------------------
        // Section: Intermediates
        // Inputs: w0_d10_c2, J_c3, J_c0, J_c2, J_c1, w0_d01_c2, w0_c0, w0_c3, w0_c4, w0_c1, w0_c5,
        // w0_c6 Outputs: fw0, fw1, fw2, fw3, fw4, fw5, fw6, fw7
        double fw0 = 0;
        double fw1 = 0;
        double fw2 = 0;
        double fw3 = 0;
        double fw4 = 0;
        double fw5 = 0;
        double fw6 = 0;
        double fw7 = 0;
        {
            double sv_4a8_0 = J_c0 * J_c3;
            double sv_4a8_1 = J_c1 * J_c2;
            double sv_4a8_2 = -sv_4a8_1;
            double sv_4a8_3 = sv_4a8_0 + sv_4a8_2;
            double sv_4a8_4 = J_c3 / sv_4a8_3;
            double sv_4a8_5 = w0_d10_c2 * sv_4a8_4;
            double sv_4a8_6 = -J_c2;
            double sv_4a8_7 = sv_4a8_6 / sv_4a8_3;
            double sv_4a8_8 = w0_d01_c2 * sv_4a8_7;
            double sv_4a8_9 = sv_4a8_5 + sv_4a8_8;
            double sv_4a8_10 = -w0_c0;
            double sv_4a8_11 = sv_4a8_9 + sv_4a8_10;
            double sv_4a8_12 = w0_c3 * sv_4a8_4;
            double sv_4a8_13 = w0_c4 * sv_4a8_7;
            double sv_4a8_14 = sv_4a8_12 + sv_4a8_13;
            double sv_4a8_15 = -sv_4a8_14;
            double sv_4a8_16 = sv_4a8_11 + sv_4a8_15;
            double sv_4a8_17 = J_c0 / sv_4a8_3;
            double sv_4a8_18 =
                quadrilateral_reference_facet_normals[entity_local_index[0]][1] * sv_4a8_17;
            double sv_4a8_19 = -J_c1;
            double sv_4a8_20 = sv_4a8_19 / sv_4a8_3;
            double sv_4a8_21 =
                quadrilateral_reference_facet_normals[entity_local_index[0]][0] * sv_4a8_20;
            double sv_4a8_22 = sv_4a8_18 + sv_4a8_21;
            double sv_4a8_23 = sv_4a8_22 * sv_4a8_22;
            double sv_4a8_24 =
                quadrilateral_reference_facet_normals[entity_local_index[0]][0] * sv_4a8_4;
            double sv_4a8_25 =
                quadrilateral_reference_facet_normals[entity_local_index[0]][1] * sv_4a8_7;
            double sv_4a8_26 = sv_4a8_24 + sv_4a8_25;
            double sv_4a8_27 = sv_4a8_26 * sv_4a8_26;
            double sv_4a8_28 = sv_4a8_23 + sv_4a8_27;
            double sv_4a8_29 = sqrt( sv_4a8_28 );
            double sv_4a8_30 = sv_4a8_22 / sv_4a8_29;
            double sv_4a8_31 = -sv_4a8_30;
            double sv_4a8_32 = sv_4a8_16 * sv_4a8_31;
            double sv_4a8_33 = w0_d01_c2 * sv_4a8_17;
            double sv_4a8_34 = w0_d10_c2 * sv_4a8_20;
            double sv_4a8_35 = sv_4a8_33 + sv_4a8_34;
            double sv_4a8_36 = -w0_c1;
            double sv_4a8_37 = sv_4a8_35 + sv_4a8_36;
            double sv_4a8_38 = w0_c4 * sv_4a8_17;
            double sv_4a8_39 = w0_c3 * sv_4a8_20;
            double sv_4a8_40 = sv_4a8_38 + sv_4a8_39;
            double sv_4a8_41 = -sv_4a8_40;
            double sv_4a8_42 = sv_4a8_37 + sv_4a8_41;
            double sv_4a8_43 = sv_4a8_26 / sv_4a8_29;
            double sv_4a8_44 = sv_4a8_42 * sv_4a8_43;
            double sv_4a8_45 = sv_4a8_32 + sv_4a8_44;
            double sv_4a8_46 = sv_4a8_31 * sv_4a8_4;
            double sv_4a8_47 = sv_4a8_31 * sv_4a8_7;
            double sv_4a8_48 = sv_4a8_20 * sv_4a8_43;
            double sv_4a8_49 = sv_4a8_17 * sv_4a8_43;
            double sv_4a8_50 = sv_4a8_46 + sv_4a8_48;
            double sv_4a8_51 = sv_4a8_47 + sv_4a8_49;
            double sv_4a8_52 = sv_4a8_45 * sv_4a8_50;
            double sv_4a8_53 = sv_4a8_45 * sv_4a8_51;
            double sv_4a8_54 = -sv_4a8_4;
            double sv_4a8_55 = -sv_4a8_7;
            double sv_4a8_56 = -sv_4a8_31;
            double sv_4a8_57 = sv_4a8_54 * sv_4a8_31;
            double sv_4a8_58 = sv_4a8_55 * sv_4a8_31;
            double sv_4a8_59 = -sv_4a8_20;
            double sv_4a8_60 = -sv_4a8_17;
            double sv_4a8_61 = sv_4a8_59 * sv_4a8_43;
            double sv_4a8_62 = sv_4a8_60 * sv_4a8_43;
            double sv_4a8_63 = -sv_4a8_43;
            double sv_4a8_64 = sv_4a8_57 + sv_4a8_61;
            double sv_4a8_65 = sv_4a8_58 + sv_4a8_62;
            double sv_4a8_66 = w0_c5 * sv_4a8_4;
            double sv_4a8_67 = w0_c6 * sv_4a8_7;
            double sv_4a8_68 = sv_4a8_66 + sv_4a8_67;
            double sv_4a8_69 = sv_4a8_68 * sv_4a8_31;
            double sv_4a8_70 = w0_c6 * sv_4a8_17;
            double sv_4a8_71 = w0_c5 * sv_4a8_20;
            double sv_4a8_72 = sv_4a8_70 + sv_4a8_71;
            double sv_4a8_73 = sv_4a8_72 * sv_4a8_43;
            double sv_4a8_74 = sv_4a8_69 + sv_4a8_73;
            double sv_4a8_75 = sv_4a8_74 * sv_4a8_50;
            double sv_4a8_76 = sv_4a8_74 * sv_4a8_51;
            double sv_4a8_77 = sv_4a8_74 * sv_4a8_56;
            double sv_4a8_78 = sv_4a8_74 * sv_4a8_64;
            double sv_4a8_79 = sv_4a8_74 * sv_4a8_65;
            double sv_4a8_80 = sv_4a8_74 * sv_4a8_63;
            double sv_4a8_81 =
                J_c0 * quadrilateral_reference_facet_jacobian[entity_local_index[0]][0][0];
            double sv_4a8_82 =
                J_c1 * quadrilateral_reference_facet_jacobian[entity_local_index[0]][1][0];
            double sv_4a8_83 = sv_4a8_81 + sv_4a8_82;
            double sv_4a8_84 = sv_4a8_83 * sv_4a8_83;
            double sv_4a8_85 =
                quadrilateral_reference_facet_jacobian[entity_local_index[0]][0][0] * J_c2;
            double sv_4a8_86 =
                quadrilateral_reference_facet_jacobian[entity_local_index[0]][1][0] * J_c3;
            double sv_4a8_87 = sv_4a8_85 + sv_4a8_86;
            double sv_4a8_88 = sv_4a8_87 * sv_4a8_87;
            double sv_4a8_89 = sv_4a8_84 + sv_4a8_88;
            double sv_4a8_90 = sqrt( sv_4a8_89 );
            double sv_4a8_91 = sv_4a8_52 * sv_4a8_90;
            double sv_4a8_92 = sv_4a8_53 * sv_4a8_90;
            double sv_4a8_93 = sv_4a8_75 * sv_4a8_90;
            double sv_4a8_94 = sv_4a8_76 * sv_4a8_90;
            double sv_4a8_95 = sv_4a8_77 * sv_4a8_90;
            double sv_4a8_96 = sv_4a8_78 * sv_4a8_90;
            double sv_4a8_97 = sv_4a8_79 * sv_4a8_90;
            double sv_4a8_98 = sv_4a8_80 * sv_4a8_90;
            fw0 = sv_4a8_95 * weights_4a8[iq];
            fw1 = sv_4a8_98 * weights_4a8[iq];
            fw2 = sv_4a8_93 * weights_4a8[iq];
            fw3 = sv_4a8_94 * weights_4a8[iq];
            fw4 = sv_4a8_96 * weights_4a8[iq];
            fw5 = sv_4a8_97 * weights_4a8[iq];
            fw6 = sv_4a8_91 * weights_4a8[iq];
            fw7 = sv_4a8_92 * weights_4a8[iq];
        }
        // ------------------------
        // ------------------------
        // Section: Tensor Computation
        // Inputs: fw1, fw0, FE4_C0_F_Q4a8
        // Outputs: A
        {
            for ( int i = 0; i < 9; ++i ) {
                A[2 * ( i )] += fw0 * FE4_C0_F_Q4a8[0][entity_local_index[0]][iq][i];
                A[( 2 * ( i ) + 1 )] += fw1 * FE4_C0_F_Q4a8[0][entity_local_index[0]][iq][i];
            }
        }
        // ------------------------
        // ------------------------
        // Section: Tensor Computation
        // Inputs: FE5_C1_D01_F_Q4a8, fw6, FE4_C2_D10_F_Q4a8, fw7, FE4_C3_F_Q4a8, fw5,
        // FE4_C4_F_Q4a8, fw3, fw4, fw2 Outputs: A
        {
            for ( int i = 0; i < 4; ++i ) {
                A[( ( i ) + 18 )] += fw2 * FE4_C2_D10_F_Q4a8[0][entity_local_index[0]][iq][i];
                A[( ( i ) + 18 )] += fw3 * FE5_C1_D01_F_Q4a8[0][entity_local_index[0]][iq][i];
                A[( ( i ) + 22 )] += fw4 * FE4_C3_F_Q4a8[0][entity_local_index[0]][iq][i];
                A[( ( i ) + 22 )] += fw5 * FE4_C4_F_Q4a8[0][entity_local_index[0]][iq][i];
                A[( ( i ) + 26 )] += fw6 * FE4_C3_F_Q4a8[0][entity_local_index[0]][iq][i];
                A[( ( i ) + 26 )] += fw7 * FE4_C4_F_Q4a8[0][entity_local_index[0]][iq][i];
            }
        }
        // ------------------------
    }
    return A;
}

VectorReal B_p4_qu9( const VectorReal &w, const VectorReal &coordinate_dofs,
                     const VectorInt &entity_local_index, const VectorReal &c ) {
    VectorReal A( 42 * 42, 0.0 ); // Vecteur 1D pour stocker les valeurs (size = nbddl)
    static const double weights_0df[4] = { 0.25, 0.25, 0.25, 0.25 };
    // Precomputed values of basis functions and precomputations
    // FE* dimensions: [permutation][entities][points][dofs]
    static const double FE1_C0_D10_Q0df[1][1][4][4] = {
        { { { -0.7886751345948126, 0.7886751345948129, -0.2113248654051871, 0.2113248654051871 },
            { -0.2113248654051872, 0.2113248654051871, -0.7886751345948129, 0.7886751345948129 },
            { -0.7886751345948126, 0.7886751345948129, -0.2113248654051871, 0.2113248654051871 },
            { -0.2113248654051872, 0.2113248654051871, -0.7886751345948129, 0.7886751345948129 } } }
    };
    static const double FE1_C1_D01_Q0df[1][1][4][4] = {
        { { { -0.7886751345948126, -0.2113248654051871, 0.7886751345948129, 0.2113248654051871 },
            { -0.7886751345948126, -0.2113248654051871, 0.7886751345948129, 0.2113248654051871 },
            { -0.2113248654051872, -0.7886751345948129, 0.2113248654051871, 0.7886751345948129 },
            { -0.2113248654051872, -0.7886751345948129, 0.2113248654051871, 0.7886751345948129 } } }
    };
    static const double FE5_C1_D01_Q0df[1][1][4][9] = {
        { { { -0.9811252243246881, 0.2628917115316043, -0.07044162180172914, 0.01887477567531192,
              -1.436467025586168, 1.051566846126417, -0.2817664872069161, -0.1031336922528346,
              1.539600717839002 },
            { 0.07044162180172903, -0.01887477567531187, 0.9811252243246882, -0.2628917115316043,
              0.1031336922528346, -1.051566846126417, 0.2817664872069161, 1.436467025586168,
              -1.539600717839002 },
            { 0.2628917115316045, -0.9811252243246883, 0.01887477567531164, -0.07044162180172886,
              -1.436467025586168, -0.2817664872069161, 1.051566846126417, -0.1031336922528346,
              1.539600717839002 },
            { -0.01887477567531198, 0.07044162180172914, -0.2628917115316041, 0.981125224324688,
              0.1031336922528346, 0.2817664872069162, -1.051566846126417, 1.436467025586168,
              -1.539600717839002 } } }
    };
    static const double FE5_C1_D10_Q0df[1][1][4][9] = {
        { { { -0.9811252243246882, -0.07044162180172919, 0.2628917115316043, 0.01887477567531193,
              1.051566846126417, -1.436467025586168, -0.1031336922528345, -0.2817664872069161,
              1.539600717839002 },
            { 0.2628917115316044, 0.01887477567531159, -0.9811252243246882, -0.07044162180172903,
              -0.2817664872069162, -1.436467025586168, -0.1031336922528345, 1.051566846126417,
              1.539600717839002 },
            { 0.07044162180172914, 0.9811252243246882, -0.01887477567531187, -0.2628917115316043,
              -1.051566846126417, 0.1031336922528346, 1.436467025586168, 0.2817664872069161,
              -1.539600717839002 },
            { -0.01887477567531198, -0.2628917115316041, 0.07044162180172925, 0.981125224324688,
              0.2817664872069162, 0.1031336922528345, 1.436467025586168, -1.051566846126417,
              -1.539600717839002 } } }
    };
    static const double FE5_C3_Q0df[1][1][4][4] = {
        { { { 0.7886751345948129, 0.0, 0.0, 0.2113248654051871 },
            { 0.2113248654051872, 0.0, 0.0, 0.7886751345948129 },
            { 0.7886751345948129, 0.0, 0.0, 0.2113248654051871 },
            { 0.2113248654051872, 0.0, 0.0, 0.7886751345948129 } } }
    };
    static const double FE5_C4_Q0df[1][1][4][4] = {
        { { { 0.0, 0.7886751345948129, 0.2113248654051871, 0.0 },
            { 0.0, 0.7886751345948129, 0.2113248654051871, 0.0 },
            { 0.0, 0.2113248654051872, 0.7886751345948129, 0.0 },
            { 0.0, 0.2113248654051872, 0.7886751345948129, 0.0 } } }
    };
    double sp_0df_0 = c[0] * c[2];
    double sp_0df_1 = c[3] * sp_0df_0;
    double sp_0df_2 = 1.0 + c[1];
    double sp_0df_3 = 4.0 * sp_0df_2;
    double sp_0df_4 = sp_0df_1 / sp_0df_3;
    double sp_0df_5 = -c[1];
    double sp_0df_6 = 1.0 + sp_0df_5;
    double sp_0df_7 = pow( c[3], 3 );
    double sp_0df_8 = c[0] * sp_0df_7;
    double sp_0df_9 = pow( c[1], 2 );
    double sp_0df_10 = -sp_0df_9;
    double sp_0df_11 = 1.0 + sp_0df_10;
    double sp_0df_12 = 24.0 * sp_0df_11;
    double sp_0df_13 = sp_0df_8 / sp_0df_12;
    for ( int iq = 0; iq < 4; ++iq ) {
        // ------------------------
        // Section: Jacobian
        // Inputs: FE1_C1_D01_Q0df, FE1_C0_D10_Q0df, coordinate_dofs
        // Outputs: J_c2, J_c1, J_c0, J_c3
        double J_c0 = 0.0;
        double J_c3 = 0.0;
        double J_c1 = 0.0;
        double J_c2 = 0.0;
        {
            for ( int ic = 0; ic < 4; ++ic ) {
                J_c0 += coordinate_dofs[(ic)*3] * FE1_C0_D10_Q0df[0][0][iq][ic];
                J_c3 += coordinate_dofs[(ic)*3 + 1] * FE1_C1_D01_Q0df[0][0][iq][ic];
                J_c1 += coordinate_dofs[(ic)*3] * FE1_C1_D01_Q0df[0][0][iq][ic];
                J_c2 += coordinate_dofs[(ic)*3 + 1] * FE1_C0_D10_Q0df[0][0][iq][ic];
            }
        }
        // ------------------------
        // ------------------------
        // Section: Intermediates
        // Inputs: J_c2, J_c1, J_c0, J_c3
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
            double sv_0df_0 = J_c0 * J_c3;
            double sv_0df_1 = J_c1 * J_c2;
            double sv_0df_2 = -sv_0df_1;
            double sv_0df_3 = sv_0df_0 + sv_0df_2;
            double sv_0df_4 = J_c0 / sv_0df_3;
            double sv_0df_5 = -J_c1;
            double sv_0df_6 = sv_0df_5 / sv_0df_3;
            double sv_0df_7 = sv_0df_4 * sv_0df_4;
            double sv_0df_8 = sv_0df_4 * sv_0df_6;
            double sv_0df_9 = sv_0df_6 * sv_0df_6;
            double sv_0df_10 = sv_0df_7 + sv_0df_7;
            double sv_0df_11 = sv_0df_8 + sv_0df_8;
            double sv_0df_12 = sv_0df_9 + sv_0df_9;
            double sv_0df_13 = J_c3 / sv_0df_3;
            double sv_0df_14 = -J_c2;
            double sv_0df_15 = sv_0df_14 / sv_0df_3;
            double sv_0df_16 = sv_0df_15 * sv_0df_15;
            double sv_0df_17 = sv_0df_13 * sv_0df_15;
            double sv_0df_18 = sv_0df_13 * sv_0df_13;
            double sv_0df_19 = sv_0df_16 + sv_0df_16;
            double sv_0df_20 = sv_0df_17 + sv_0df_17;
            double sv_0df_21 = sv_0df_18 + sv_0df_18;
            double sv_0df_22 = sv_0df_10 + sv_0df_19;
            double sv_0df_23 = sv_0df_11 + sv_0df_20;
            double sv_0df_24 = sv_0df_21 + sv_0df_12;
            double sv_0df_25 = sv_0df_22 * sp_0df_4;
            double sv_0df_26 = sv_0df_23 * sp_0df_4;
            double sv_0df_27 = sv_0df_24 * sp_0df_4;
            double sv_0df_28 = sv_0df_4 + sv_0df_4;
            double sv_0df_29 = sv_0df_6 + sv_0df_6;
            double sv_0df_30 = sv_0df_28 / 2;
            double sv_0df_31 = sv_0df_29 / 2;
            double sv_0df_32 = sv_0df_30 * sv_0df_30;
            double sv_0df_33 = sv_0df_30 * sv_0df_31;
            double sv_0df_34 = sv_0df_31 * sv_0df_31;
            double sv_0df_35 = sv_0df_32 + sv_0df_32;
            double sv_0df_36 = sv_0df_33 + sv_0df_33;
            double sv_0df_37 = sv_0df_34 + sv_0df_34;
            double sv_0df_38 = sv_0df_15 / 2;
            double sv_0df_39 = sv_0df_13 / 2;
            double sv_0df_40 = sv_0df_4 / 2;
            double sv_0df_41 = sv_0df_6 / 2;
            double sv_0df_42 = sv_0df_38 * sv_0df_38;
            double sv_0df_43 = sv_0df_39 * sv_0df_38;
            double sv_0df_44 = sv_0df_40 * sv_0df_38;
            double sv_0df_45 = sv_0df_41 * sv_0df_38;
            double sv_0df_46 = sv_0df_39 * sv_0df_39;
            double sv_0df_47 = sv_0df_40 * sv_0df_39;
            double sv_0df_48 = sv_0df_39 * sv_0df_41;
            double sv_0df_49 = sv_0df_40 * sv_0df_40;
            double sv_0df_50 = sv_0df_40 * sv_0df_41;
            double sv_0df_51 = sv_0df_41 * sv_0df_41;
            double sv_0df_52 = sv_0df_42 + sv_0df_42;
            double sv_0df_53 = sv_0df_43 + sv_0df_43;
            double sv_0df_54 = sv_0df_44 + sv_0df_44;
            double sv_0df_55 = sv_0df_45 + sv_0df_45;
            double sv_0df_56 = sv_0df_46 + sv_0df_46;
            double sv_0df_57 = sv_0df_47 + sv_0df_47;
            double sv_0df_58 = sv_0df_48 + sv_0df_48;
            double sv_0df_59 = sv_0df_49 + sv_0df_49;
            double sv_0df_60 = sv_0df_50 + sv_0df_50;
            double sv_0df_61 = sv_0df_51 + sv_0df_51;
            double sv_0df_62 = sv_0df_35 + sv_0df_52;
            double sv_0df_63 = sv_0df_36 + sv_0df_53;
            double sv_0df_64 = sv_0df_37 + sv_0df_56;
            double sv_0df_65 = sv_0df_15 + sv_0df_15;
            double sv_0df_66 = sv_0df_13 + sv_0df_13;
            double sv_0df_67 = sv_0df_65 / 2;
            double sv_0df_68 = sv_0df_66 / 2;
            double sv_0df_69 = sv_0df_67 * sv_0df_67;
            double sv_0df_70 = sv_0df_68 * sv_0df_67;
            double sv_0df_71 = sv_0df_68 * sv_0df_68;
            double sv_0df_72 = sv_0df_69 + sv_0df_69;
            double sv_0df_73 = sv_0df_70 + sv_0df_70;
            double sv_0df_74 = sv_0df_71 + sv_0df_71;
            double sv_0df_75 = sv_0df_72 + sv_0df_59;
            double sv_0df_76 = sv_0df_73 + sv_0df_60;
            double sv_0df_77 = sv_0df_74 + sv_0df_61;
            double sv_0df_78 = sv_0df_62 + sv_0df_52;
            double sv_0df_79 = sv_0df_63 + sv_0df_53;
            double sv_0df_80 = sv_0df_54 + sv_0df_54;
            double sv_0df_81 = sv_0df_55 + sv_0df_55;
            double sv_0df_82 = sv_0df_64 + sv_0df_56;
            double sv_0df_83 = sv_0df_57 + sv_0df_57;
            double sv_0df_84 = sv_0df_58 + sv_0df_58;
            double sv_0df_85 = sv_0df_75 + sv_0df_59;
            double sv_0df_86 = sv_0df_76 + sv_0df_60;
            double sv_0df_87 = sv_0df_77 + sv_0df_61;
            double sv_0df_88 = sv_0df_78 * sp_0df_6;
            double sv_0df_89 = sv_0df_79 * sp_0df_6;
            double sv_0df_90 = sv_0df_80 * sp_0df_6;
            double sv_0df_91 = sv_0df_81 * sp_0df_6;
            double sv_0df_92 = sv_0df_82 * sp_0df_6;
            double sv_0df_93 = sv_0df_83 * sp_0df_6;
            double sv_0df_94 = sv_0df_84 * sp_0df_6;
            double sv_0df_95 = sv_0df_85 * sp_0df_6;
            double sv_0df_96 = sv_0df_86 * sp_0df_6;
            double sv_0df_97 = sv_0df_87 * sp_0df_6;
            double sv_0df_98 = 2 * sv_0df_30;
            double sv_0df_99 = 2 * sv_0df_31;
            double sv_0df_100 = 2 * sv_0df_67;
            double sv_0df_101 = 2 * sv_0df_68;
            double sv_0df_102 = sv_0df_98 * sv_0df_30;
            double sv_0df_103 = sv_0df_99 * sv_0df_30;
            double sv_0df_104 = sv_0df_100 * sv_0df_30;
            double sv_0df_105 = sv_0df_101 * sv_0df_30;
            double sv_0df_106 = sv_0df_98 * sv_0df_31;
            double sv_0df_107 = sv_0df_99 * sv_0df_31;
            double sv_0df_108 = sv_0df_100 * sv_0df_31;
            double sv_0df_109 = sv_0df_101 * sv_0df_31;
            double sv_0df_110 = sv_0df_98 * sv_0df_67;
            double sv_0df_111 = sv_0df_99 * sv_0df_67;
            double sv_0df_112 = sv_0df_100 * sv_0df_67;
            double sv_0df_113 = sv_0df_101 * sv_0df_67;
            double sv_0df_114 = sv_0df_98 * sv_0df_68;
            double sv_0df_115 = sv_0df_99 * sv_0df_68;
            double sv_0df_116 = sv_0df_100 * sv_0df_68;
            double sv_0df_117 = sv_0df_101 * sv_0df_68;
            double sv_0df_118 = c[1] * sv_0df_102;
            double sv_0df_119 = c[1] * sv_0df_106;
            double sv_0df_120 = c[1] * sv_0df_110;
            double sv_0df_121 = c[1] * sv_0df_114;
            double sv_0df_122 = c[1] * sv_0df_103;
            double sv_0df_123 = c[1] * sv_0df_107;
            double sv_0df_124 = c[1] * sv_0df_111;
            double sv_0df_125 = c[1] * sv_0df_115;
            double sv_0df_126 = c[1] * sv_0df_104;
            double sv_0df_127 = c[1] * sv_0df_105;
            double sv_0df_128 = c[1] * sv_0df_108;
            double sv_0df_129 = c[1] * sv_0df_109;
            double sv_0df_130 = c[1] * sv_0df_112;
            double sv_0df_131 = c[1] * sv_0df_116;
            double sv_0df_132 = c[1] * sv_0df_113;
            double sv_0df_133 = c[1] * sv_0df_117;
            double sv_0df_134 = sv_0df_88 + sv_0df_118;
            double sv_0df_135 = sv_0df_89 + sv_0df_119;
            double sv_0df_136 = sv_0df_90 + sv_0df_120;
            double sv_0df_137 = sv_0df_91 + sv_0df_121;
            double sv_0df_138 = sv_0df_89 + sv_0df_122;
            double sv_0df_139 = sv_0df_92 + sv_0df_123;
            double sv_0df_140 = sv_0df_93 + sv_0df_124;
            double sv_0df_141 = sv_0df_94 + sv_0df_125;
            double sv_0df_142 = sv_0df_90 + sv_0df_126;
            double sv_0df_143 = sv_0df_91 + sv_0df_127;
            double sv_0df_144 = sv_0df_93 + sv_0df_128;
            double sv_0df_145 = sv_0df_94 + sv_0df_129;
            double sv_0df_146 = sv_0df_95 + sv_0df_130;
            double sv_0df_147 = sv_0df_96 + sv_0df_131;
            double sv_0df_148 = sv_0df_96 + sv_0df_132;
            double sv_0df_149 = sv_0df_97 + sv_0df_133;
            double sv_0df_150 = sv_0df_134 * sp_0df_13;
            double sv_0df_151 = sv_0df_135 * sp_0df_13;
            double sv_0df_152 = sv_0df_136 * sp_0df_13;
            double sv_0df_153 = sv_0df_137 * sp_0df_13;
            double sv_0df_154 = sv_0df_138 * sp_0df_13;
            double sv_0df_155 = sv_0df_139 * sp_0df_13;
            double sv_0df_156 = sv_0df_140 * sp_0df_13;
            double sv_0df_157 = sv_0df_141 * sp_0df_13;
            double sv_0df_158 = sv_0df_142 * sp_0df_13;
            double sv_0df_159 = sv_0df_143 * sp_0df_13;
            double sv_0df_160 = sv_0df_144 * sp_0df_13;
            double sv_0df_161 = sv_0df_145 * sp_0df_13;
            double sv_0df_162 = sv_0df_146 * sp_0df_13;
            double sv_0df_163 = sv_0df_147 * sp_0df_13;
            double sv_0df_164 = sv_0df_148 * sp_0df_13;
            double sv_0df_165 = sv_0df_149 * sp_0df_13;
            double sv_0df_166 = fabs( sv_0df_3 );
            double sv_0df_167 = sv_0df_25 * sv_0df_166;
            double sv_0df_168 = sv_0df_26 * sv_0df_166;
            double sv_0df_169 = sv_0df_27 * sv_0df_166;
            double sv_0df_170 = sv_0df_150 * sv_0df_166;
            double sv_0df_171 = sv_0df_151 * sv_0df_166;
            double sv_0df_172 = sv_0df_152 * sv_0df_166;
            double sv_0df_173 = sv_0df_153 * sv_0df_166;
            double sv_0df_174 = sv_0df_154 * sv_0df_166;
            double sv_0df_175 = sv_0df_155 * sv_0df_166;
            double sv_0df_176 = sv_0df_156 * sv_0df_166;
            double sv_0df_177 = sv_0df_157 * sv_0df_166;
            double sv_0df_178 = sv_0df_158 * sv_0df_166;
            double sv_0df_179 = sv_0df_159 * sv_0df_166;
            double sv_0df_180 = sv_0df_160 * sv_0df_166;
            double sv_0df_181 = sv_0df_161 * sv_0df_166;
            double sv_0df_182 = sv_0df_162 * sv_0df_166;
            double sv_0df_183 = sv_0df_163 * sv_0df_166;
            double sv_0df_184 = sv_0df_164 * sv_0df_166;
            double sv_0df_185 = sv_0df_165 * sv_0df_166;
            fw0 = sv_0df_185 * weights_0df[iq];
            fw1 = sv_0df_184 * weights_0df[iq];
            fw2 = sv_0df_183 * weights_0df[iq];
            fw3 = sv_0df_182 * weights_0df[iq];
            fw4 = sv_0df_181 * weights_0df[iq];
            fw5 = sv_0df_179 * weights_0df[iq];
            fw6 = sv_0df_180 * weights_0df[iq];
            fw7 = sv_0df_178 * weights_0df[iq];
            fw8 = sv_0df_177 * weights_0df[iq];
            fw9 = sv_0df_176 * weights_0df[iq];
            fw10 = sv_0df_173 * weights_0df[iq];
            fw11 = sv_0df_172 * weights_0df[iq];
            fw12 = sv_0df_175 * weights_0df[iq];
            fw13 = sv_0df_174 * weights_0df[iq];
            fw14 = sv_0df_171 * weights_0df[iq];
            fw15 = sv_0df_170 * weights_0df[iq];
            fw16 = sv_0df_169 * weights_0df[iq];
            fw17 = sv_0df_168 * weights_0df[iq];
            fw18 = sv_0df_167 * weights_0df[iq];
        }
        // ------------------------
        // ------------------------
        // Section: Tensor Computation
        // Inputs: fw11, fw2, fw5, fw3, fw1, fw14, FE5_C1_D10_Q0df, FE5_C1_D01_Q0df, fw7, fw12, fw4,
        // fw8, fw10, fw9, fw15, fw0, fw13, fw6 Outputs: A
        {
            double temp_0[9] = { 0 };
            for ( int j = 0; j < 9; ++j ) {
                temp_0[j] = fw0 * FE5_C1_D10_Q0df[0][0][iq][j];
            }
            double temp_1[9] = { 0 };
            for ( int j = 0; j < 9; ++j ) {
                temp_1[j] = fw1 * FE5_C1_D01_Q0df[0][0][iq][j];
            }
            double temp_2[9] = { 0 };
            for ( int j = 0; j < 9; ++j ) {
                temp_2[j] = fw2 * FE5_C1_D10_Q0df[0][0][iq][j];
            }
            double temp_3[9] = { 0 };
            for ( int j = 0; j < 9; ++j ) {
                temp_3[j] = fw3 * FE5_C1_D01_Q0df[0][0][iq][j];
            }
            double temp_4[9] = { 0 };
            for ( int j = 0; j < 9; ++j ) {
                temp_4[j] = fw4 * FE5_C1_D10_Q0df[0][0][iq][j];
            }
            double temp_5[9] = { 0 };
            for ( int j = 0; j < 9; ++j ) {
                temp_5[j] = fw5 * FE5_C1_D01_Q0df[0][0][iq][j];
            }
            double temp_6[9] = { 0 };
            for ( int j = 0; j < 9; ++j ) {
                temp_6[j] = fw6 * FE5_C1_D10_Q0df[0][0][iq][j];
            }
            double temp_7[9] = { 0 };
            for ( int j = 0; j < 9; ++j ) {
                temp_7[j] = fw7 * FE5_C1_D01_Q0df[0][0][iq][j];
            }
            double temp_8[9] = { 0 };
            for ( int j = 0; j < 9; ++j ) {
                temp_8[j] = fw8 * FE5_C1_D10_Q0df[0][0][iq][j];
            }
            double temp_9[9] = { 0 };
            for ( int j = 0; j < 9; ++j ) {
                temp_9[j] = fw9 * FE5_C1_D01_Q0df[0][0][iq][j];
            }
            double temp_10[9] = { 0 };
            for ( int j = 0; j < 9; ++j ) {
                temp_10[j] = fw10 * FE5_C1_D10_Q0df[0][0][iq][j];
            }
            double temp_11[9] = { 0 };
            for ( int j = 0; j < 9; ++j ) {
                temp_11[j] = fw11 * FE5_C1_D01_Q0df[0][0][iq][j];
            }
            double temp_12[9] = { 0 };
            for ( int j = 0; j < 9; ++j ) {
                temp_12[j] = fw12 * FE5_C1_D10_Q0df[0][0][iq][j];
            }
            double temp_13[9] = { 0 };
            for ( int j = 0; j < 9; ++j ) {
                temp_13[j] = fw13 * FE5_C1_D01_Q0df[0][0][iq][j];
            }
            double temp_14[9] = { 0 };
            for ( int j = 0; j < 9; ++j ) {
                temp_14[j] = fw14 * FE5_C1_D10_Q0df[0][0][iq][j];
            }
            double temp_15[9] = { 0 };
            for ( int j = 0; j < 9; ++j ) {
                temp_15[j] = fw15 * FE5_C1_D01_Q0df[0][0][iq][j];
            }
            for ( int j = 0; j < 9; ++j ) {
                for ( int i = 0; i < 9; ++i ) {
                    A[30 * ( 2 * ( i ) ) + 2 * ( j )] += FE5_C1_D10_Q0df[0][0][iq][i] * temp_0[j];
                    A[30 * ( 2 * ( i ) ) + 2 * ( j )] += FE5_C1_D10_Q0df[0][0][iq][i] * temp_1[j];
                    A[30 * ( 2 * ( i ) ) + 2 * ( j )] += FE5_C1_D01_Q0df[0][0][iq][i] * temp_2[j];
                    A[30 * ( 2 * ( i ) ) + 2 * ( j )] += FE5_C1_D01_Q0df[0][0][iq][i] * temp_3[j];
                    A[30 * ( 2 * ( i ) ) + ( 2 * ( j ) + 1 )] +=
                        FE5_C1_D10_Q0df[0][0][iq][i] * temp_4[j];
                    A[30 * ( 2 * ( i ) ) + ( 2 * ( j ) + 1 )] +=
                        FE5_C1_D10_Q0df[0][0][iq][i] * temp_5[j];
                    A[30 * ( 2 * ( i ) ) + ( 2 * ( j ) + 1 )] +=
                        FE5_C1_D01_Q0df[0][0][iq][i] * temp_6[j];
                    A[30 * ( 2 * ( i ) ) + ( 2 * ( j ) + 1 )] +=
                        FE5_C1_D01_Q0df[0][0][iq][i] * temp_7[j];
                    A[30 * ( 2 * ( i ) + 1 ) + 2 * ( j )] +=
                        FE5_C1_D10_Q0df[0][0][iq][i] * temp_8[j];
                    A[30 * ( 2 * ( i ) + 1 ) + 2 * ( j )] +=
                        FE5_C1_D10_Q0df[0][0][iq][i] * temp_9[j];
                    A[30 * ( 2 * ( i ) + 1 ) + 2 * ( j )] +=
                        FE5_C1_D01_Q0df[0][0][iq][i] * temp_10[j];
                    A[30 * ( 2 * ( i ) + 1 ) + 2 * ( j )] +=
                        FE5_C1_D01_Q0df[0][0][iq][i] * temp_11[j];
                    A[30 * ( 2 * ( i ) + 1 ) + ( 2 * ( j ) + 1 )] +=
                        FE5_C1_D10_Q0df[0][0][iq][i] * temp_12[j];
                    A[30 * ( 2 * ( i ) + 1 ) + ( 2 * ( j ) + 1 )] +=
                        FE5_C1_D10_Q0df[0][0][iq][i] * temp_13[j];
                    A[30 * ( 2 * ( i ) + 1 ) + ( 2 * ( j ) + 1 )] +=
                        FE5_C1_D01_Q0df[0][0][iq][i] * temp_14[j];
                    A[30 * ( 2 * ( i ) + 1 ) + ( 2 * ( j ) + 1 )] +=
                        FE5_C1_D01_Q0df[0][0][iq][i] * temp_15[j];
                }
            }
        }
        // ------------------------
        // ------------------------
        // Section: Tensor Computation
        // Inputs: FE5_C3_Q0df, fw18, fw16, FE5_C4_Q0df, fw17
        // Outputs: A
        {
            double temp_0[4] = { 0 };
            for ( int j = 0; j < 4; ++j ) {
                temp_0[j] = fw16 * FE5_C3_Q0df[0][0][iq][j];
            }
            double temp_1[4] = { 0 };
            for ( int j = 0; j < 4; ++j ) {
                temp_1[j] = fw17 * FE5_C4_Q0df[0][0][iq][j];
            }
            double temp_2[4] = { 0 };
            for ( int j = 0; j < 4; ++j ) {
                temp_2[j] = fw17 * FE5_C3_Q0df[0][0][iq][j];
            }
            double temp_3[4] = { 0 };
            for ( int j = 0; j < 4; ++j ) {
                temp_3[j] = fw18 * FE5_C4_Q0df[0][0][iq][j];
            }
            for ( int j = 0; j < 4; ++j ) {
                for ( int i = 0; i < 4; ++i ) {
                    A[30 * ( ( i ) + 22 ) + ( ( j ) + 22 )] += FE5_C3_Q0df[0][0][iq][i] * temp_0[j];
                    A[30 * ( ( i ) + 22 ) + ( ( j ) + 22 )] += FE5_C3_Q0df[0][0][iq][i] * temp_1[j];
                    A[30 * ( ( i ) + 22 ) + ( ( j ) + 22 )] += FE5_C4_Q0df[0][0][iq][i] * temp_2[j];
                    A[30 * ( ( i ) + 22 ) + ( ( j ) + 22 )] += FE5_C4_Q0df[0][0][iq][i] * temp_3[j];
                }
            }
        }
        // ------------------------
    }
    return A;
}

VectorReal B_p5_qu9( const VectorReal &w, const VectorReal &coordinate_dofs,
                     const VectorInt &entity_local_index, const VectorReal &c ) {
    VectorReal A( 42 * 42, 0.0 ); // Vecteur 1D pour stocker les valeurs (size = nbddl)

    // Quadrature rules
    static const double weights_4a8[2] = { 0.5, 0.5 };
    // Precomputed values of basis functions and precomputations
    // FE* dimensions: [permutation][entities][points][dofs]
    static const double FE1_C1_D01_F_Q4a8[1][4][2][4] = {
        { { { -0.7886751345948126, -0.2113248654051871, 0.7886751345948129, 0.2113248654051871 },
            { -0.2113248654051872, -0.7886751345948129, 0.2113248654051871, 0.7886751345948129 } },
          { { -1.0, 0.0, 1.0, 0.0 }, { -1.0, 0.0, 1.0, 0.0 } },
          { { 0.0, -1.0, 0.0, 1.0 }, { 0.0, -1.0, 0.0, 1.0 } },
          { { -0.7886751345948126, -0.2113248654051871, 0.7886751345948129, 0.2113248654051871 },
            { -0.2113248654051872, -0.7886751345948129, 0.2113248654051871, 0.7886751345948129 } } }
    };
    static const double FE5_C0_F_Q4a8[1][4][2][9] = {
        { { { 0.4553418012614795, -0.1220084679281462, 0.0, 0.0, 0.6666666666666667, 0.0, 0.0, 0.0,
              0.0 },
            { -0.1220084679281463, 0.4553418012614796, 0.0, 0.0, 0.6666666666666667, 0.0, 0.0, 0.0,
              0.0 } },
          { { 0.4553418012614795, 0.0, -0.1220084679281462, 0.0, 0.0, 0.6666666666666667, 0.0, 0.0,
              0.0 },
            { -0.1220084679281462, 0.0, 0.4553418012614796, 0.0, 0.0, 0.6666666666666667, 0.0, 0.0,
              0.0 } },
          { { 0.0, 0.4553418012614796, 0.0, -0.1220084679281462, 0.0, 0.0, 0.6666666666666667, 0.0,
              0.0 },
            { 0.0, -0.1220084679281461, 0.0, 0.4553418012614795, 0.0, 0.0, 0.6666666666666667, 0.0,
              0.0 } },
          { { 0.0, 0.0, 0.4553418012614796, -0.1220084679281462, 0.0, 0.0, 0.0, 0.6666666666666667,
              0.0 },
            { 0.0, 0.0, -0.1220084679281461, 0.4553418012614794, 0.0, 0.0, 0.0, 0.6666666666666667,
              0.0 } } }
    };
    static const double FE5_C2_D10_F_Q4a8[1][4][2][4] = {
        { { { -1.0, 1.0, 0.0, 0.0 }, { -1.0, 1.0, 0.0, 0.0 } },
          { { -0.7886751345948126, 0.7886751345948129, -0.2113248654051871, 0.2113248654051871 },
            { -0.2113248654051872, 0.2113248654051871, -0.7886751345948129, 0.7886751345948129 } },
          { { -0.7886751345948126, 0.7886751345948129, -0.2113248654051871, 0.2113248654051871 },
            { -0.2113248654051872, 0.2113248654051871, -0.7886751345948129, 0.7886751345948129 } },
          { { 0.0, 0.0, -1.0, 1.0 }, { 0.0, 0.0, -1.0, 1.0 } } }
    };
    static const double FE5_C3_F_Q4a8[1][4][2][4] = {
        { { { 1.0, 0.0, 0.0, 0.0 }, { 1.0, 0.0, 0.0, 0.0 } },
          { { 0.7886751345948129, 0.0, 0.0, 0.2113248654051871 },
            { 0.2113248654051872, 0.0, 0.0, 0.7886751345948129 } },
          { { 0.7886751345948129, 0.0, 0.0, 0.2113248654051871 },
            { 0.2113248654051872, 0.0, 0.0, 0.7886751345948129 } },
          { { 0.0, 0.0, 0.0, 1.0 }, { 0.0, 0.0, 0.0, 1.0 } } }
    };
    static const double FE5_C4_F_Q4a8[1][4][2][4] = {
        { { { 0.0, 0.7886751345948129, 0.2113248654051871, 0.0 },
            { 0.0, 0.2113248654051872, 0.7886751345948129, 0.0 } },
          { { 0.0, 1.0, 0.0, 0.0 }, { 0.0, 1.0, 0.0, 0.0 } },
          { { 0.0, 0.0, 1.0, 0.0 }, { 0.0, 0.0, 1.0, 0.0 } },
          { { 0.0, 0.7886751345948129, 0.2113248654051871, 0.0 },
            { 0.0, 0.2113248654051872, 0.7886751345948129, 0.0 } } }
    };
    static const double quadrilateral_reference_facet_jacobian[4][2][1] = {
        { { 1.0 }, { 0.0 } }, { { 0.0 }, { 1.0 } }, { { 0.0 }, { 1.0 } }, { { 1.0 }, { 0.0 } }
    };
    static const double quadrilateral_reference_facet_normals[4][2] = {
        { 0.0, -1.0 }, { -1.0, -0.0 }, { 1.0, 0.0 }, { -0.0, 1.0 }
    };
    for ( int iq = 0; iq < 2; ++iq ) {
        // ------------------------
        // Section: Jacobian
        // Inputs: FE5_C2_D10_F_Q4a8, FE1_C1_D01_F_Q4a8, coordinate_dofs
        // Outputs: J_c2, J_c1, J_c0, J_c3
        double J_c3 = 0.0;
        double J_c0 = 0.0;
        double J_c1 = 0.0;
        double J_c2 = 0.0;
        {
            for ( int ic = 0; ic < 4; ++ic ) {
                J_c3 += coordinate_dofs[(ic)*3 + 1] *
                        FE1_C1_D01_F_Q4a8[0][entity_local_index[0]][iq][ic];
                J_c0 +=
                    coordinate_dofs[(ic)*3] * FE5_C2_D10_F_Q4a8[0][entity_local_index[0]][iq][ic];
                J_c1 +=
                    coordinate_dofs[(ic)*3] * FE1_C1_D01_F_Q4a8[0][entity_local_index[0]][iq][ic];
                J_c2 += coordinate_dofs[(ic)*3 + 1] *
                        FE5_C2_D10_F_Q4a8[0][entity_local_index[0]][iq][ic];
            }
        }
        // ------------------------
        // ------------------------
        // Section: Intermediates
        // Inputs: J_c2, J_c1, J_c0, J_c3
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
            double sv_4a8_0 = J_c0 * J_c3;
            double sv_4a8_1 = J_c1 * J_c2;
            double sv_4a8_2 = -sv_4a8_1;
            double sv_4a8_3 = sv_4a8_0 + sv_4a8_2;
            double sv_4a8_4 = J_c3 / sv_4a8_3;
            double sv_4a8_5 = -J_c2;
            double sv_4a8_6 = sv_4a8_5 / sv_4a8_3;
            double sv_4a8_7 = -sv_4a8_4;
            double sv_4a8_8 = -sv_4a8_6;
            double sv_4a8_9 = J_c0 / sv_4a8_3;
            double sv_4a8_10 =
                quadrilateral_reference_facet_normals[entity_local_index[0]][1] * sv_4a8_9;
            double sv_4a8_11 = -J_c1;
            double sv_4a8_12 = sv_4a8_11 / sv_4a8_3;
            double sv_4a8_13 =
                quadrilateral_reference_facet_normals[entity_local_index[0]][0] * sv_4a8_12;
            double sv_4a8_14 = sv_4a8_10 + sv_4a8_13;
            double sv_4a8_15 = sv_4a8_14 * sv_4a8_14;
            double sv_4a8_16 =
                quadrilateral_reference_facet_normals[entity_local_index[0]][0] * sv_4a8_4;
            double sv_4a8_17 =
                quadrilateral_reference_facet_normals[entity_local_index[0]][1] * sv_4a8_6;
            double sv_4a8_18 = sv_4a8_16 + sv_4a8_17;
            double sv_4a8_19 = sv_4a8_18 * sv_4a8_18;
            double sv_4a8_20 = sv_4a8_15 + sv_4a8_19;
            double sv_4a8_21 = sqrt( sv_4a8_20 );
            double sv_4a8_22 = sv_4a8_14 / sv_4a8_21;
            double sv_4a8_23 = -sv_4a8_22;
            double sv_4a8_24 = sv_4a8_23 * sv_4a8_4;
            double sv_4a8_25 = sv_4a8_23 * sv_4a8_6;
            double sv_4a8_26 = -sv_4a8_23;
            double sv_4a8_27 = sv_4a8_7 * sv_4a8_23;
            double sv_4a8_28 = sv_4a8_8 * sv_4a8_23;
            double sv_4a8_29 = -sv_4a8_12;
            double sv_4a8_30 = -sv_4a8_9;
            double sv_4a8_31 = sv_4a8_18 / sv_4a8_21;
            double sv_4a8_32 = sv_4a8_12 * sv_4a8_31;
            double sv_4a8_33 = sv_4a8_9 * sv_4a8_31;
            double sv_4a8_34 = sv_4a8_29 * sv_4a8_31;
            double sv_4a8_35 = sv_4a8_30 * sv_4a8_31;
            double sv_4a8_36 = -sv_4a8_31;
            double sv_4a8_37 = sv_4a8_24 + sv_4a8_32;
            double sv_4a8_38 = sv_4a8_25 + sv_4a8_33;
            double sv_4a8_39 = sv_4a8_27 + sv_4a8_34;
            double sv_4a8_40 = sv_4a8_28 + sv_4a8_35;
            double sv_4a8_41 = sv_4a8_37 * sv_4a8_37;
            double sv_4a8_42 = sv_4a8_38 * sv_4a8_37;
            double sv_4a8_43 = sv_4a8_38 * sv_4a8_38;
            double sv_4a8_44 = sv_4a8_37 * sv_4a8_26;
            double sv_4a8_45 = sv_4a8_38 * sv_4a8_26;
            double sv_4a8_46 = sv_4a8_39 * sv_4a8_37;
            double sv_4a8_47 = sv_4a8_39 * sv_4a8_38;
            double sv_4a8_48 = sv_4a8_40 * sv_4a8_37;
            double sv_4a8_49 = sv_4a8_40 * sv_4a8_38;
            double sv_4a8_50 = sv_4a8_37 * sv_4a8_36;
            double sv_4a8_51 = sv_4a8_38 * sv_4a8_36;
            double sv_4a8_52 =
                J_c0 * quadrilateral_reference_facet_jacobian[entity_local_index[0]][0][0];
            double sv_4a8_53 =
                J_c1 * quadrilateral_reference_facet_jacobian[entity_local_index[0]][1][0];
            double sv_4a8_54 = sv_4a8_52 + sv_4a8_53;
            double sv_4a8_55 = sv_4a8_54 * sv_4a8_54;
            double sv_4a8_56 =
                quadrilateral_reference_facet_jacobian[entity_local_index[0]][0][0] * J_c2;
            double sv_4a8_57 =
                quadrilateral_reference_facet_jacobian[entity_local_index[0]][1][0] * J_c3;
            double sv_4a8_58 = sv_4a8_56 + sv_4a8_57;
            double sv_4a8_59 = sv_4a8_58 * sv_4a8_58;
            double sv_4a8_60 = sv_4a8_55 + sv_4a8_59;
            double sv_4a8_61 = sqrt( sv_4a8_60 );
            double sv_4a8_62 = sv_4a8_41 * sv_4a8_61;
            double sv_4a8_63 = sv_4a8_42 * sv_4a8_61;
            double sv_4a8_64 = sv_4a8_43 * sv_4a8_61;
            double sv_4a8_65 = sv_4a8_44 * sv_4a8_61;
            double sv_4a8_66 = sv_4a8_45 * sv_4a8_61;
            double sv_4a8_67 = sv_4a8_46 * sv_4a8_61;
            double sv_4a8_68 = sv_4a8_47 * sv_4a8_61;
            double sv_4a8_69 = sv_4a8_48 * sv_4a8_61;
            double sv_4a8_70 = sv_4a8_49 * sv_4a8_61;
            double sv_4a8_71 = sv_4a8_50 * sv_4a8_61;
            double sv_4a8_72 = sv_4a8_51 * sv_4a8_61;
            fw0 = sv_4a8_65 * weights_4a8[iq];
            fw1 = sv_4a8_66 * weights_4a8[iq];
            fw2 = sv_4a8_71 * weights_4a8[iq];
            fw3 = sv_4a8_72 * weights_4a8[iq];
            fw4 = sv_4a8_62 * weights_4a8[iq];
            fw5 = sv_4a8_63 * weights_4a8[iq];
            fw6 = sv_4a8_64 * weights_4a8[iq];
            fw7 = sv_4a8_67 * weights_4a8[iq];
            fw8 = sv_4a8_68 * weights_4a8[iq];
            fw9 = sv_4a8_69 * weights_4a8[iq];
            fw10 = sv_4a8_70 * weights_4a8[iq];
        }
        // ------------------------
        // ------------------------
        // Section: Tensor Computation
        // Inputs: FE5_C3_F_Q4a8, fw2, fw3, FE5_C0_F_Q4a8, fw1, FE5_C4_F_Q4a8, fw0
        // Outputs: A
        {
            double temp_0[4] = { 0 };
            for ( int j = 0; j < 4; ++j ) {
                temp_0[j] = fw0 * FE5_C3_F_Q4a8[0][entity_local_index[0]][iq][j];
            }
            double temp_1[4] = { 0 };
            for ( int j = 0; j < 4; ++j ) {
                temp_1[j] = fw1 * FE5_C4_F_Q4a8[0][entity_local_index[0]][iq][j];
            }
            double temp_2[4] = { 0 };
            for ( int j = 0; j < 4; ++j ) {
                temp_2[j] = fw2 * FE5_C3_F_Q4a8[0][entity_local_index[0]][iq][j];
            }
            double temp_3[4] = { 0 };
            for ( int j = 0; j < 4; ++j ) {
                temp_3[j] = fw3 * FE5_C4_F_Q4a8[0][entity_local_index[0]][iq][j];
            }
            for ( int j = 0; j < 4; ++j ) {
                for ( int i = 0; i < 9; ++i ) {
                    A[30 * ( 2 * ( i ) ) + ( ( j ) + 26 )] +=
                        FE5_C0_F_Q4a8[0][entity_local_index[0]][iq][i] * temp_0[j];
                    A[30 * ( 2 * ( i ) ) + ( ( j ) + 26 )] +=
                        FE5_C0_F_Q4a8[0][entity_local_index[0]][iq][i] * temp_1[j];
                    A[30 * ( 2 * ( i ) + 1 ) + ( ( j ) + 26 )] +=
                        FE5_C0_F_Q4a8[0][entity_local_index[0]][iq][i] * temp_2[j];
                    A[30 * ( 2 * ( i ) + 1 ) + ( ( j ) + 26 )] +=
                        FE5_C0_F_Q4a8[0][entity_local_index[0]][iq][i] * temp_3[j];
                }
            }
        }
        // ------------------------
        // ------------------------
        // Section: Tensor Computation
        // Inputs: FE5_C2_D10_F_Q4a8, fw4, FE5_C3_F_Q4a8, fw5, FE1_C1_D01_F_Q4a8, FE5_C4_F_Q4a8,
        // fw7, fw8, fw10, fw9, fw6 Outputs: A
        {
            double temp_0[4] = { 0 };
            for ( int j = 0; j < 4; ++j ) {
                temp_0[j] = fw4 * FE5_C3_F_Q4a8[0][entity_local_index[0]][iq][j];
            }
            double temp_1[4] = { 0 };
            for ( int j = 0; j < 4; ++j ) {
                temp_1[j] = fw5 * FE5_C4_F_Q4a8[0][entity_local_index[0]][iq][j];
            }
            double temp_2[4] = { 0 };
            for ( int j = 0; j < 4; ++j ) {
                temp_2[j] = fw5 * FE5_C3_F_Q4a8[0][entity_local_index[0]][iq][j];
            }
            double temp_3[4] = { 0 };
            for ( int j = 0; j < 4; ++j ) {
                temp_3[j] = fw6 * FE5_C4_F_Q4a8[0][entity_local_index[0]][iq][j];
            }
            double temp_4[4] = { 0 };
            for ( int j = 0; j < 4; ++j ) {
                temp_4[j] = fw7 * FE5_C3_F_Q4a8[0][entity_local_index[0]][iq][j];
            }
            double temp_5[4] = { 0 };
            for ( int j = 0; j < 4; ++j ) {
                temp_5[j] = fw8 * FE5_C4_F_Q4a8[0][entity_local_index[0]][iq][j];
            }
            double temp_6[4] = { 0 };
            for ( int j = 0; j < 4; ++j ) {
                temp_6[j] = fw9 * FE5_C3_F_Q4a8[0][entity_local_index[0]][iq][j];
            }
            double temp_7[4] = { 0 };
            for ( int j = 0; j < 4; ++j ) {
                temp_7[j] = fw10 * FE5_C4_F_Q4a8[0][entity_local_index[0]][iq][j];
            }
            double temp_8[4] = { 0 };
            for ( int j = 0; j < 4; ++j ) {
                temp_8[j] = fw4 * FE5_C2_D10_F_Q4a8[0][entity_local_index[0]][iq][j];
            }
            double temp_9[4] = { 0 };
            for ( int j = 0; j < 4; ++j ) {
                temp_9[j] = fw5 * FE1_C1_D01_F_Q4a8[0][entity_local_index[0]][iq][j];
            }
            double temp_10[4] = { 0 };
            for ( int j = 0; j < 4; ++j ) {
                temp_10[j] = fw5 * FE5_C2_D10_F_Q4a8[0][entity_local_index[0]][iq][j];
            }
            double temp_11[4] = { 0 };
            for ( int j = 0; j < 4; ++j ) {
                temp_11[j] = fw6 * FE1_C1_D01_F_Q4a8[0][entity_local_index[0]][iq][j];
            }
            double temp_12[4] = { 0 };
            for ( int j = 0; j < 4; ++j ) {
                temp_12[j] = fw7 * FE5_C3_F_Q4a8[0][entity_local_index[0]][iq][j];
            }
            double temp_13[4] = { 0 };
            for ( int j = 0; j < 4; ++j ) {
                temp_13[j] = fw9 * FE5_C4_F_Q4a8[0][entity_local_index[0]][iq][j];
            }
            double temp_14[4] = { 0 };
            for ( int j = 0; j < 4; ++j ) {
                temp_14[j] = fw8 * FE5_C3_F_Q4a8[0][entity_local_index[0]][iq][j];
            }
            double temp_15[4] = { 0 };
            for ( int j = 0; j < 4; ++j ) {
                temp_15[j] = fw10 * FE5_C4_F_Q4a8[0][entity_local_index[0]][iq][j];
            }
            for ( int j = 0; j < 4; ++j ) {
                for ( int i = 0; i < 4; ++i ) {
                    A[30 * ( ( i ) + 18 ) + ( ( j ) + 26 )] +=
                        FE5_C2_D10_F_Q4a8[0][entity_local_index[0]][iq][i] * temp_0[j];
                    A[30 * ( ( i ) + 18 ) + ( ( j ) + 26 )] +=
                        FE5_C2_D10_F_Q4a8[0][entity_local_index[0]][iq][i] * temp_1[j];
                    A[30 * ( ( i ) + 18 ) + ( ( j ) + 26 )] +=
                        FE1_C1_D01_F_Q4a8[0][entity_local_index[0]][iq][i] * temp_2[j];
                    A[30 * ( ( i ) + 18 ) + ( ( j ) + 26 )] +=
                        FE1_C1_D01_F_Q4a8[0][entity_local_index[0]][iq][i] * temp_3[j];
                    A[30 * ( ( i ) + 22 ) + ( ( j ) + 26 )] +=
                        FE5_C3_F_Q4a8[0][entity_local_index[0]][iq][i] * temp_4[j];
                    A[30 * ( ( i ) + 22 ) + ( ( j ) + 26 )] +=
                        FE5_C3_F_Q4a8[0][entity_local_index[0]][iq][i] * temp_5[j];
                    A[30 * ( ( i ) + 22 ) + ( ( j ) + 26 )] +=
                        FE5_C4_F_Q4a8[0][entity_local_index[0]][iq][i] * temp_6[j];
                    A[30 * ( ( i ) + 22 ) + ( ( j ) + 26 )] +=
                        FE5_C4_F_Q4a8[0][entity_local_index[0]][iq][i] * temp_7[j];
                    A[30 * ( ( i ) + 26 ) + ( ( j ) + 18 )] +=
                        FE5_C3_F_Q4a8[0][entity_local_index[0]][iq][i] * temp_8[j];
                    A[30 * ( ( i ) + 26 ) + ( ( j ) + 18 )] +=
                        FE5_C3_F_Q4a8[0][entity_local_index[0]][iq][i] * temp_9[j];
                    A[30 * ( ( i ) + 26 ) + ( ( j ) + 18 )] +=
                        FE5_C4_F_Q4a8[0][entity_local_index[0]][iq][i] * temp_10[j];
                    A[30 * ( ( i ) + 26 ) + ( ( j ) + 18 )] +=
                        FE5_C4_F_Q4a8[0][entity_local_index[0]][iq][i] * temp_11[j];
                    A[30 * ( ( i ) + 26 ) + ( ( j ) + 22 )] +=
                        FE5_C3_F_Q4a8[0][entity_local_index[0]][iq][i] * temp_12[j];
                    A[30 * ( ( i ) + 26 ) + ( ( j ) + 22 )] +=
                        FE5_C3_F_Q4a8[0][entity_local_index[0]][iq][i] * temp_13[j];
                    A[30 * ( ( i ) + 26 ) + ( ( j ) + 22 )] +=
                        FE5_C4_F_Q4a8[0][entity_local_index[0]][iq][i] * temp_14[j];
                    A[30 * ( ( i ) + 26 ) + ( ( j ) + 22 )] +=
                        FE5_C4_F_Q4a8[0][entity_local_index[0]][iq][i] * temp_15[j];
                }
            }
        }
        // ------------------------
        // ------------------------
        // Section: Tensor Computation
        // Inputs: FE5_C3_F_Q4a8, fw2, fw3, FE5_C0_F_Q4a8, fw1, FE5_C4_F_Q4a8, fw0
        // Outputs: A
        {
            double temp_0[9] = { 0 };
            for ( int j = 0; j < 9; ++j ) {
                temp_0[j] = fw0 * FE5_C0_F_Q4a8[0][entity_local_index[0]][iq][j];
            }
            double temp_1[9] = { 0 };
            for ( int j = 0; j < 9; ++j ) {
                temp_1[j] = fw1 * FE5_C0_F_Q4a8[0][entity_local_index[0]][iq][j];
            }
            double temp_2[9] = { 0 };
            for ( int j = 0; j < 9; ++j ) {
                temp_2[j] = fw2 * FE5_C0_F_Q4a8[0][entity_local_index[0]][iq][j];
            }
            double temp_3[9] = { 0 };
            for ( int j = 0; j < 9; ++j ) {
                temp_3[j] = fw3 * FE5_C0_F_Q4a8[0][entity_local_index[0]][iq][j];
            }
            for ( int j = 0; j < 9; ++j ) {
                for ( int i = 0; i < 4; ++i ) {
                    A[30 * ( ( i ) + 26 ) + 2 * ( j )] +=
                        FE5_C3_F_Q4a8[0][entity_local_index[0]][iq][i] * temp_0[j];
                    A[30 * ( ( i ) + 26 ) + 2 * ( j )] +=
                        FE5_C4_F_Q4a8[0][entity_local_index[0]][iq][i] * temp_1[j];
                    A[30 * ( ( i ) + 26 ) + ( 2 * ( j ) + 1 )] +=
                        FE5_C3_F_Q4a8[0][entity_local_index[0]][iq][i] * temp_2[j];
                    A[30 * ( ( i ) + 26 ) + ( 2 * ( j ) + 1 )] +=
                        FE5_C4_F_Q4a8[0][entity_local_index[0]][iq][i] * temp_3[j];
                }
            }
        }
        // ------------------------
    }
    return A;
}

extern "C" {
void BP1_qu9_Fortran( const double *w, const int nw, const double *coordinate_dofs, const int ncd,
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
    VectorReal A_vec = B_p1_qu9( w_vec, coordinate_dofs_vec, entities_vec, cst_vec );

    // Copie du résultat dans le tableau C fourni par Fortran
    for ( size_t i = 0; i < A_vec.size(); ++i ) {
        A[i] = A_vec[i];
    }
}
void BP2_qu9_Fortran( const double *w, const int nw, const double *coordinate_dofs, const int ncd,
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
    VectorReal A_vec = B_p2_qu9( w_vec, coordinate_dofs_vec, entities_vec, cst_vec );

    // Copie du résultat dans le tableau C fourni par Fortran
    for ( size_t i = 0; i < A_vec.size(); ++i ) {
        A[i] = A_vec[i];
    }
}
void BP4_qu9_Fortran( const double *w, const int nw, const double *coordinate_dofs, const int ncd,
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
    VectorReal A_vec = B_p4_qu9( w_vec, coordinate_dofs_vec, entities_vec, cst_vec );

    // Copie du résultat dans le tableau C fourni par Fortran
    for ( size_t i = 0; i < A_vec.size(); ++i ) {
        A[i] = A_vec[i];
    }
}
void BP5_qu9_Fortran( const double *w, const int nw, const double *coordinate_dofs, const int ncd,
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
    VectorReal A_vec = B_p5_qu9( w_vec, coordinate_dofs_vec, entities_vec, cst_vec );

    // Copie du résultat dans le tableau C fourni par Fortran
    for ( size_t i = 0; i < A_vec.size(); ++i ) {
        A[i] = A_vec[i];
    }
}
}
