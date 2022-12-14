! --------------------------------------------------------------------
! Copyright (C) 1991 - 2022 - EDF R&D - www.code-aster.org
! This file is part of code_aster.
!
! code_aster is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! code_aster is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with code_aster.  If not, see <http://www.gnu.org/licenses/>.
! --------------------------------------------------------------------

subroutine load_neut_data(i_type_neum    , nb_type_neumz, type_calc_,&
                          load_type_ligr_, load_opti_r_ , load_opti_f_, load_para_r_,&
                          load_para_f_   , load_keyw_   , load_obje_  , nb_obje_)
!
implicit none
!
#include "asterfort/assert.h"
!
! person_in_charge: mickael.abbas at edf.fr
!
    integer, intent(in) :: i_type_neum
    integer, intent(in) :: nb_type_neumz
    character(len=4), optional, intent(in) :: type_calc_
    character(len=6), optional, intent(out) :: load_type_ligr_
    character(len=16), optional, intent(out) :: load_opti_r_
    character(len=16), optional, intent(out) :: load_opti_f_
    character(len=8), optional, intent(out) :: load_para_r_(2)
    character(len=8), optional, intent(out) :: load_para_f_(2)
    character(len=24), optional, intent(out) :: load_keyw_
    character(len=10), optional, intent(out) :: load_obje_(2)
    integer, optional, intent(out) :: nb_obje_
!
! --------------------------------------------------------------------------------------------------
!
! Neumann loads computation - Thermic
!
! Get information about load (Neumann)
!
! --------------------------------------------------------------------------------------------------
!
! In  i_type_neum      : index for Neumann load type
! In  nb_type_neumz    : maximum number of Neumann load type
! In  type_calc        : type of option to compute
!                        '2MBR' for second member (vector)
!                        'RESI' for residual (vector)
!                        'MRIG' for rigidity (matrix)
!                        'MTAN' for tangent matrix
! Out load_type_ligr   : type of LIGREL for current load
! Out load_opti_r      : option for real parameter
! Out load_opti_f      : option for function parameter
! Out load_para_r      : name of parameterS (real)
! Out load_para_f      : name of parameterS (function)
! Out load_keyw        : keyword in AFFE_CHAR_THER for this load
! Out load_obje        : name of objectS (cart in AFFE_CHAR_THER)
! Out nb_obje          : number of objects
!
! --------------------------------------------------------------------------------------------------
!
    integer :: nb_type_neum
    parameter (nb_type_neum = 11)
    character(len=10) :: object1(nb_type_neum)
    character(len=10) :: object2(nb_type_neum)
    character(len=6) :: ligrel(nb_type_neum)
    character(len=24) :: keyw(nb_type_neum)
    character(len=8) :: para_2mbr_r1(nb_type_neum), para_2mbr_f1(nb_type_neum)
    character(len=8) :: para_2mbr_r2(nb_type_neum), para_2mbr_f2(nb_type_neum)
    character(len=16) :: opti_2mbr_f(nb_type_neum), opti_2mbr_r(nb_type_neum)
    character(len=8) :: para_resi_r1(nb_type_neum), para_resi_f1(nb_type_neum)
    character(len=8) :: para_resi_r2(nb_type_neum), para_resi_f2(nb_type_neum)
    character(len=16) :: opti_resi_f(nb_type_neum), opti_resi_r(nb_type_neum)
    character(len=8) :: para_mrig_r(nb_type_neum), para_mrig_f(nb_type_neum)
    character(len=16) :: opti_mrig_f(nb_type_neum), opti_mrig_r(nb_type_neum)
    character(len=8) :: para_mtan_r(nb_type_neum), para_mtan_f(nb_type_neum)
    character(len=16) :: opti_mtan_f(nb_type_neum), opti_mtan_r(nb_type_neum)
!
! - Keyword in AFFE_CHAR_THER
!
    data keyw     /'ECHANGE'      ,'FLUX_REP_XYZ' ,'FLUX_REP_NORM','SOURCE' ,&
                   'ECHANGE_PAROI','PRE_GRAD_TEMP','FLUX_NL'      ,'SOUR_NL',&
                   'RAYONNEMENT'  ,'EVOL_CHAR'    ,'SOURCE_CALCULEE'/
!
! - Object name construct in AFFE_CHAR_THER
!
    data object1  /'.COEFH','.FLURE','.FLUR2','.SOURE',&
                   '.HECHP','.GRAIN','.FLUNL','.SOUNL',&
                   '.RAYO ','.EVOL.CHAR','.SOURC'/
    data object2  /'.T_EXT',' '     ,' '     ,' '     ,&
                   ' '     ,' '     ,' '     ,' '     ,&
                   ' '     ,' '     ,' '/
!
! - Type of LIGREL to compute
!
    data ligrel   /'Model','Model','Model','Model',&
                   'Load' ,'Model','Model','Model',&
                   'Model','Model','Model'/
!
! - Name of option for second member (real coefficient)
!
    data opti_2mbr_r /'CHAR_THER_TEXT_R','CHAR_THER_FLUN_R','CHAR_THER_FLUX_R','CHAR_THER_SOUR_R',&
                      'CHAR_THER_PARO_R','CHAR_THER_GRAI_R','No_Load'         ,'No_Load'         ,&
                      'No_Load'         ,'CHAR_EVOL_CHAR'  ,'CHAR_THER_SOUR_R'/
!
! - Name of option for second member (function coefficient)
!
    data opti_2mbr_f /'CHAR_THER_TEXT_F','CHAR_THER_FLUN_F','CHAR_THER_FLUX_F','CHAR_THER_SOUR_F',&
                      'CHAR_THER_PARO_F','CHAR_THER_GRAI_F','No_Load'         ,'No_Load'         ,&
                      'No_Load'         ,'No_Load'         ,'No_Load'/
!
! - Name of input parameter field for second member  (real coefficient)
!
    data para_2mbr_r1 /'PCOEFHR','PFLUXNR','PFLUXVR','PSOURCR',&
                       'PHECHPR','PGRAINR','       ','       ',&
                       '       ','       ','PSOURCR'/
!
! - Name of input parameter field for second member  (function coefficient)
!
    data para_2mbr_f1 /'PCOEFHF','PFLUXNF','PFLUXVF','PSOURCF',&
                       'PHECHPF','PGRAINF','       ','       ',&
                       '       ','       ','       '/
!
! - Name of input parameter field for second member  (real coefficient)
!
    data para_2mbr_r2 /'PT_EXTR','       ','       ','       ',&
                       '       ','       ','       ','       ',&
                       '       ','       ','       '/
!
! - Name of input parameter field for second member  (function coefficient)
!
    data para_2mbr_f2 /'PT_EXTF','       ','       ','       ',&
                       '       ','       ','       ','       ',&
                       '       ','       ','       '/
!
! - Name of option for residual (real coefficient)
!
    data opti_resi_r /'RESI_THER_COEF_R','No_Load'       ,'No_Load'         ,'No_Load'         ,&
                      'RESI_THER_PARO_R','No_Load'       ,'RESI_THER_FLUXNL','RESI_THER_SOURNL',&
                      'RESI_THER_RAYO_R','RESI_EVOL_CHAR','No_Load'/
!
! - Name of option for residual (function coefficient)
!
    data opti_resi_f /'RESI_THER_COEF_F','No_Load'       ,'No_Load'         ,'No_Load',&
                      'RESI_THER_PARO_F','No_Load'       ,'RESI_THER_FLUXNL','RESI_THER_SOURNL',&
                      'RESI_THER_RAYO_F','No_Load'       ,'No_Load'/
!
! - Name of input parameter field for residual (real coefficient)
!
    data para_resi_r1 /'PCOEFHR','PFLUXNR','PFLUXVR','PSOURCR',&
                       'PHECHPR','PGRAINR','PFLUXNL','PSOURNL',&
                       'PRAYONR','       ','PSOURCR'/
!
! - Name of input parameter field for residual (function coefficient)
!
    data para_resi_f1 /'PCOEFHF','PFLUXNF','PFLUXVF','PSOURCF',&
                       'PHECHPF','PGRAINF','PFLUXNL','PSOURNL',&
                       'PRAYONF','       ','       '/
!
! - Name of input parameter field for second member  (real coefficient)
!
    data para_resi_r2 /'PT_EXTR','       ','       ','       ',&
                       '       ','       ','       ','       ',&
                       '       ','       ','       '/
!
! - Name of input parameter field for second member  (function coefficient)
!
    data para_resi_f2 /'PT_EXTF','       ','       ','       ',&
                       '       ','       ','       ','       ',&
                       '       ','       ','       '/
!
! - Name of option for rigidity matrix (real coefficient)
!
    data opti_mrig_r /'RIGI_THER_COEH_R','No_Load','No_Load','No_Load',&
                      'RIGI_THER_PARO_R','No_Load','No_Load','No_Load',&
                      'No_Load'         ,'No_Load','No_Load'/
!
! - Name of option for rigidity matrix (function coefficient)
!
    data opti_mrig_f /'RIGI_THER_COEH_F','No_Load','No_Load','No_Load',&
                      'RIGI_THER_PARO_F','No_Load','No_Load','No_Load',&
                      'No_Load'         ,'No_Load','No_Load'/
!
! - Name of input parameter field for rigidity matrix  (real coefficient)
!
    data para_mrig_r /'PCOEFHR','       ','       ','       ',&
                      'PHECHPR','       ','       ','       ',&
                      '       ','       ','       '/
!
! - Name of input parameter field for rigidity matrix  (function coefficient)
!
    data para_mrig_f /'PCOEFHF','       ','       ','       ',&
                      'PHECHPF','       ','       ','       ',&
                      '       ','       ','       '/
!
! - Name of option for tangent matrix (real coefficient)
!
    data opti_mtan_r /'MTAN_THER_COEF_R','No_Load','No_Load','No_Load',&
                      'MTAN_THER_PARO_R','No_Load','No_Load','No_Load',&
                      'MTAN_THER_RAYO_R','No_Load','No_Load'/
!
! - Name of option for tangent matrix (function coefficient)
!
    data opti_mtan_f /'MTAN_THER_COEF_F','No_Load','No_Load'         ,'No_Load'         ,&
                      'MTAN_THER_PARO_F','No_Load','MTAN_THER_FLUXNL','MTAN_THER_SOURNL',&
                      'MTAN_THER_RAYO_F','No_Load','No_Load'/
!
! - Name of input parameter field for tangent matrix  (real coefficient)
!
    data para_mtan_r /'PCOEFHR','       ','       ','       ',&
                      'PHECHPR','       ','       ','       ',&
                      'PRAYONR','       ','       '/
!
! - Name of input parameter field for tangent matrix  (function coefficient)
!
    data para_mtan_f /'PCOEFHF','       ','       ','       ',&
                      'PHECHPF','       ','PFLUXNL','PSOURNL',&
                      'PRAYONF','       ','       '/
!
! --------------------------------------------------------------------------------------------------
!
    ASSERT(i_type_neum.le.nb_type_neum)
    ASSERT(nb_type_neumz.eq.nb_type_neum)

    if (present(load_type_ligr_)) then
        load_type_ligr_ = ligrel(i_type_neum)
    endif
    if (present(load_opti_r_)) then
        if (type_calc_.eq.'2MBR') then
            load_opti_r_ = opti_2mbr_r(i_type_neum)
        elseif (type_calc_.eq.'RESI') then
            load_opti_r_ = opti_resi_r(i_type_neum)
        elseif (type_calc_.eq.'MRIG') then
            load_opti_r_ = opti_mrig_r(i_type_neum)
        elseif (type_calc_.eq.'MTAN') then
            load_opti_r_ = opti_mtan_r(i_type_neum)
        else
            ASSERT(.false.)
        endif
    endif
    if (present(load_opti_f_)) then
        if (type_calc_.eq.'2MBR') then
            load_opti_f_ = opti_2mbr_f(i_type_neum)
        elseif (type_calc_.eq.'RESI') then
            load_opti_f_ = opti_resi_f(i_type_neum)
        elseif (type_calc_.eq.'MRIG') then
            load_opti_f_ = opti_mrig_f(i_type_neum)
        elseif (type_calc_.eq.'MTAN') then
            load_opti_f_ = opti_mtan_f(i_type_neum)
        else
            ASSERT(.false.)
        endif
    endif
    if (present(load_para_r_)) then
        if (type_calc_.eq.'2MBR') then
            load_para_r_(1) = para_2mbr_r1(i_type_neum)
            load_para_r_(2) = para_2mbr_r2(i_type_neum)
        elseif (type_calc_.eq.'RESI') then
            load_para_r_(1) = para_resi_r1(i_type_neum)
            load_para_r_(2) = para_resi_r2(i_type_neum)
        elseif (type_calc_.eq.'MRIG') then
            load_para_r_(1) = para_mrig_r(i_type_neum)
            load_para_r_(2) = ' '
        elseif (type_calc_.eq.'MTAN') then
            load_para_r_(1) = para_mtan_r(i_type_neum)
            load_para_r_(2) = ' '
        else
            ASSERT(.false.)
        endif
    endif
    if (present(load_para_f_)) then
        if (type_calc_.eq.'2MBR') then
            load_para_f_(1) = para_2mbr_f1(i_type_neum)
            load_para_f_(2) = para_2mbr_f2(i_type_neum)
        elseif (type_calc_.eq.'RESI') then
            load_para_f_(1) = para_resi_f1(i_type_neum)
            load_para_f_(2) = para_resi_f2(i_type_neum)
        elseif (type_calc_.eq.'MRIG') then
            load_para_f_(1) = para_mrig_f(i_type_neum)
            load_para_f_(2) = ' '
        elseif (type_calc_.eq.'MTAN') then
            load_para_f_(1) = para_mtan_f(i_type_neum)
            load_para_f_(2) = ' '
        else
            ASSERT(.false.)
        endif
    endif
    if (present(load_keyw_)) then
        load_keyw_ = keyw(i_type_neum)
    endif
    if (present(load_obje_)) then
        load_obje_(1) = object1(i_type_neum)
        load_obje_(2) = object2(i_type_neum)
    endif
    if (present(nb_obje_)) then
        nb_obje_ = 1
        if (load_obje_(2).ne.' ') then
            nb_obje_ = 2
        endif
    endif
!
end subroutine
