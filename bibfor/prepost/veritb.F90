! --------------------------------------------------------------------
! Copyright (C) 1991 - 2023 - EDF R&D - www.code-aster.org
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
!
subroutine veritb(nk1d, ndim, oridef, deklag, profil)
!
    implicit none
#include "asterf_types.h"
#include "asterfort/getvid.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/tbexp2.h"
#include "asterfort/tbexv1.h"
#include "asterfort/utmess.h"
#include "asterfort/verinr.h"
    integer :: nk1d, ndim
    real(kind=8) :: deklag
    character(len=8) :: oridef
    character(len=12) :: profil
! --- BUT : VERIFICATION DE LA PRESENCE DES CHAMPS NECESSAIRES ---------
! ======================================================================
! IN  : NK1D   : NOMBRE D'OCCURENCE DE K1D -----------------------------
! --- : NDIM   : DIMENSION DE L'ESPACE ---------------------------------
! --- : ORIDEF : TYPE D'ORIENTATION DU DEFAUT --------------------------
! ======================================================================
    aster_logical :: teste
    integer :: i, ibid, nbval1, nbval2, irev
    character(len=8) :: motfac, k8b, tabrev, tabmdb, tabthr
    character(len=19) :: tbins1, tbins2
! ======================================================================
    call jemarq()
! ======================================================================
! --- DEFINITION DES TABLES --------------------------------------------
! ======================================================================
    motfac = 'K1D'
    tbins1 = '&&VERITB.TBINS1'
    tbins2 = '&&VERITB.TBINS2'
! ======================================================================
! --- CAS DE LA PREMIERE OCCURENCE DE K1D ------------------------------
! ======================================================================
    call getvid(motfac, 'TABL_MECA_REV', iocc=1, scal=tabrev, nbret=irev)
    call getvid(motfac, 'TABL_MECA_MDB', iocc=1, scal=tabmdb, nbret=ibid)
    call getvid(motfac, 'TABL_THER', iocc=1, scal=tabthr, nbret=ibid)
    if (profil(1:7) .eq. 'ELLIPSE') then
        if (deklag .lt. 0.D0 .and. irev .eq. 0) call utmess('F', 'PREPOST_7')
    end if
    if (irev .eq. 0) tabrev = tabmdb
! ======================================================================
! --- VERIFICATION DE LA PRESENCE DE LISTE D'INSTANT -------------------
! ======================================================================
    call tbexp2(tabrev, 'INST')
    call tbexp2(tabmdb, 'INST')
    call tbexp2(tabthr, 'INST')
    call tbexp2(tabrev, 'ABSC_CURV')
    call tbexp2(tabthr, 'ABSC_CURV')
!
! ======================================================================
! --- VERIFICATION DE LA COHERENCE DES LISTES D'INSTANT POUR -----------
! --- LES CHAMPS MECANIQUES --------------------------------------------
! ======================================================================
    call tbexv1(tabrev, 'INST', tbins1, 'V', nbval1, &
                k8b)
    call tbexv1(tabmdb, 'INST', tbins2, 'V', nbval2, &
                k8b)
    if (nbval1 .ne. nbval2) then
        call utmess('F', 'PREPOST4_90')
    end if
    teste = verinr(nbval1, tbins1, tbins2)
    if (teste) then
        call utmess('F', 'PREPOST4_91')
    end if
! ======================================================================
! --- DESTRUCTIONS DES VECTEURS INUTILES -------------------------------
! ======================================================================
    call jedetr(tbins2)
    call jedetr(tabrev)
    call jedetr(tabmdb)
    call jedetr(tabthr)
! ======================================================================
! --- ITERATIONS SUR LES OCCURENCES DE K1D -----------------------------
! ======================================================================
    do i = 2, nk1d
! ======================================================================
! --- RECUPERATION DES TABLES ASSOCIEES A K1D POUR L'ITERATION COURANTE-
! ======================================================================
        call getvid(motfac, 'TABL_MECA_REV', iocc=i, scal=tabrev, nbret=ibid)
        call getvid(motfac, 'TABL_MECA_MDB', iocc=i, scal=tabmdb, nbret=ibid)
        call getvid(motfac, 'TABL_THER', iocc=i, scal=tabthr, nbret=ibid)
! ======================================================================
! --- VERIFICATION DE LA PRESENCE DE LISTE D'INSTANT -------------------
! ======================================================================
        call tbexp2(tabrev, 'INST')
        call tbexp2(tabmdb, 'INST')
        call tbexp2(tabthr, 'INST')
! ======================================================================
! --- VERIFICATION DE LA COHERENCE DES LISTES D'INSTANT POUR -----------
! --- LES CHAMPS MECANIQUES --------------------------------------------
! ======================================================================
        call tbexv1(tabrev, 'INST', tbins2, 'V', nbval2, &
                    k8b)
        if (nbval1 .ne. nbval2) then
            call utmess('F', 'PREPOST4_92')
        end if
        teste = verinr(nbval1, tbins1, tbins2)
        if (teste) then
            call utmess('F', 'PREPOST4_91')
        end if
        call jedetr(tbins2)
        call tbexv1(tabmdb, 'INST', tbins2, 'V', nbval2, &
                    k8b)
        if (nbval1 .ne. nbval2) then
            call utmess('F', 'PREPOST4_92')
        end if
        teste = verinr(nbval1, tbins1, tbins2)
        if (teste) then
            call utmess('F', 'PREPOST4_91')
        end if
        call jedetr(tbins2)
        call jedetr(tabrev)
        call jedetr(tabmdb)
        call jedetr(tabthr)
    end do
    call jedetr(tbins1)
! ======================================================================
! --- VERIFICATION DE LA PRESENCE DES BONNES COMPOSANTES POUR LE -------
! --- CALCUL DES FACTEURS D'INTENSITE DE CONTRAINTE --------------------
! ======================================================================
    do i = 1, nk1d
! ======================================================================
! --- RECUPERATION DES TABLES ASSOCIEES A K1D POUR L'ITERATION COURANTE-
! ======================================================================
        call getvid(motfac, 'TABL_MECA_REV', iocc=i, scal=tabrev, nbret=ibid)
        call getvid(motfac, 'TABL_MECA_MDB', iocc=i, scal=tabmdb, nbret=ibid)
        call getvid(motfac, 'TABL_THER', iocc=i, scal=tabthr, nbret=ibid)
        if (ndim .eq. 2) then
! ======================================================================
! --- CAS D'UNE DIMENSION D'ORDRE 2 ------------------------------------
! ======================================================================
            if (oridef .eq. 'CIRC') then
! ======================================================================
! --- CAS D'UN DEFAUT CIRCONFERENTIEL ----------------------------------
! ======================================================================
                call tbexp2(tabrev, 'SIYY')
                call tbexp2(tabmdb, 'SIYY')
            else
! ======================================================================
! --- CAS D'UN DEFAUT LONGITUDINAL -------------------------------------
! ======================================================================
                call tbexp2(tabrev, 'SIZZ')
                call tbexp2(tabmdb, 'SIZZ')
            end if
        else
! ======================================================================
! --- CAS D'UNE DIMENSION D'ORDRE 3 ------------------------------------
! ======================================================================
            if (oridef .eq. 'CIRC') then
! ======================================================================
! --- CAS D'UN DEFAUT CIRCONFERENTIEL ----------------------------------
! ======================================================================
                call tbexp2(tabrev, 'SIZZ')
                call tbexp2(tabmdb, 'SIZZ')
            else
! ======================================================================
! --- CAS D'UN DEFAUT LONGITUDINAL -------------------------------------
! ======================================================================
                call tbexp2(tabrev, 'SIXX')
                call tbexp2(tabmdb, 'SIXX')
                call tbexp2(tabrev, 'SIYY')
                call tbexp2(tabmdb, 'SIYY')
                call tbexp2(tabrev, 'SIZZ')
                call tbexp2(tabmdb, 'SIZZ')
                call tbexp2(tabrev, 'SIXY')
                call tbexp2(tabmdb, 'SIXY')
                call tbexp2(tabrev, 'COOR_X')
                call tbexp2(tabrev, 'COOR_Y')
            end if
        end if
! ======================================================================
! --- VERIFICATION DU PARAMETRE TEMP POUR LA TABLE DONNEE THERMIQUE ----
! ======================================================================
        call tbexp2(tabthr, 'TEMP')
        call jedetr(tabrev)
        call jedetr(tabmdb)
        call jedetr(tabthr)
    end do
! ======================================================================
    call jedema()
! ======================================================================
end subroutine
