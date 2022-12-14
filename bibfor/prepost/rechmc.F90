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
!
subroutine rechmc(ndim, temps, oridef, tabrev, tabmdb,&
                  norev, sigmrv, nomdb, sigmdb)
!
    implicit none
#include "jeveux.h"
#include "asterfort/detrsd.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/tbextb.h"
#include "asterfort/tbexve.h"
#include "asterfort/utmess.h"
#include "asterfort/wkvect.h"
    integer :: ndim, norev, nomdb
    real(kind=8) :: temps
    character(len=8) :: oridef, tabrev, tabmdb
    character(len=19) :: sigmrv, sigmdb
! --- BUT : CALCUL DES CHAMPS DES CONTRAINTES MECANIQUES ASSOCIES AU ---
! ------- : REVETEMENT ET AU METAL DE BASE -----------------------------
! ======================================================================
! IN  : NDIM   : DIMENSION DE L'ESPACE ---------------------------------
! --- : NBNO   : NOMBRE DE NOEUDS --------------------------------------
! --- : ORIDEF : ORIENTATION DU DEFAUT ---------------------------------
! --- : TABREV : TABLE ASSOCIEE AU REVETEMENT --------------------------
! --- : TABMDB : TABLE ASSOCIEE AU METAL DE BASE -----------------------
! OUT : NOREV  : NOMBRE DE NOEUDS COTE REVETEMENT ----------------------
! --- : SIGMRV : CHAMP DES CONTRAINTES ASSOCIE AU REVETEMENT -----------
! --- : NOMDB  : NOMBRE DE NOEUDS COTE METAL DE BASE -------------------
! --- : SIGMDB : CHAMP DES CONTRAINTES ASSOCIE AU METAL DE BASE --------
! ======================================================================
! ======================================================================
    integer :: jsigmr, jsigmb, ibid, ii, jcoorx, jcoory, iret
    integer :: jrevxx, jrevyy, jrevxy, jmdbxx, jmdbyy, jmdbxy
    real(kind=8) :: lprec, rt, cost, sint, sicot
    complex(kind=8) :: cbid
    character(len=8) :: lcrit, k8b
    character(len=19) :: tmprev, tmpmdb, coorxx, cooryy
    character(len=19) :: revxx, revyy, revxy, mdbxx, mdbyy, mdbxy
    character(len=24) :: valk(2)
! ======================================================================
    call jemarq()
! ======================================================================
! --- INITIALISATIONS --------------------------------------------------
! ======================================================================
    cbid=(0.d0,0.d0)
    ibid=0
    lcrit = 'RELATIF'
    lprec = 1.0d-06
    tmprev = '&&RECHMC.TMPREV'
    tmpmdb = '&&RECHMC.TMPMDB'
    revxx = '&&RECHMC.REVXX'
    revyy = '&&RECHMC.REVYY'
    revxy = '&&RECHMC.REVXY'
    mdbxx = '&&RECHMC.MDBXX'
    mdbyy = '&&RECHMC.MDBYY'
    mdbxy = '&&RECHMC.MDBXY'
    coorxx = '&&RECHMC.COORXX'
    cooryy = '&&RECHMC.COORYY'
! ======================================================================
! --- RECUPERATION DES SOUS-TABLES ASSOCIEES A L'INSTANT COURANT -------
! ======================================================================
    call tbextb(tabrev, 'V', tmprev, 1, 'INST',&
                'EQ', [ibid], [temps], [cbid], k8b,&
                [lprec], lcrit, iret)
    if (iret .eq. 10) then
        valk(1) = 'INST'
        valk(2) = tabrev
        call utmess('F', 'UTILITAI7_1', nk=2, valk=valk)
    else if (iret .eq. 20) then
        valk(1) = tabrev
        valk(2) = 'INST'
        call utmess('F', 'UTILITAI7_3', nk=2, valk=valk)
    endif
    call tbextb(tabmdb, 'V', tmpmdb, 1, 'INST',&
                'EQ', [ibid], [temps], [cbid], k8b,&
                [lprec], lcrit, iret)
    if (iret .eq. 10) then
        valk(1) = 'INST'
        valk(2) = tabmdb
        call utmess('F', 'UTILITAI7_1', nk=2, valk=valk)
    else if (iret .eq. 20) then
        valk(1) = tabmdb
        valk(2) = 'INST'
        call utmess('F', 'UTILITAI7_3', nk=2, valk=valk)
    endif
! ======================================================================
! --- PROBLEME EN DIMENSION 2 ------------------------------------------
! ======================================================================
    if (ndim .eq. 2) then
        if (oridef .eq. 'CIRC') then
! ======================================================================
! --- RECUPERATION DE LA LISTE DE CONTRAINTE SIYY COTE REVETEMENT ------
! ======================================================================
            call tbexve(tmprev, 'SIYY', sigmrv, 'V', norev,&
                        k8b)
! ======================================================================
! --- RECUPERATION DE LA LISTE DE CONTRAINTE SIYY COTE METAL DE BASE ---
! ======================================================================
            call tbexve(tmpmdb, 'SIYY', sigmdb, 'V', nomdb,&
                        k8b)
        else
! ======================================================================
! --- RECUPERATION DE LA LISTE DE CONTRAINTE SIZZ COTE REVETEMENT ------
! ======================================================================
            call tbexve(tmprev, 'SIZZ', sigmrv, 'V', norev,&
                        k8b)
! ======================================================================
! --- RECUPERATION DE LA LISTE DE CONTRAINTE SIZZ COTE METAL DE BASE ---
! ======================================================================
            call tbexve(tmpmdb, 'SIZZ', sigmdb, 'V', nomdb,&
                        k8b)
        endif
! ======================================================================
! --- PROBLEME EN DIMENSION 3 ------------------------------------------
! ======================================================================
    else
        if (oridef .eq. 'CIRC') then
! ======================================================================
! --- RECUPERATION DE LA LISTE DE CONTRAINTE SIZZ COTE REVETEMENT ------
! ======================================================================
            call tbexve(tmprev, 'SIZZ', sigmrv, 'V', norev,&
                        k8b)
! ======================================================================
! --- RECUPERATION DE LA LISTE DE CONTRAINTE SIZZ COTE METAL DE BASE ---
! ======================================================================
            call tbexve(tmpmdb, 'SIZZ', sigmdb, 'V', nomdb,&
                        k8b)
        else
! ======================================================================
! --- PASSAGE DE LA BASE CARTESIENNE (MODELE 3D) A LA BASE -------------
! --- CYLINDRIQUE ------------------------------------------------------
! ======================================================================
            call wkvect(sigmrv, 'V V R8', norev, jsigmr)
            call wkvect(sigmdb, 'V V R8', nomdb, jsigmb)
            call tbexve(tmprev, 'COOR_X', coorxx, 'V', norev,&
                        k8b)
            call tbexve(tmprev, 'COOR_Y', cooryy, 'V', norev,&
                        k8b)
            call jeveuo(coorxx, 'L', jcoorx)
            call jeveuo(cooryy, 'L', jcoory)
            rt = zr(jcoorx)*zr(jcoorx) + zr(jcoory)*zr(jcoory)
            cost = zr(jcoorx)*zr(jcoorx) / rt
            sint = zr(jcoory)*zr(jcoory) / rt
            sicot = 2 * zr(jcoorx)*zr(jcoory) / rt
! ======================================================================
! --- RECUPERATION DES LISTES DE CONTRAINTE SIXX - SIYY - SIXY ---------
! --- COTE REVETEMENT --------------------------------------------------
! ======================================================================
            call tbexve(tmprev, 'SIXX', revxx, 'V', norev,&
                        k8b)
            call tbexve(tmprev, 'SIYY', revyy, 'V', norev,&
                        k8b)
            call tbexve(tmprev, 'SIXY', revxy, 'V', norev,&
                        k8b)
            call jeveuo(revxx, 'L', jrevxx)
            call jeveuo(revyy, 'L', jrevyy)
            call jeveuo(revxy, 'L', jrevxy)
            do ii = 1, norev
                zr(jsigmr-1+ii) = sint * zr(jrevxx-1+ii) + cost * zr( jrevyy-1+ii) - sicot* zr(jr&
                                  &evxy-1+ii)
            end do
! ======================================================================
! --- RECUPERATION DES LISTES DE CONTRAINTE SIXX - SIYY - SIXY ---------
! --- COTE METAL DE BASE -----------------------------------------------
! ======================================================================
            call tbexve(tmpmdb, 'SIXX', mdbxx, 'V', nomdb,&
                        k8b)
            call tbexve(tmpmdb, 'SIYY', mdbyy, 'V', nomdb,&
                        k8b)
            call tbexve(tmpmdb, 'SIXY', mdbxy, 'V', nomdb,&
                        k8b)
            call jeveuo(mdbxx, 'L', jmdbxx)
            call jeveuo(mdbyy, 'L', jmdbyy)
            call jeveuo(mdbxy, 'L', jmdbxy)
            do ii = 1, nomdb
                zr(jsigmb-1+ii) = sint * zr(jmdbxx-1+ii) + cost * zr( jmdbyy-1+ii) - sicot* zr(jm&
                                  &dbxy-1+ii)
            end do
        endif
    endif
! ======================================================================
! --- DESTRUCTION DES TABLES INUTILES ----------------------------------
! ======================================================================
    call detrsd('TABLE', tmprev)
    call detrsd('TABLE', tmpmdb)
    call jedetr(revxx)
    call jedetr(revxy)
    call jedetr(revyy)
    call jedetr(mdbxx)
    call jedetr(mdbxy)
    call jedetr(mdbyy)
    call jedetr(coorxx)
    call jedetr(cooryy)
! ======================================================================
    call jedema()
! ======================================================================
end subroutine
