! --------------------------------------------------------------------
! Copyright (C) 1991 - 2025 - EDF R&D - www.code-aster.org
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
subroutine fointa(ipif, nbpu, nompu, valpu, resu)
    implicit none
#include "jeveux.h"
#include "asterc/r8prem.h"
#include "asterc/r8vide.h"
#include "asterfort/assert.h"
#include "asterfort/fiintf.h"
#include "asterfort/focoli.h"
#include "asterfort/fointn.h"
#include "asterfort/folocx.h"
#include "asterfort/tecael.h"
#include "asterfort/utmess.h"
!
    integer(kind=8) :: ipif, nbpu
    real(kind=8) :: valpu(*), resu
    character(len=*) :: nompu(*)
!     INTERPOLATION POUR CALCULER RESU = F(X,Y,Z,...)
! ----------------------------------------------------------------------
! IN  IPIF  : POINTEUR DANS LE MATERIAU CODE (FONCTION OU NAPPE)
! IN  NBPU  : NOMBRE DE PARAMETRES DANS NOMPU ET VALPU
! IN  NOMPU : NOMS DES PARAMETRES "UTILISATEUR"
! IN  VALPU : VALEURS DES PARAMETRES "UTILISATEUR"
! OUT RESU  : R : RESULTAT DE L'INTERPOLATION
! ----------------------------------------------------------------------
!
!
!
    integer(kind=8) :: npar(2), indfct, jpro, jpar, lpara, nbvn, nbpara, i
    integer(kind=8) :: nupar, nbpt, jval, inume, ier, iret
    real(kind=8) :: tab(4), rpar, rvar, epsi, tresu(1)
    character(len=1) :: coli
    character(len=19) :: nomf
    character(len=24) :: nompf(2)
!     ------------------------------------------------------------------
    integer(kind=8) :: iadzi, iazk24
    character(len=24) :: valk(3)
! ----------------------------------------------------------------------
! PARAMETER ASSOCIE AU MATERIAU CODE
!
    parameter(indfct=7)
! ----------------------------------------------------------------------
!     FONCTION EN LIGNE
!
#define linlin(x,x1,y1,x2,y2) y1+(x-x1)*(y2-y1)/(x2-x1)
#define linlog(x,x1,y1,x2,y2) exp(log(y1)+(x-x1)*(log(y2)-log(y1)) \
    /(x2-x1))
#define loglog(x,x1,y1,x2,y2) exp(log(y1)+(log(x)-log(x1))*(log(y2) \
    -log(y1))/(log(x2)-log(x1)))
#define loglin(x,x1,y1,x2,y2) y1+(log(x)-log(x1))*(y2-y1) \
    /(log(x2)-log(x1))
!     ------------------------------------------------------------------
    npar(1) = 0
    npar(2) = 0
    nomf = ' '
    jpro = zi(ipif+1)
    jpar = zi(ipif+2)
    resu = r8vide()
    ier = 0
!
    epsi = sqrt(r8prem())
!
! --- FONCTION "CONSTANT"
!
    if (zk24(jpro) .eq. 'CONSTANT') then
!         ------------------------
        resu = zr(jpar+1)
        goto 999
!
!
! --- FONCTION "INTERPRE" : FORMULE
!
    else if (zk24(jpro) .eq. 'INTERPRE') then
!             ------------------------
        nomf = zk24(jpro+5) (1:19)
        call fiintf(nomf, nbpu, nompu, valpu, iret, &
                    'E', tresu)
        resu = tresu(1)
        if (iret .ne. 0) then
            call tecael(iadzi, iazk24)
            call utmess('F+', 'FONCT0_9', sk=nomf)
            call utmess('F', 'FONCT0_10', sk=zk24(iazk24-1+3))
        end if
        goto 999
!
! --- AUTRES TYPES DE FONCTION
!
    else if (zk24(jpro) .eq. 'FONCTION') then
        nbpara = 1
        nompf(1) = zk24(jpro+2)
        nomf = zk24(jpro+5) (1:19)
!
    else if (zk24(jpro) .eq. 'NAPPE') then
        nbpara = 2
        nompf(1) = zk24(jpro+2)
        nompf(2) = zk24(jpro+6)
        nomf = zk24(jpro+5) (1:19)
!
    else
        ASSERT(ASTER_FALSE)
    end if
!
    do i = 1, nbpara
        do nupar = 1, nbpu
            if (nompu(nupar) .eq. nompf(i)) then
!           -- SI UN PARAMETRE EST FOURNI PLUSIEURS FOIS
!              ON PREND LE DERNIER (VOIR RCVALB)
                npar(i) = nupar
            end if
        end do
        if (npar(i) .eq. 0) then
            valk(1) = nomf
            valk(2) = nompf(i)
            call tecael(iadzi, iazk24)
            valk(3) = zk24(iazk24-1+3)
            if (valk(2) .eq. 'X' .or. valk(2) .eq. 'Y' .or. valk(2) .eq. 'Z') then
                call utmess('F', 'CALCULEL6_65', sk=nomf)
            else
                call utmess('F', 'CALCULEL6_62', nk=3, valk=valk)
            end if
        end if
    end do
!
! =====================================================================
!                          F O N C T I O N
! =====================================================================
!
    if (zk24(jpro) .eq. 'FONCTION') then
        nbpt = zi(ipif)
        jval = jpar+nbpt
        rvar = valpu(npar(1))
        call folocx(zr(jpar), nbpt, rvar, zk24(jpro+4), zi(ipif+indfct), &
                    epsi, coli, ier)
        if (ier .ne. 0) goto 999
        call focoli(zi(ipif+indfct), coli, zk24(jpro+1), zr(jpar), zr(jval), &
                    rvar, resu, ier)
        if (ier .ne. 0) goto 999
!
! =====================================================================
!                            N A P P E
! =====================================================================
!
    else if (zk24(jpro) .eq. 'NAPPE   ') then
        rpar = valpu(npar(1))
        rvar = valpu(npar(2))
        lpara = zi(ipif+4)
        nbvn = zi(ipif+5)
        call folocx(zr(lpara), nbvn, rpar, zk24(jpro+4), zi(ipif+indfct), &
                    epsi, coli, ier)
        if (ier .ne. 0) goto 999
        inume = zi(ipif+indfct)
!
        if (coli .eq. 'C') then
            call fointn(ipif, nomf, rvar, inume, epsi, &
                        resu, ier)
            if (ier .ne. 0) goto 999
!
        else if (coli .eq. 'I') then
            call fointn(ipif, nomf, rvar, inume, epsi, &
                        tab(3), ier)
            if (ier .ne. 0) goto 999
            call fointn(ipif, nomf, rvar, inume+1, epsi, &
                        tab(4), ier)
            if (ier .ne. 0) goto 999
!
! ------- INTERPOLATION FINALE SUR LES PARAMETRES
!
            tab(1) = zr(lpara+inume-1)
            tab(2) = zr(lpara+inume)
            if (zk24(jpro+1) .eq. 'LIN LIN ') then
                resu = linlin(rpar, tab(1), tab(3), tab(2), tab(4))
            else if (zk24(jpro+1) .eq. 'LIN LOG ') then
                resu = linlog(rpar, tab(1), tab(3), tab(2), tab(4))
            else if (zk24(jpro+1) .eq. 'LOG LOG ') then
                resu = loglog(rpar, tab(1), tab(3), tab(2), tab(4))
            else if (zk24(jpro+1) .eq. 'LOG LIN ') then
                resu = loglin(rpar, tab(1), tab(3), tab(2), tab(4))
            end if
!
        else if (coli .eq. 'E') then
            call fointn(ipif, nomf, rvar, inume, epsi, &
                        tab(3), ier)
            if (ier .ne. 0) goto 999
            call fointn(ipif, nomf, rvar, inume+1, epsi, &
                        tab(4), ier)
            if (ier .ne. 0) goto 999
            tab(1) = zr(lpara+inume-1)
            tab(2) = zr(lpara+inume)
            resu = linlin(rpar, tab(1), tab(3), tab(2), tab(4))
!
        else
            call utmess('F', 'UTILITAI2_13', sk=coli)
        end if
!
    else
        call utmess('F', 'UTILITAI2_14', sk=zk24(jpro))
    end if
!
999 continue
    if (ier .ne. 0) then
        call utmess('F', 'CALCULEL6_63', sk=nomf, si=ier)
    end if
!
!
end subroutine
