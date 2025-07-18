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
subroutine fonno62(resu, noma, ndim, iseg, noe, &
                   indr, nbnoel, vnor, vdir, basseg, &
                   vect, sens)
! aslint: disable=W1306
    implicit none
#include "jeveux.h"
#include "asterc/r8rddg.h"
#include "asterfort/assert.h"
#include "asterfort/dismoi.h"
#include "asterfort/getvtx.h"
#include "asterfort/jedema.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexnum.h"
#include "asterfort/trigom.h"
#include "asterfort/utmess.h"
#include "asterfort/char8_to_int.h"
#include "blas/ddot.h"
!
    character(len=8) :: resu, noma
    integer(kind=8) :: ndim, iseg, noe(4, 4)
    integer(kind=8) :: indr(2), nbnoel
    real(kind=8) :: vdir(2, 3), vnor(2, 3), vect(3), sens
    character(len=19) :: basseg
!
!     BUTS :
!        - VERIFIER LA COHERENCE DES 2 VECTEURS DIRECTION
!        - DANS LE CAS SYME : DETERMINER LA BONNE DIRECTION
!
!     ----------------------------------------------------
!
!  ENTREES
!       RESU   : NOM DU CONCEPT RESULTAT
!       NOMA   : NOM DU MAILLAGE
!       NDIM   : DIMENSION DU MODELE
!       ISEG   : INDICE DU SEGMENT DU FOND DE FISSURE COURANT
!       NOE    : NOEUDS DES FACES CONTENANT NA et NB ET APPARTENANT AUX
!                MAILLES CONNECTEES AU NOEUD SOMMET COURANT
!                ET AUX LEVRES
!       INDR   : INDICES DES FACES LIBRES DANS LA LISTE DES FACES
!                DES MAILLES CONNECTEES AU FOND DE FISSURE
!       NBNOEL : NOMBRE DE NOEUDS SOMMETS PAR ELEMENTS
!       VNOR   : VECTEUR NORMAL A LA SURFACE DE LA FISSURE
!       VDIR   : VECTEUR DANS LA DIRECTION DE PROPAGATION
!  ENTREE/SORTIE
!       BAFEFO : BASES LOCALES PAR SEGMENT DU FOND (NORMEE)
!
!       ----------------------------------------------------
!
    integer(kind=8) :: iamase, ityp, iatyma, jbasse
    integer(kind=8) :: i, j, iret, inp, compt, ino, ifl
    integer(kind=8) :: ilev, inor
    integer(kind=8) :: nblev, nn
    real(kind=8) :: s, ndir, nnor, alpha, angmax, beta
    real(kind=8) :: vecdir(ndim), vecnor(ndim), vnprec(ndim)
    real(kind=8) :: p, prvd1, prvd2
    character(len=6) :: syme
    character(len=8) :: type
    character(len=8), pointer :: mail(:) => null()
    blas_int :: b_incx, b_incy, b_n
    parameter(angmax=2.5d0)
!
!     -----------------------------------------------------------------
!
    call jemarq()
!
    call getvtx(' ', 'SYME', scal=syme, nbret=iret)
!
!     INDICE DE LA LEVRE A CONSIDERER
    ifl = 0
!
!     VERIFICATION DE LA PRESENCE DE LEVRE_SUP
    call jeexin(resu//'.LEVRESUP.MAIL', ilev)
!
!     VERIFICATION DE LA PRESENCE DE NORMALE
    call jeexin(resu//'.NORMALE', inor)
!
!     RECUPERATION DE L'ADRESSE DE LA BASE PAR SEGMENT DU FOND
    call jeveuo(basseg, 'E', jbasse)
!
!
!
!
!     1) VERIFICATION DE LA COHERENCE DES 2 VECTEURS DIRECTION
!        A FAIRE UNIQUEMENT SI LES LEVRES SONT COLLEES
!     --------------------------------------------------------
!
    if (inor .eq. 0) then
!        ALPHA = ANGLE ENTRE LES 2 VECTEURS (EN DEGRES)
        s = vdir(1, 1)*vdir(2, 1)+vdir(1, 2)*vdir(2, 2)+vdir(1, 3)*vdir(2, 3)
!
!        ATTENTION, NE JAMAIS UTILISER LA FONCTION FORTRAN ACOS
        alpha = trigom('ACOS', s)*r8rddg()
!
!        CAS SYMETRIQUE
        if (syme .eq. 'OUI') then
!
!          ANGLE DOIT ETRE EGAL A 180+-2,5 DEGRES, SINON CA VEUT DIRE
!          QUE L'HYPOTHESE DE LEVRES COLLEES EST FAUSSE : ON PLANTE
            if (abs(alpha-180.d0) .gt. angmax) then
                call utmess('F', 'RUPTURE0_34')
            end if
!
        else if (syme .eq. 'NON') then
!
!          ANGLE DOIT ETRE EGAL A 0+-5 DEGRES, SINON CA VEUT DIRE
!          QUE L'HYPOTHESE DE LEVRES COLLEES EST FAUSSE : ON PLANTE
            if (abs(alpha) .gt. 2.d0*angmax) then
                call utmess('F', 'RUPTURE0_34')
            end if
!
        end if
    end if
!
!
!     2) DANS LE CAS SYMETRIQUE, RECHERCHE DU BON VECTEUR DIRECTION
!     CETTE OPERATION N'EST FAITE QU'UNE SEULE FOIS POUR ISEG=1
!     -------------------------------------------------------------
!
    if (syme .eq. 'OUI' .and. iseg .eq. 1) then
!
!       CAS OU LES LEVRES SONT DONNEES
        if (ilev .ne. 0) then
!
            call jeveuo(resu//'.LEVRESUP.MAIL', 'L', vk8=mail)
            call jelira(resu//'.LEVRESUP.MAIL', 'LONUTI', nblev)
!
            call jeveuo(noma//'.TYPMAIL', 'L', iatyma)
!
!         BOUCLE SUR LES MAILLES DES LEVRES POUR TROUVER LE BON COTE
            do i = 1, nblev
                iret = char8_to_int(mail(i))
                call jeveuo(jexnum(noma//'.CONNEX', iret), 'L', iamase)
                ityp = iatyma-1+iret
                call jenuno(jexnum('&CATA.TM.NOMTM', zi(ityp)), type)
                call dismoi('NBNO_TYPMAIL', type, 'TYPE_MAILLE', repi=nn)
                do inp = 1, 2
                    compt = 0
                    do j = 1, nn
                        do ino = 1, nbnoel
                            if (zi(iamase-1+j) .eq. noe(indr(inp), ino)) then
                                compt = compt+1
                            end if
                        end do
                    end do
!             ON A TROUVE UNE FACE COINCIDENTE A UNE LEVRE, ON SORT
                    if (compt .eq. nbnoel) then
                        ifl = inp
                        goto 300
                    end if
                end do
            end do
        else
!
!       SINON, PB CAR LA DEFINITION DES LEVRE EST OBLIGATOIRE EN SYNTAXE
!
            ASSERT(.FALSE.)
        end if
!
    end if
!
300 continue
!
!     3) CALCUL DES VRAIS VECTEURS DIRECTION ET NORMAL
!     ------------------------------------------------
!
!     CAS OU IL FAUT PRENDRE LA MOYENNE DES 2 VECTEURS
    if (syme .eq. 'NON') then
!
        ndir = sqrt( &
               (vdir(1, 1)+vdir(2, 1))**2+(vdir(1, 2)+vdir(2, 2))**2+(vdir(1, 3)+vdir(2, 3))**2)
!
        nnor = sqrt( &
               (vnor(1, 1)+vnor(2, 1))**2+(vnor(1, 2)+vnor(2, 2))**2+(vnor(1, 3)+vnor(2, 3))**2)
!
        do i = 1, ndim
            vecdir(i) = (vdir(1, i)+vdir(2, i))/ndir
            vecnor(i) = sens*(vnor(1, i)+vnor(2, i))/nnor
        end do
!
!       LE VECTEUR NORMAL DOIT ALLER DE LA LEVRE INF
!       VERS LA LEVRE SUP
        if ((iseg .eq. 1) .and. (ilev .ne. 0)) then
            b_n = to_blas_int(ndim)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            p = ddot(b_n, vecnor, b_incx, vect, b_incy)
            if (p .lt. 0.d0) then
                sens = -1.d0
                do i = 1, ndim
                    vecnor(i) = sens*vecnor(i)
                end do
            end if
        end if
!
!     CAS OU IL NE FAUT PRENDRE QU'UN SEUL VECTEUR
    else if (syme .eq. 'OUI') then
!
!       POUR LE 1ER VECTEUR, ON CONNAIT DEJA IFL
        if (iseg .eq. 1) then
!
            ifl = ifl
!
!       POUR LES AUTRES SEGMENTS, ON PROCEDE PAR CONTINUITE
!       DES VDIR AVEC LE VECTEUR PRECEDENT
        else if (iseg .gt. 1) then
!
            prvd1 = vdir(1, 1)*zr(jbasse-1+2*ndim*(iseg-1-1)+1+ndim)+vdir(1, 2)*zr(jbasse-1+2*ndi&
                    &m*(iseg-1-1)+2+ndim)+vdir(1, 3)*zr(jbasse-1+2*ndim*(iseg-1-1)+3+ndim)
!
            prvd2 = vdir(2, 1)*zr(jbasse-1+2*ndim*(iseg-1-1)+1+ndim)+vdir(2, 2)*zr(jbasse-1+2*ndi&
                    &m*(iseg-1-1)+2+ndim)+vdir(2, 3)*zr(jbasse-1+2*ndim*(iseg-1-1)+3+ndim)
!
            if (prvd1 .gt. 0) then
!           LE VECTEUR VDIR PRECEDENT EST EN CONFORMITE AVEC IFL=1
                ifl = 1
            else if (prvd2 .gt. 0) then
!           LE VECTEUR VDIR PRECEDENT EST EN CONFORMITE AVEC IFL=2
                ifl = 2
            else
                ASSERT(.false.)
            end if
!
        end if
!
        ASSERT(ifl .ne. 0)
!
        ndir = sqrt(vdir(ifl, 1)**2+vdir(ifl, 2)**2+vdir(ifl, 3)**2)
!
        nnor = sqrt(vnor(ifl, 1)**2+vnor(ifl, 2)**2+vnor(ifl, 3)**2)
!
        do i = 1, ndim
            vecdir(i) = vdir(ifl, i)/ndir
            vecnor(i) = sens*vnor(ifl, i)/nnor
        end do
!
!       LE VECTEUR NORMAL DOIT ALLER DE LA LEVRE INF
!       VERS LA LEVRE SUP
        if ((iseg .eq. 1) .and. (ilev .ne. 0)) then
            b_n = to_blas_int(ndim)
            b_incx = to_blas_int(1)
            b_incy = to_blas_int(1)
            p = ddot(b_n, vecnor, b_incx, vect, b_incy)
            if (p .lt. 0.d0) then
                sens = -1.d0
                do i = 1, ndim
                    vecnor(i) = sens*vecnor(i)
                end do
            end if
        end if
!
    end if
!
!     4) ECRITURE DE LA BASE PAR SEGMENT DU FOND
!     -------------------------------------------
!
    do i = 1, ndim
        zr(jbasse-1+2*ndim*(iseg-1)+i) = vecnor(i)
        zr(jbasse-1+2*ndim*(iseg-1)+i+ndim) = vecdir(i)
    end do
!
!
!     5) VERIF QUE LE VECTEUR NORMAL N'EST PAS TROP DIFFERENT DU
!     VECTEUR NORMAL DU SEGMENT PRECEDENT (ON TOLERE 10 DEGRES)
!     ---------------------------------------------------------
!
    if (iseg .gt. 1) then
!       RECUP DU VECTEUR NORMAL PRECEDENT
        do i = 1, ndim
            vnprec(i) = zr(jbasse-1+2*ndim*(iseg-2)+i)
        end do
        b_n = to_blas_int(ndim)
        b_incx = to_blas_int(1)
        b_incy = to_blas_int(1)
        s = ddot(b_n, vecnor, b_incx, vnprec, b_incy)
        beta = trigom('ACOS', s)*r8rddg()
        if (abs(beta) .gt. 10.d0) then
            call utmess('A', 'RUPTURE0_61')
        end if
    end if
!
    call jedema()
end subroutine
