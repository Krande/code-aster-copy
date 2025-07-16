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

subroutine elim75(nomres, resgen, matgen, massgen)
    use aster_petsc_module
    use elg_data_module
    implicit none
#include "asterf_types.h"
#include "asterf_petsc.h"
#include "jeveux.h"
#include "asterfort/assert.h"
#include "asterfort/dcapno.h"
#include "asterfort/dismoi.h"
#include "asterfort/elg_allocvr.h"
#include "asterfort/elg_calc_solu.h"
#include "asterfort/getvis.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jelibe.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/refdcp.h"
#include "asterfort/rsadpa.h"
#include "asterfort/rscrsd.h"
#include "asterfort/rsexch.h"
#include "asterfort/rsnoch.h"
#include "asterfort/rsorac.h"
#include "asterfort/titre.h"
#include "asterfort/vtcreb.h"
#include "asterfort/wkvect.h"
!
    character(len=8) :: nomres, resgen
    character(len=19) :: matgen, massgen
!-----------------------------------------------------------------------
!
!  BUT : < RESTITUTION GENERALISEE > D'UNE MATRICE REDUITE PAR
!  LA COMMANDE ELIM_LAGR
!
!-----------------------------------------------------------------------!
!  LE CONCEPT RESULTAT EST UN RESULTAT COMPOSE "MODE_MECA"
!               OU "MODE_MECA_C"
!-----------------------------------------------------------------------
!
! NOMRES /I/ : NOM K8 DU CONCEPT MODE MECA RESULTAT
! RESGEN /I/ : NOM K8 DU MODE_GENE AMONT
!
!
#ifdef ASTER_HAVE_PETSC
    integer(kind=8) :: i, j, ibid, ier, iord
    integer(kind=8) :: jbid, nbmod, neq, nsecm
    integer(kind=8) :: nno, iadpar(13), iadpas(13), tmod(1)
    real(kind=8) :: rbid, omega2
    complex(kind=8) :: cbid
    character(len=1) :: typsca
    character(len=8) :: kbid
    character(len=14) :: numddl
    character(len=16) :: depl, nompar(13)
    character(len=19) :: chamne, raid, mass
    character(len=24) :: chamol
    aster_logical :: zcmplx
    complex(kind=8), pointer :: vgene_c(:) => null(), vnew_c(:) => null()
    real(kind=8), pointer :: vgene_r(:) => null(), vnew_r(:) => null()
    real(kind=8), pointer :: deplr(:) => null(), deplc(:) => null()
    character(len=24), pointer :: refa(:) => null()
!
    PetscErrorCode  :: ierr
    PetscInt :: n1, n2, ke_mass
    !
!-----------------------------------------------------------------------
    data depl/'DEPL            '/
    data nompar/'FREQ', 'RIGI_GENE', 'MASS_GENE', 'OMEGA2',&
    &           'AMOR_REDUIT',&
    &           'FACT_PARTICI_DX', 'FACT_PARTICI_DY', 'FACT_PARTICI_DZ',&
    &           'MASS_EFFE_DX', 'MASS_EFFE_DY', 'MASS_EFFE_DZ',&
    &           'NUME_MODE',&
    &           'TYPE_MODE'/
!-----------------------------------------------------------------------
    call jemarq()
!
!-----RECUPERATION NOMBRE DE MODES PROPRES CALCULES---------------------
!
    zcmplx = .false.
!
    call dcapno(resgen, depl, 1_8, chamol)
    call jelira(chamol, 'TYPE', cval=typsca)
    if (typsca .eq. 'C') zcmplx = .true.
!
    call rsorac(resgen, 'LONUTI', 0_8, rbid, kbid, &
                cbid, rbid, kbid, tmod, 1_8, &
                ibid)
    nbmod = tmod(1)
!
! --- ON RESTITUE SUR TOUS LES MODES OU SUR QUELQUES MODES:
!
    call getvis(' ', 'NUME_ORDRE', nbval=0_8, nbret=nno)
    if (nno .ne. 0) then
        nbmod = -nno
        call wkvect('&&ELIM75.NUME', 'V V I', nbmod, jbid)
        call getvis(' ', 'NUME_ORDRE', nbval=nbmod, vect=zi(jbid), nbret=nno)
    else
        call wkvect('&&ELIM75.NUME', 'V V I', nbmod, jbid)
        do i = 1, nbmod
            zi(jbid+i-1) = i
        end do
    end if
!
! --- ALLOCATION STRUCTURE DE DONNEES RESULTAT
!
    if (zcmplx) then
        call rscrsd('G', nomres, 'MODE_MECA_C', nbmod)
    else
        call rscrsd('G', nomres, 'MODE_MECA', nbmod)
    end if

!
! ---- RESTITUTION PROPREMENT DITE
!
!  on recupère le nom de la matrice de rigidité d'origine
    call jeveuo(matgen//'.REFA', 'L', vk24=refa)
    raid = refa(20) (1:19)
!  on recupère le nom de la matrice de masse d'origine
    call jeveuo(massgen//'.REFA', 'L', vk24=refa)
    mass = refa(20) (1:19)
!  on récupère le numéro d'instance ke_mass associé à la matrice de masse
    call elg_gest_data('CHERCHE', mass, massgen, ' ')
    ke_mass = ke
!  on récupère le numéro d'instance ke associé à la matrice de rigidité
    call elg_gest_data('CHERCHE', raid, matgen, ' ')
!  on construit les vecteurs nécessaires pour la restitution
    call MatGetSize(elg_context(ke)%matc, n2, n1, ierr)
    ASSERT(ierr == 0)
    call elg_allocvr(elg_context(ke)%vecb, int(n1, 8))
    call elg_allocvr(elg_context(ke)%vecc, int(n2, 8))
!  on les remplie de 0 car c'est une résolution modale
    call VecZeroEntries(elg_context(ke)%vecb, ierr)
    ASSERT(ierr == 0)
    call VecDuplicate(elg_context(ke)%vecb, elg_context(ke)%vx0, ierr)
    ASSERT(ierr == 0)
    call VecZeroEntries(elg_context(ke)%vecc, ierr)
    ASSERT(ierr == 0)

    call dismoi('NB_EQUA', raid, 'MATR_ASSE', repi=neq)
    call dismoi('NOM_NUME_DDL', raid, 'MATR_ASSE', repk=numddl)
    call wkvect('&&ELIM75.DEPL_R', 'V V R', neq, vr=deplr)
    call wkvect('&&ELIM75.DEPL_C', 'V V R', neq, vr=deplc)
!
!
! ------ BOUCLE SUR LES MODES A RESTITUER
!
    do i = 1, nbmod
        iord = zi(jbid+i-1)
!
! --------- REQUETE NOM ET ADRESSE CHAMNO GENERALISE
!
        call dcapno(resgen, depl, iord, chamol)
        if (zcmplx) then
            call jeveuo(chamol, 'L', vc=vgene_c)
        else
            call jeveuo(chamol, 'L', vr=vgene_r)
        end if
!
! --------- REQUETE NOM ET ADRESSE NOUVEAU CHAMNO
!
        call rsexch(' ', nomres, depl, i, chamne, ier)
        if (zcmplx) then
            call vtcreb(chamne, 'G', 'C', nume_ddlz=numddl)
            call jeveuo(chamne//'.VALE', 'E', vc=vnew_c)
        else
            call vtcreb(chamne, 'G', 'R', nume_ddlz=numddl)
            call jeveuo(chamne//'.VALE', 'E', vr=vnew_r)
        end if
!
        call rsadpa(resgen, 'L', 13_8, nompar, iord, 0_8, tjv=iadpar, styp=kbid, istop=0_8)
!
        omega2 = zr(iadpar(4))
!
        nsecm = 1
        if (zcmplx) then
            ! dans le cas complexe, on restitue les parties réelle et imaginaire
            ! séparemment et on les somme
            call elg_calc_solu(raid, nsecm, real(vgene_c), deplr, omega2, ke_mass)
            call elg_calc_solu(raid, nsecm, aimag(vgene_c), deplc, omega2, ke_mass)
            vnew_c = dcmplx(deplr, deplc)
        else
            call elg_calc_solu(raid, nsecm, vgene_r, vnew_r)
        end if
!
        call rsnoch(nomres, depl, i)
        call rsadpa(nomres, 'E', 13_8, nompar, i, 0_8, tjv=iadpas, styp=kbid)
        do j = 1, 11
            zr(iadpas(j)) = zr(iadpar(j))
        end do
        zi(iadpas(12)) = zi(iadpar(12))
        zk16(iadpas(13)) = 'MODE_DYN'
!
        call jelibe(chamol)
    end do

! --- ENREGISTREMENT DE LA REFERENCE DYNAMIQUE
    call refdcp(resgen, nomres)
!
! --- MENAGE
    call jedetr('&&ELIM75.NUME')
    call jedetr('&&ELIM75.DEPL_R')
    call jedetr('&&ELIM75.DEPL_C')
!
    call titre()
    call jedema()

#else
    character(len=19) :: kbid
    kbid = matgen
#endif
end subroutine
