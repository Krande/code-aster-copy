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
subroutine assvec(jvBase, vectAsseZ,&
                  nbVectElem, listVectElem, coefVectElem,&
                  numeDofZ_, vectAsseForNumeZ_,&
                  vectScalType_)
!
implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/indik8.h"
#include "asterfort/asmpi_barrier.h"
#include "asterfort/asmpi_comm_jev.h"
#include "asterfort/assert.h"
#include "asterfort/corddl.h"
#include "asterfort/crelil.h"
#include "asterfort/dbgobj.h"
#include "asterfort/detrsd.h"
#include "asterfort/digdel.h"
#include "asterfort/dismoi.h"
#include "asterfort/infniv.h"
#include "asterfort/jaexin.h"
#include "asterfort/jecreo.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jeecra.h"
#include "asterfort/jeexin.h"
#include "asterfort/jelibe.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jenonu.h"
#include "asterfort/jenuno.h"
#include "asterfort/jeveuo.h"
#include "asterfort/jexatr.h"
#include "asterfort/jexnom.h"
#include "asterfort/jexnum.h"
#include "asterfort/nbec.h"
#include "asterfort/nbno.h"
#include "asterfort/utmess.h"
#include "asterfort/uttcpr.h"
#include "asterfort/uttcpu.h"
#include "asterfort/vtcreb.h"
#include "asterfort/wkvect.h"
#include "asterfort/asseVectSuper.h"
#include "asterfort/asseVectField.h"
#include "asterfort/getDistributionParameters.h"
!
character(len=1), intent(in) :: jvBase
character(len=*), intent(in) :: vectAsseZ
integer, intent(in) :: nbVectElem
character(len=*), intent(in) :: listVectElem(nbVectElem)
real(kind=8), intent(in) :: coefVectElem(nbVectElem)
character(len=*), optional, intent(in) :: vectAsseForNumeZ_, numeDofZ_
integer, optional, intent(in) :: vectScalType_
!
! --------------------------------------------------------------------------------------------------
!
! OUT K19 VEC   : NOM DU CHAM_NO RESULTAT
!                CHAM_NO ::= CHAM_NO_GD + OBJETS PROVISOIRES POUR L'ASS.
! IN  K* BASE   : NOM DE LA BASE SUR LAQUELLE ON VEUT CREER LE CHAM_NO
! IN  I  NBVEC  : NOMBRE DE VECT_ELEM A ASSEMBLER DANS VEC
! IN  K* TLIVEC : LISTE DES VECT_ELEM A ASSEMBLER
! IN  R  LICOEF : LISTE DES COEF. MULTIPLICATEURS DES VECT_ELEM
! IN  K14 NU    : NOM D'UN NUME_DDL (LE STOCKAGE N'EST PAS NECESSAIRE)
!
! IN  K* VECPRO: NOM D'UN CHAM_NO MODELE(NU OU VECPRO EST OBLIGATOIRE)
! IN  I  TYPE   : TYPE DU VECTEUR ASSEMBLE : 1 --> REEL
!                                            2 --> COMPLEXE
!
! --------------------------------------------------------------------------------------------------
! - Convention: first LIGREL (model) is on mesh
    integer, parameter :: ligrelMeshIndx = 1
    character(len=24), parameter :: ligrelMesh = '&MAILLA'
    aster_logical, parameter :: dbg = ASTER_FALSE
    integer :: physQuan, nec, nlili
    integer, parameter :: nbecmx =10
    character(len=8) :: mesh, model, vectElemModel, nogdsi, nogdco
    character(len=14) :: numeDof, answer
    character(len=24) :: vectRefeJv, vectDescJv, vectValeJv
    character(len=24), pointer :: vectRefe(:) => null()
    integer, pointer :: vectDesc(:) => null()
    character(len=19) :: vectAsse, vectAsseForNume
    character(len=19) :: vectElem, resuElem, numeEqua
    character(len=24) :: numePrnoJv, numeNueqJv, numeNequJv
    character(len=24) :: numeLiliJv, vectAsseLili, ligrelName
    aster_logical :: ldist, ldgrel, compSuperElement,lparallel_mesh
    integer :: iDofMode, iVectElem
    integer :: iancmp, ianueq, iapsdl, iad1
    integer :: icmp, iconx2
    integer :: idprn1, idprn2, jresl, iElem
    integer :: iGrel, iDof, ilim
    integer :: liliNume, ligrelNume, jvVectElem
    integer :: iResuElem, iret, jec, jvale, iNodeMode
    integer :: lgncmp, mode, nbNode, meshNbCell
    integer :: nbResuElem, nbSuperElement, nbCmp, nbCmpMode, nbDofMode, nbElem, nbEqua, nbDof
    integer :: meshNbNode, nmxcmp, nbNodeMode, nugd, elemNume, iexi
    integer :: icodla(nbecmx), icodge(nbecmx), lshift
    integer :: admodl, lcmodl, ifm, niv, rang, nbproc
    real(kind=8) :: temps(7)
    character(len=24), pointer :: refe(:) => null(), noli(:) => null()
    integer, pointer :: adli(:) => null(), adne(:) => null()
    character(len=24), pointer :: relr(:) => null()
    integer, pointer :: numsd(:) => null()
    integer, pointer :: desc(:) => null(), nequ(:) => null()
    integer, pointer :: connex(:) => null()
    character(len=8), pointer :: nomacr(:) => null()
    real(kind=8) :: vectElemCoef, elemCoef
    integer :: vectScalType
! --------------------------------------------------------------------------------------------------
#define zzngel(ligrelNume) adli(1+3*(ligrelNume-1))
#define zznelg(ligrelNume,iGrel) zi(adli(1+3*(ligrelNume-1)+2)+iGrel)-\
                                 zi(adli(1+3*(ligrelNume-1)+2)+iGrel-1)-1
#define zzliel(ligrelNume,iGrel,iElem) zi(adli(1+3*(ligrelNume-1)+1)-1+\
                                       zi(adli(1+3*(ligrelNume-1)+2)+iGrel-1)+iElem-1)
!
! --------------------------------------------------------------------------------------------------
!
    call jemarq()
    call infniv(ifm, niv)
    vectAsse = vectAsseZ
    numeDof = ' '
    if (present(numeDofZ_)) then
        numeDof = numeDofZ_
    endif
    vectAsseForNume =  ' '
    if (present(vectAsseForNumeZ_)) then
        vectAsseForNume = vectAsseForNumeZ_
    endif
    vectScalType = 1
    if (present(vectScalType_)) then
        vectScalType = vectScalType_
    endif

! - Acces to description of local mode
    call jeveuo(jexatr('&CATA.TE.MODELOC', 'LONCUM'), 'L', lcmodl)
    call jeveuo(jexnum('&CATA.TE.MODELOC', 1), 'L', admodl)

! - Prepare list of elementary vectors
    call detrsd('CHAM_NO', vectAsse)
    call wkvect(vectAsse//'.LIVE', jvBase//' V K24 ', nbVectElem, jvVectElem)
    do iVectElem = 1, nbVectElem
        zk24(jvVectElem-1+iVectElem) = listVectElem(iVectElem)
    end do

! - Get numbering
    if (numeDof(1:1) .eq. ' ') then
        call jeveuo(vectAsseForNume//'.REFE', 'L', vk24=refe)
        numeDof = refe(2)(1:14)
    endif
    numeEqua = numeDof(1:14)//'.NUME'

! - Get model
    call dismoi('NOM_MODELE', numeDof, 'NUME_DDL', repk=model)

! - Get mesh
    call dismoi('NOM_MAILLA', numeDof, 'NUME_DDL', repk=mesh)
    call dismoi('PARALLEL_MESH', mesh, 'MAILLAGE', repk=answer)
    lparallel_mesh=(answer.eq.'OUI')
    if (.not.lparallel_mesh) then
        call asmpi_barrier()
    endif
    call jeexin(mesh(1:8)//'.CONNEX', iret)
    if (iret .gt. 0) then
        call jeveuo(mesh(1:8)//'.CONNEX', 'L', vi=connex)
        call jeveuo(jexatr(mesh(1:8)//'.CONNEX', 'LONCUM'), 'L', iconx2)
    endif
    call dismoi('NB_NO_MAILLA', model, 'MODELE', repi=meshNbNode)

! - Start timers for calcul.F90 monitoring
    call uttcpu('CPU.CALC.1', 'DEBUT', ' ')
    call uttcpu('CPU.ASSE.1', 'DEBUT', ' ')
    call uttcpu('CPU.ASSE.3', 'DEBUT', ' ')

! - Create list of ligrel (LILI) and objects ADNE and ADLI
    vectAsseLili = vectAsse//'.LILI'
    call crelil('C', nbVectElem, jvVectElem, vectAsseLili, 'V',&
                ligrelMesh, vectAsse, physQuan, mesh, nec,&
                nbCmp, ilim, nlili, meshNbCell)

! - No elementary terms to assemble, but maybe a nodal field
    if (nlili .eq. 1) then
        call vtcreb(vectAsse, jvBase, 'R',&
                    nume_ddlz = numeDof,&
                    nb_equa_outz = nbEqua)
        nbDof  = nbEqua
        goto 270
    endif
    call jeveuo(vectAsse(1:19)//'.ADLI', 'E', vi=adli)
    call jeveuo(vectAsse(1:19)//'.ADNE', 'E', vi=adne)

! - Get parameters for distribution of elementary vectors
    call getDistributionParameters(nbVectElem, listVectElem,&
                                   ldist, ldgrel,&
                                   rang, nbproc,&
                                   numsd)

! - Get parameters about physical quantity
    call dismoi('NOM_GD', numeDof, 'NUME_DDL', repk=nogdco)
    call dismoi('NOM_GD_SI', nogdco, 'GRANDEUR', repk=nogdsi)
    call dismoi('NB_CMP_MAX', nogdsi, 'GRANDEUR', repi=nmxcmp)
    call dismoi('NUM_GD_SI', nogdsi, 'GRANDEUR', repi=nugd)
    nec = nbec(nugd)
    nbCmp = nmxcmp

! - Objects for DOF
    icodla = 0
    icodge = 0
    call wkvect('&&ASSVEC.POSDDL', 'V V I', nmxcmp, iapsdl)

! - Preparation for super elements
    call jeexin(mesh//'.NOMACR', iret)
    if (iret .gt. 0) then
        call jeveuo(mesh//'.NOMACR', 'L', vk8=nomacr)
        call jeveuo(jexnom('&CATA.GD.NOMCMP', nogdsi), 'L', iancmp)
        call jelira(jexnom('&CATA.GD.NOMCMP', nogdsi), 'LONMAX', lgncmp)
        icmp=indik8(zk8(iancmp),'LAGR',1,lgncmp)
        ASSERT(icmp.ne.0)
        ASSERT(icmp.le.30)
!       -- ICODLA EST L'ENTIER CODE CORRESPONDANT A LA CMP "LAGR"
        jec=(icmp-1)/30+1
        icodla(jec)=lshift(1,icmp)
    endif

! - Start local timers if required
    if (niv .ge. 2) then
        call uttcpu('CPU.ASSVEC', 'INIT ', ' ')
        call uttcpu('CPU.ASSVEC', 'DEBUT', ' ')
    endif

! - Acces to numbering objects
    numePrnoJv = numeEqua(1:19)//'.PRNO'
    numeLiliJv = numeEqua(1:19)//'.LILI'
    numeNueqJv = numeEqua(1:19)//'.NUEQ'
    numeNequJv = numeEqua(1:19)//'.NEQU'
    call jeveuo(numePrnoJv, 'L', idprn1)
    call jeveuo(jexatr(numePrnoJv, 'LONCUM'), 'L', idprn2)
    call jeveuo(numeNueqJv, 'L', ianueq)

! - Get number of equations
    call jeexin(numeNequJv, iexi)
    if (iexi.eq.0) then
        call jelira(numeNueqJv, 'LONMAX', nbEqua)
        nbDof  = nbEqua
    else
        call jeveuo(numeNequJv, 'L', vi = nequ)
        nbEqua = nequ(1)
        nbDof  = nequ(2)
    endif
    if (nbDof .eq. 0) then
        nbDof = nbEqua
    endif

! - Access to vector
    vectRefeJv = vectAsse//'.REFE'
    vectValeJv = vectAsse//'.VALE'
    vectDescJv = vectAsse//'.DESC'

! - Create vector objeccts
    call jecreo(vectRefeJv, jvBase//' V K24')
    call jeecra(vectRefeJv, 'LONMAX', 4)
    call jeecra(vectRefeJv, 'LONUTI', 4)
    call jeveuo(vectRefeJv, 'E', vk24 = vectRefe)
    call jecreo(vectDescJv, jvBase//' V I')
    call jeecra(vectDescJv, 'LONMAX', 2)
    call jeecra(vectDescJv, 'LONUTI', 2)
    call jeecra(vectDescJv, 'DOCU', cval='CHNO')
    call jeveuo(vectDescJv, 'E', vi = vectDesc)
    vectRefe(1) = mesh
    vectRefe(2) = numePrnoJv(1:14)//'.NUME'
    vectDesc(1) = physQuan
    vectDesc(2) = 1
    if (vectScalType .eq. 1) then
        call jecreo(vectValeJv, jvBase//' V R8')
    else if (vectScalType.eq.2) then
        call jecreo(vectValeJv, jvBase//' V C16')
    else
        ASSERT(ASTER_FALSE)
    endif
    call jeecra(vectValeJv, 'LONMAX', nbEqua)
    call jeecra(vectValeJv, 'LONUTI', nbEqua)
    call jeveuo(vectValeJv, 'E', jvale)

! - Loop on elementary vectors
    do iVectElem = 1, nbVectElem
        vectElemCoef = coefVectElem(iVectElem)
        vectElem = zk24(jvVectElem+iVectElem-1)(1:19)
        call dismoi('NOM_MODELE', vectElem, 'VECT_ELEM', repk=vectElemModel)
        if (vectElemModel .ne. model) then
            call utmess('F', 'ASSEMBLA_5')
        endif

! ----- Add super elements
        call dismoi('NB_SS_ACTI', vectElem, 'VECT_ELEM', repi=nbSuperElement)
        if (nbSuperElement.gt.0) then
            compSuperElement = ASTER_TRUE
            if (ldist .and. rang .ne. 0) compSuperElement = ASTER_FALSE
        else
            compSuperElement = ASTER_FALSE
        endif
        if (compSuperElement) then
            call asseVectSuper(model, mesh, vectElem,&
                               vectScalType, vectElemCoef,&
                               nomacr, meshNbNode,&
                               nec, nbecmx, nbCmp,&
                               icodla, icodge,&
                               idprn1, idprn2,&
                               iapsdl, ianueq, jvale, jresl)
        endif

! ----- Add standard finite elemnts
        call jeexin(vectElem//'.RELR', iret)
        if (iret .gt. 0) then
            call jelira(vectElem//'.RELR', 'LONUTI', nbResuElem)
            if (nbResuElem .gt. 0) call jeveuo(vectElem//'.RELR', 'L', vk24=relr)

! --------- Loop on elementary terms
            do iResuElem = 1, nbResuElem
                resuElem = relr(iResuElem)(1:19)
                call jeexin(resuElem//'.NOLI', iexi)
                if (iexi .eq. 0) cycle
                call jeveuo(resuElem//'.NOLI', 'L', vk24 = noli)
                ligrelName = noli(1)
                call jenonu(jexnom(vectAsseLili, ligrelName), ligrelNume)
                call jenonu(jexnom(numeLiliJv, ligrelName), liliNume)

! ------------- Loop on GREL
                do iGrel = 1, zzngel(ligrelNume)
                    if (ldgrel .and. mod(iGrel,nbproc) .ne. rang) cycle
                    call jaexin(jexnum(resuElem//'.RESL', iGrel), iexi)
                    if (iexi .eq. 0) cycle
                    call jeveuo(resuElem//'.DESC', 'L', vi=desc)
                    mode = desc(1+iGrel+1)
!
                    if (mode .gt. 0) then
                        nbNodeMode = nbno(mode)
                        nbElem = zznelg(ligrelNume,iGrel)
                        call jeveuo(jexnum(resuElem//'.RESL', iGrel), 'L', jresl)
                        nbCmpMode = digdel(mode)

! --------------------- Loop on elements
                        do iElem = 1, nbElem
                            elemNume = zzliel(ligrelNume,iGrel,iElem)
                            elemCoef = vectElemCoef

                            if (ldist .and. .not. ldgrel) then
                                if (elemNume .gt. 0) then
                                    if (numsd(elemNume) .ne. rang) cycle
                                else
                                    if (rang .ne. 0) cycle
                                endif
                            endif
                            if (elemNume .gt. 0) then
! ----------------------------- Physical element
                                iDof = 0
                                do iNodeMode = 1, nbNodeMode
                                    nbNode = connex(zi(iconx2+elemNume-1)+iNodeMode-1)
                                    iad1 = zi(idprn1-1+zi(idprn2+ ligrelMeshIndx-1)+&
                                              (nbNode-1)*(nec+2)+1-1)
                                    call corddl(admodl, lcmodl, idprn1, idprn2, ligrelMeshIndx,&
                                                mode, nec, nbCmp, nbNode, iNodeMode,&
                                                nbDofMode, zi(iapsdl))
                                    if (nbDofMode .eq. 0) cycle
                                    ASSERT(iad1 .ne. 0)
                                    ASSERT(iad1 .le. nbDof)
                                    ASSERT(nbDofMode .le. 100)
                                    if (vectScalType .eq. 1) then
                                        do iDofMode = 1, nbDofMode
                                            iDof = iDof+1
                                            zr(jvale-1+&
                                               zi(ianueq-1+iad1+zi(iapsdl-1+iDofMode)-1)) = &
                                            zr(jvale-1+zi(ianueq-1+iad1+zi(iapsdl-1+iDofMode)-1))+&
                                            zr(jresl+(iElem-1)*nbCmpMode+iDof-1)*elemCoef
                                        end do
!
                                    elseif (vectScalType .eq. 2) then
                                        do iDofMode = 1, nbDofMode
                                            iDof = iDof+1
                                            zc(jvale-1+&
                                                zi(ianueq-1+iad1+zi(iapsdl-1+iDofMode)-1)) = &
                                            zc(jvale-1+ zi(ianueq-1+iad1+zi(iapsdl-1+iDofMode)-1))+&
                                            zc(jresl+(iElem-1)*nbCmpMode+ iDof-1)*elemCoef
                                        end do
                                    else
                                        ASSERT(ASTER_FALSE)
                                    endif
                                end do
                            else
! ----------------------------- Lagrange element
                                elemNume = -elemNume
                                nbNode = zi(adne(1+3*(ligrelNume-1)+2)+elemNume)-&
                                         zi(adne(1+3*(ligrelNume-1)+2)+elemNume-1)-1
                                ASSERT(nbNodeMode .eq. nbNode)
                                iDof=0
                                do iNodeMode = 1, nbNodeMode
                                    nbNode = zi(adne(1+3*(ligrelNume-1)+1)-1+&
                                                zi(adne(1+3*(ligrelNume-1)+2)+elemNume-1)+&
                                                iNodeMode-1)
                                    if (nbNode .lt. 0) then
                                        nbNode = -nbNode
                                        if (liliNume .eq. 0) then
                                            call utmess('F', 'ASSEMBLA_45')
                                        endif
                                        ASSERT(liliNume .ne. 0)
                                        iad1 = zi(idprn1-1+&
                                                  zi(idprn2+liliNume-1)+&
                                                  (nbNode-1)*(nec+2)+1-1)
                                        call corddl(admodl, lcmodl, idprn1, idprn2, liliNume,&
                                                    mode, nec, nbCmp, nbNode, iNodeMode,&
                                                    nbDofMode, zi(iapsdl))
                                        ASSERT(nbDofMode .le. 100)
                                    else
                                        iad1 = zi(idprn1-1+&
                                                  zi(idprn2+ligrelMeshIndx-1)+&
                                                 (nbNode-1)*(nec+2)+1-1)
                                        call corddl(admodl, lcmodl, idprn1, idprn2,&
                                                    ligrelMeshIndx,&
                                                    mode, nec, nbCmp, nbNode, iNodeMode,&
                                                    nbDofMode, zi( iapsdl))
                                        ASSERT(nbDofMode .le. 100)
                                    endif
                                    ASSERT(iad1 .ne. 0)
                                    ASSERT(iad1 .le. nbDof)
                                    if (vectScalType .eq. 1) then
                                        do iDofMode = 1, nbDofMode
                                            iDof = iDof + 1
                                            zr(jvale-1+&
                                               zi(ianueq-1+iad1+zi(iapsdl-1+ iDofMode)-1)) =&
                                            zr(jvale-1+zi(ianueq-1+iad1+zi(iapsdl-1+iDofMode)-1))+&
                                            zr(jresl+(iElem-1)*nbCmpMode+iDof-1)*elemCoef
                                        end do
                                    else
                                        do iDofMode = 1, nbDofMode
                                            iDof = iDof + 1
                                            zc(jvale-1+&
                                               zi(ianueq-1+iad1+zi(iapsdl-1+ iDofMode)-1)) =&
                                            zc(jvale-1+zi(ianueq-1+iad1+zi(iapsdl-1+iDofMode)-1))+&
                                            zc(jresl+(iElem-1)*nbCmpMode+iDof-1)*elemCoef
                                        end do
                                    endif
                                end do
                            endif
                        end do
                        call jelibe(jexnum(resuElem//'.RESL', iGrel))
                    endif
                end do
            end do
        endif
    end do
!
    if (niv .ge. 2) then
        call uttcpu('CPU.ASSVEC', 'FIN', ' ')
        call uttcpr('CPU.ASSVEC', 7, temps)
        write (ifm, '(A44,D11.4,D11.4,D11.4)') &
            'TEMPS CPU/SYS/ELAPSED ASSEMBLAGE V        : ',&
            temps(5), temps(6), temps(7)
    endif

!
    if (ldist) call asmpi_comm_jev('MPI_SUM', vectValeJv)
!
270 continue

! - Elementary vector is a nodal field
    call asseVectField(vectAsse, numeDof, vectScalType,&
                       nbVectElem, listVectElem)

!
    if (dbg) then
        call dbgobj(vectValeJv, 'OUI', 6, '&&ASSVEC')
    endif
!
    call jedetr(vectAsse//'.LILI')
    call jedetr(vectAsse//'.LIVE')
    call jedetr(vectAsse//'.ADNE')
    call jedetr(vectAsse//'.ADLI')
    call jedetr('&&ASSVEC.POSDDL')
!
    if (.not.lparallel_mesh) then
        call asmpi_barrier()
    endif
    call uttcpu('CPU.CALC.1', 'FIN', ' ')
    call uttcpu('CPU.ASSE.1', 'FIN', ' ')
    call uttcpu('CPU.ASSE.3', 'FIN', ' ')
    call jedema()
end subroutine
