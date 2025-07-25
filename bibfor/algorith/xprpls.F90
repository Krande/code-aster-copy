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
subroutine xprpls(ligrel, dnoma, dcnsln, dcnslt, noma, &
                  cnsln, cnslt, grln, grlt, corres, &
                  ndim, ndomp, edomg)
    implicit none
!
#include "asterf_types.h"
#include "jeveux.h"
#include "asterc/r8maem.h"
#include "asterfort/assert.h"
#include "asterfort/calcul.h"
#include "asterfort/celces.h"
#include "asterfort/cescns.h"
#include "asterfort/cnscno.h"
#include "asterfort/cnsprj.h"
#include "asterfort/detrsd.h"
#include "asterfort/infmaj.h"
#include "asterfort/infniv.h"
#include "asterfort/jedema.h"
#include "asterfort/jedetr.h"
#include "asterfort/jelira.h"
#include "asterfort/jemarq.h"
#include "asterfort/jeveuo.h"
#include "asterfort/pj2dco.h"
#include "asterfort/pj3dco.h"
    character(len=8) :: dnoma, noma
    character(len=16) :: corres
    character(len=19) :: ligrel
    character(len=19) :: dcnsln, dcnslt, cnsln, cnslt, grln, grlt, ndomp, edomg
    integer(kind=8) :: ndim
!
! person_in_charge: samuel.geniaut at edf.fr
!
!     ------------------------------------------------------------------
!
!       XPRPLS   : X-FEM PROPAGATION : PROJECTION DES LEVEL SETS
!       ------     -     --            -              -     -
!
!  DANS LE CADRE DE LA PROPAGATION X-FEM AVEC LA METHODE FAST MARCHING UPWIND, ON PEUT
!  UTILISER DEUX MAILLAGES DIFFERENTS POUR LE MODELE PHYSIQUE ET POUR LA
!  REPRESENTATION DES LEVEL SETS. DANS CE CAS, ON DOIT PROJECTER
!  LES LEVEL SETS SUR LE MAILLAGE PHYSIQUE A PARTIR DU MAILLAGE UTILISE
!  POUR LA REPRESENTATION DES LEVEL SETS.
!
!  ENTREE
!  ------
!
!     ligrel = nom du ligrel construit sur l'ensemble des mailles
!              principales du maillage physique pour appel au calcul
!              du gradient des levels-set
!
!    * MODELE POUR LA REPRESENTATION DES LEVEL SETS
!      --------------------------------------------
!     DNOMA  = NOM DU MAILLAGE
!     DCNSLN = CHAMP_NO_S DES VALEURS DE LA LEVEL SET NORMALE
!     DCNSLT = CHAMP_NO_S DES VALEURS DE LA LEVEL SET TANGENTE
!     EDOMG  = VOIR XPRDOM.F
!
!
!    * MODELE PHYSIQUE
!      ---------------
!     NOMA   = NOME DU MAILLAGE
!     CNSLN  = CHAMP_NO_S DES VALEURS DE LA LEVEL SET NORMALE
!     CNSLT  = CHAMP_NO_S DES VALEURS DE LA LEVEL SET TANGENTE
!     GRLN   = CHAMP_NO_S DES VALEURS DU GRADIENT DE CNSLN
!     GRLT   = CHAMP_NO_S DES VALEURS DU GRADIENT DE CNSLT
!     NDOMP  = VOIR XPRDOM.F
!
!     NDIM   = DIMENSION DU MAILLAGE
!
!
!  SORTIE
!  ------
!     CORRES = NOM DU OBJET JEVEUX OU ON PEUT STOCKER LA CORRESPONDANCE
!              ENTRE LES DEUX MAILLAGES
!     CNSLN  = CHAMP_NO_S DES NOUVELLES VALEURS DE LA LEVEL SET NORMALE
!              POUR LE MODELE PHYSIQUE
!     CNSLT  = CHAMP_NO_S DES NOUVELLES VALEURS DE LA LEVEL SET TANGENTE
!              POUR LE MODELE PHYSIQUE
!     GRLN   = CHAMP_NO_S DES NOUVELLES VALEURS DU GRADIENT DE CNSLN
!              POUR LE MODELE PHYSIQUE
!     GRLT   = CHAMP_NO_S DES NOUVELLES VALEURS DU GRADIENT DE CNSLT
!              POUR LE MODELE PHYSIQUE
!
!     ------------------------------------------------------------------
!
!
!     CHARACTERISTICS OF THE MESHES
!
!     PROJECTION LEVEL SETS MESH
    integer(kind=8) :: nbelpr, jefrom, jcnlnl, jcnltl
!
!     PROJECTION PHYSICAL MESH
    character(len=19) :: tmplsn, tmplst
    integer(kind=8) :: jnto, nunopr
!
!     PROJECTION CODE
    aster_logical :: ldmax
    real(kind=8) :: distma
    character(len=8) :: lpain(4), lpaout(2)
    character(len=19) :: cnols, celgls, chams
    character(len=24) :: lchin(4), lchout(2)
!
!     GENERAL PURPOSE
    integer(kind=8) :: i, ibid
    integer(kind=8) :: ifm, niv
    real(kind=8), pointer :: cnlnv(:) => null()
    real(kind=8), pointer :: cnltv(:) => null()
    real(kind=8), pointer :: tmpln(:) => null()
    real(kind=8), pointer :: tmplt(:) => null()
!
!-----------------------------------------------------------------------
!     DEBUT
!-----------------------------------------------------------------------
    call jemarq()
    call infmaj()
    call infniv(ifm, niv)
!
    if (niv .ge. 0) then
        write (ifm, *)
        write (ifm, *) '   PROJECTION DES LEVEL SETS SUR LE MAILLAGE'//&
     &                  ' PHYSIQUE'
    end if
!
!     RETREIVE THE LIST OF THE NODES IN THE DOMAIN ON THE PHYSICAL MESH
    call jelira(ndomp, 'LONMAX', nunopr)
    call jeveuo(ndomp, 'L', jnto)
!
!     RETREIVE THE LIST OF THE ELEMENTS IN THE DOMAIN ON THE AUXILIARY
!     GRID
    call jelira(edomg, 'LONMAX', nbelpr)
    call jeveuo(edomg, 'L', jefrom)
!
!        PROJECT THE LEVELSETS FROM THE AUXILIARY MESH TO THE PHYSICAL
!        MESH. THE NODAL FIELD IS EXTRAPOLATED OUTSIDE THE AUXILIARY
!        MESH.
    ldmax = .false.
    distma = r8maem()
!
!        CREATE THE "CONNECTION" TABLE BETWEEN THE PHYSICAL AND
!        AUXILIARY MESHES
    if (ndim .eq. 2) then
        call pj2dco('PARTIE', dnoma, noma, nbelpr, zi(jefrom), &
                    nunopr, zi(jnto), ' ', ' ', corres, &
                    ldmax, distma, 0.d0)
    else
        call pj3dco('PARTIE', dnoma, noma, nbelpr, zi(jefrom), &
                    nunopr, zi(jnto), ' ', ' ', corres, &
                    ldmax, distma, 0.d0)
    end if
!
!        CREATE TWO TEMPORARY FIELDS WHERE THE PROJECTED VALUES WILL BE
!        STORED
    tmplsn = '&&OP0010.TMPLSN'
    tmplst = '&&OP0010.TMPLST'
!
!        PROJECTION OF THE NORMAL LEVELSET. THE EXISTING FIELD DATA
!        STRUCTURES ARE AUTOMATICALLY DESTROYED BY THE SUBROUTINE
!        "CNSPRJ"
    call cnsprj(dcnsln, corres, 'G', tmplsn, ibid)
    ASSERT(ibid .eq. 0)
!
!        PROJECTION OF THE TANGENTIAL LEVELSET. THE EXISTING FIELD DATA
!        STRUCTURES ARE AUTOMATICALLY DESTROYED BY THE SUBROUTINE
!        "CNSPRJ"
    call cnsprj(dcnslt, corres, 'G', tmplst, ibid)
    ASSERT(ibid .eq. 0)
!
    call jedetr(corres)
!
!        RETREIVE THE EXISTING NORMAL LEVEL SET FIELD
    call jeveuo(cnsln//'.CNSV', 'E', vr=cnlnv)
    call jeveuo(cnsln//'.CNSL', 'E', jcnlnl)
!
!        RETREIVE THE EXISTING TANGENTIAL LEVEL SET FIELD
    call jeveuo(cnslt//'.CNSV', 'E', vr=cnltv)
    call jeveuo(cnslt//'.CNSL', 'E', jcnltl)
!
!        RETREIVE THE TEMPORARY FIELDS WHERE THE PROJECTED VALUES OF THE
!        LEVEL SET HAVE BEEN STORED
    call jeveuo(tmplsn//'.CNSV', 'L', vr=tmpln)
    call jeveuo(tmplst//'.CNSV', 'L', vr=tmplt)
!
!        SUBSTITUTE THE PROJECTED VALUES OF THE LEVEL SETS INTO THE
!        EXISTING LEVEL SET FIELDS ONLY FOR THE NODES IN THE PHYSICAL
!        MESH INVOLVED IN THE PROJECTION
    do i = 1, nunopr
!
!           SUBSTITUTE THE PROJECTED VALUES OF THE NORMAL LEVEL SET INTO
!           THE EXISTING NORMAL LEVEL SET OF THE PHYSICAL MESH
        cnlnv(zi(jnto-1+i)) = tmpln(zi(jnto-1+i))
        zl(jcnlnl-1+zi(jnto-1+i)) = .true.
!
!           SUBSTITUTE THE PROJECTED VALUES OF THE TANGENTIAL LEVEL SET
!           INTO THE EXISTING TANGENTIAL LEVEL SET OF THE PHYSICAL MESH
        cnltv(zi(jnto-1+i)) = tmplt(zi(jnto-1+i))
        zl(jcnltl-1+zi(jnto-1+i)) = .true.
!
    end do
!
!        DESTROY THE TEMPORARY PROJECTED LEVEL SETS
    call detrsd('CHAM_NO_S', tmplsn)
    call detrsd('CHAM_NO_S', tmplst)
!
! ----------------------------------------------------------------------
!        CALCULATE THE GRADIENTS OF THE LEVEL SETS OF THE PHYSICAL MESH
! ----------------------------------------------------------------------
!
    if (niv .ge. 0) then
        write (ifm, *) '   CALCUL DES GRADIENTS DES LEVEL SETS SUR'//&
     &                  ' LE MAILLAGE PHYSIQUE'
    end if
!
!        NORMAL LEVEL SET
    cnols = '&&OP0010.GR.CNOLS'
    celgls = '&&OP0010.GR.CELGLS'
    chams = '&&OP0010.GR.CHAMS'
!
    call cnscno(cnsln, ' ', 'NON', 'V', cnols, &
                'F', ibid)
    lpain(1) = 'PGEOMER'
    lchin(1) = noma//'.COORDO'
    lpain(2) = 'PNEUTER'
    lchin(2) = cnols
    lpaout(1) = 'PGNEUTR'
    celgls = '&&OP0010.GR.CELGLS'
    lchout(1) = celgls
!
    call calcul('S', 'GRAD_NEUT_R', ligrel, 2, lchin, &
                lpain, 1, lchout, lpaout, 'V', &
                'OUI')
!
    call celces(celgls, 'V', chams)
    call cescns(chams, ' ', 'V', grln, ' ', &
                ibid)
!
    call detrsd('CHAM_NO', cnols)
    call detrsd('CHAM_ELEM', celgls)
    call detrsd('CHAM_ELEM_S', chams)
!
!        TANGENTIAL LEVEL SET
    call cnscno(cnslt, ' ', 'NON', 'V', cnols, &
                'F', ibid)
    lpain(1) = 'PGEOMER'
    lchin(1) = noma//'.COORDO'
    lpain(2) = 'PNEUTER'
    lchin(2) = cnols
    lpaout(1) = 'PGNEUTR'
    celgls = '&&OP0010.GR.CELGLS'
    lchout(1) = celgls
!
    call calcul('S', 'GRAD_NEUT_R', ligrel, 2, lchin, &
                lpain, 1, lchout, lpaout, 'V', &
                'OUI')
!
    call celces(celgls, 'V', chams)
    call cescns(chams, ' ', 'V', grlt, ' ', &
                ibid)
!
    call detrsd('CHAM_NO', cnols)
    call detrsd('CHAM_ELEM', celgls)
    call detrsd('CHAM_ELEM_S', chams)
!
!-----------------------------------------------------------------------
!     FIN
!-----------------------------------------------------------------------
    call jedema()
end subroutine
