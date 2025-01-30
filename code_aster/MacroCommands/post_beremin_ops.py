# coding=utf-8
# --------------------------------------------------------------------
# Copyright (C) 1991 - 2025 - EDF R&D - www.code-aster.org
# This file is part of code_aster.
#
# code_aster is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# code_aster is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with code_aster.  If not, see <http://www.gnu.org/licenses/>.
# --------------------------------------------------------------------
"""
Compute:
  - Weibull stress
  - Weibull stress**m
  - Probability of failure

  - Strain model -> Stress to consider (6 components: XX, YY, ZZ, XY, XZ, YZ)
      - "PETIT"      -> "SIEF"
      - "PETIT_REAC" -> "SIEF"
      - "GDEF_LOG"   -> "VARI_ELGA"
      - "SIMO_MIEHE" -> To retrieve Kirchhof stresses ? Not obvious ...
                        Normally, we can do without
                        Use of HHO models instead ?

Integral must be computed on the initial configuration

_biblio : CR-T66-2017-132, Lorentz
"""
from math import exp
import os
import numpy as np
import tempfile

import medcoupling as mc

from .Fracture.post_beremin_utils import CELL_TO_POINT
from ..Cata.Syntax import _F
from ..CodeCommands import CALC_CHAM_ELEM, CALC_CHAMP, CREA_TABLE, CALC_TABLE, DEFI_FICHIER
from ..CodeCommands import IMPR_RESU, AFFE_MODELE
from ..Objects import FieldOnCellsReal, NonLinearResult, Table
from ..Utilities import logger, disable_fpe, no_new_attributes
from ..Messages import UTMESS

DEBUG = bool(os.environ.get("DEBUG_POST_BEREMIN", ""))


class PostBeremin:
    """
    Beremin post-processor

    Arguments:
        resultat (NonLinearResult): Mechanical code_aster result
        group_ma (str): Mesh cells group on which Beremin post-treatment
            is carried out
            Only 1 mesh cells group authorized
        deformation (str): "PETIT", "PETIT_REAC", "GDEF_LOG"
        filtre_sigm (str): "SIGM_ELGA" or "SIGM_ELMOY".
        coef_mult (float): As explained in doc aster u4.81.22 (weibull) p.15

    Returns:
        Table:

            - Weibull stress
            - Weibull stress**m
            - Failure probability

    if SIGM_MAXI in command file, MED file containing:

        - where the plasticity is active : maximum of major principal
          stress
        - where the plasticity is not active : 0
          Result of type ELGA (apply ElgaFieldToSurface filter in
          Paraview to see it)
    """

    # data
    _result = _zone = _zone_ids = _stress_option = _strain_type = None
    _method_2D = _proj_3D_2D = _mesh_proj_mc = _mesh_proj = _mesh_3D_cells = _prec_proj = None
    _medfilename_temp = _unite_temp = None
    _intvar_idx = _stress_idx = None
    _use_hist = _use_indiplas = _use_function = None
    _coef_mult = None
    _weib_params = None
    # result to compute weibull stress
    _reswb = _rsieq = _rout = None
    __setattr__ = no_new_attributes(object.__setattr__)

    def __init__(self, result: NonLinearResult, strain_type: str, stress_option: str) -> None:
        """Initialization and definition of main arguments.

        Args:
            result (*Result*): Result object of the calculation.
            strain_type (str): Type of strain.
            stress_option (str): Type of stress to be used or calculated.
        """
        self._result = result
        assert stress_option in ("SIGM_ELGA", "SIGM_ELMOY"), f"unknown value: {stress_option}"
        self._stress_option = stress_option
        self._strain_type = strain_type

    def set_zone(self, group: str) -> None:
        """Define the location where the calculation will be done.

        Args:
            group (str): Mesh group of cells to be used.
        """
        self._zone = group
        self._zone_ids = self._result.getMesh().getCells(group)

    def set_projection_parameters(self, args):
        """Define projection parameters when METHODE_2D is used.

        Args:
            args (dict): Keywords for POST_BEREMIN.
        """
        if "METHODE_2D" in args:
            self._method_2D = True
            self._mesh_proj = args["METHODE_2D"]["MAILLAGE"]
            self._prec_proj = args["METHODE_2D"]["PRECISION"]

            ##ToFIX : ajout d'une verification que le résultat donné en entrée est bien 3D

    def set_indexes(self, intvar_idx: list = None, stress_idx: list = None) -> None:
        """Define the indexes used to extract the relevant components depending
        of the behaviour.

        Args:
            intvar_idx (list[int], optional): Indexes to access to
                INDIPLAS internal variable.
            stress_idx (list[int], optional): Indexes to extract log stresses
                tensor components.
        """
        if intvar_idx:
            if intvar_idx > 0:
                self._use_indiplas = True
                self._intvar_idx = intvar_idx
            else:
                self._use_indiplas = False
        if stress_idx:
            assert len(stress_idx) in (4, 6), stress_idx
            self._stress_idx = stress_idx

    def use_history(self, value: bool) -> None:
        """Enable the use of history or just the current timestep.

        Args:
            value (bool): *True* to use the history since the beginning,
              or *False* to use values on the current timestep.
        """
        self._use_hist = value

    def set_coef_mult(self, value: float = 1.0) -> None:
        """Define the multiplicative factor to take symmetric into account.

        Args:
            value (float): Multiplicative factor (see documentation).
        """
        self._coef_mult = value

    def set_weibull_parameters(self, args):
        """Define Weibull parameters.

        Args:
            args (dict): Keywords for POST_BEREMIN.
        """
        assert not (
            ("WEIBULL_FO" in args) and ("WEIBULL" in args)
        ), "Cannot use both WEIBULL and WEIBULL_FO"

        if "WEIBULL" in args:
            params = args["WEIBULL"]
            self._weib_params = params.copy()
        else:
            params = args["WEIBULL_FO"]
            self._use_function = True
            self._weib_params = params.copy()
            self._weib_params["SIGM_REFE"] = [params["SIGM_REFE"]]

        if self._rout:
            for param in ["SIGM_REFE", "M", "SIGM_SEUIL"]:
                assert len(self._weib_params[param]) in [
                    0,
                    1,
                ], "SIGM_MAXI can only be used with one set of WEIBULL parameter"

    def setup_calc_result(self):
        """Setup the result object to be used to compute the Beremin stress."""
        if self._strain_type == "GDEF_LOG":
            result = self._result
            indexes = result.getIndexes()
            # store VARI_ELGA + log stresses (from VARI_ELGA) as SIEF_ELGA
            dim = result.getMesh().getDimension()
            varn = [f"V{i}" for i in self._stress_idx]
            cmps = ("SIXX", "SIYY", "SIZZ", "SIXY", "SIXZ", "SIYZ")
            if dim == 2:
                cmps = ("SIXX", "SIYY", "SIZZ", "SIXY")
            dconv = dict(zip(varn, cmps))

            rvga = result
            rvga.userName = "rvga____"

            for idx in indexes:
                chvari = result.getField("VARI_ELGA", idx)
                rvga.setField(chvari.asPhysicalQuantity("SIEF_R", dconv), "SIEF_ELGA", idx)
            self._reswb = rvga
        else:
            self._reswb = self._result

    def get_internal_variables(self, idx):
        """Return internal variables.

        Args:
            idx (int): storage index.

        Returns:
            tuple: *ComponentOnCells* objects for INDIPLAS.
        """
        v1 = f"V{self._intvar_idx}"
        chvari = self._result.getField("VARI_ELGA", idx)
        chvari = chvari.toSimpleFieldOnCells()
        indip = chvari.getComponentOnCells(v1)

        return indip

    def get_major_stress(self, idx):
        """Return the major principal stress.

        Args:
            idx (int): storage index.

        Returns:
            *ComponentOnCells*: Major principal stress.
        """
        if not self._rsieq:
            self._rsieq = self.compute_sieq()
        sieq = self._rsieq.getField("SIEQ_ELGA", idx)
        sieq = sieq.toSimpleFieldOnCells()
        return sieq

    def compute_sieq(self):
        """Compute SIEQ_ELGA stress."""
        # compute fields on all timesteps, add INST/NUME_ORDRE to be limited
        input = self._reswb

        d_group_ma = {}
        if not self._rout:
            d_group_ma["GROUP_MA"] = self._zone

        if self._stress_option == "SIGM_ELMOY":
            input = CALC_CHAMP(RESULTAT=self._reswb, CONTRAINTE="SIMY_ELGA", **d_group_ma)
            for idx in input.getIndexes():
                simy = input.getField("SIMY_ELGA", idx)
                input.setField(simy, "SIEF_ELGA", idx)

        rsieq = CALC_CHAMP(RESULTAT=input, CRITERES="SIEQ_ELGA", **d_group_ma)

        if self._method_2D:

            self._medfilename_temp = tempfile.NamedTemporaryFile(dir=".", suffix=".med").name
            self._unite_temp = DEFI_FICHIER(
                FICHIER=self._medfilename_temp, ACTION="ASSOCIER", TYPE="LIBRE", ACCES="NEW"
            )
            IMPR_RESU(
                UNITE=self._unite_temp,
                PROC0="NON",
                RESU=_F(
                    RESULTAT=rsieq,
                    NOM_CHAM="SIEQ_ELGA",
                    NOM_CMP="PRIN_3",
                    NOM_CHAM_MED="TMP",
                    NUME_ORDRE=1,
                ),
            )
            DEFI_FICHIER(ACTION="LIBERER", UNITE=self._unite_temp)

        if DEBUG:
            for i in rsieq.getIndexes():
                print("NEW: rsieq:", i, np.sum(rsieq.getField("SIEQ_ELGA", i).getValues()))

        return rsieq

    def apply_stress_correction(self, sig1, idx, sigma_thr, sigma_refe):
        """Compute the major principal stress, eventually apply SIGM_CNV/SIGM_REFE
        correction and threshold by SIGM_SEUIL.

        *The values are changed in place.*

        Args:
            sig1 (*ComponentOnCells*): PRIN_3 component of SIEQ_ELGA, changed in place.
            idx (int): Timestep index
            cells_ids (list[int]): Cells list.
            sigma_thr (float): Stress threshold
            sigma_refe (float): sigm_refe parameter
        """

        if self._use_function and sig1.size > 0:
            # PRIN_3 / SIGM_REFE(TEMP) * SIGM_CNV
            temp = self.get_temperature_field(idx)
            sigref = temp.apply(sigma_refe)
            assert abs(sigref).min() > 0.0
            sig1 /= sigref
            sig1 *= self._weib_params["SIGM_CNV"]

        if sigma_thr > 0.0:
            # apply threshold
            sig1 = sig1 - sigma_thr

            def sig_filter(array):
                return np.where(array < 0.0, 0.0, array)

            sig1 = sig1.apply(sig_filter)

        return sig1

    def apply_threshold(self, sig1, indiplas):
        """Apply plastic strain threshold onto major stress.

        NB: *ComponentOnCells* changed in place.

        Args:
            sig1 (*ComponentOnCells*): Major stress.
            indiplas (*ComponentOnCells*): Plastic indicator.
        """
        indip = indiplas.onSupportOf(sig1)

        # * (1 if indip != 0 else 0)
        def toint_sign(array):
            return np.sign(np.array(array, dtype=int))

        indip = indip.apply(toint_sign)
        sig1 *= indip

        return sig1

    def get_temperature_field(self, idx):
        """Return the temperature field.

        Args:
            idx (idx): storage index.

        Returns:
            *ComponentOnCells*: Temperature values.
        """
        varc_elga = CALC_CHAMP(
            RESULTAT=self._result, VARI_INTERNE="VARC_ELGA", NUME_ORDRE=idx, GROUP_MA=self._zone
        )
        varc_elga = varc_elga.getField("VARC_ELGA", idx)
        varc_elga = varc_elga.toSimpleFieldOnCells()

        return varc_elga.TEMP

    def build_projector(self):
        """Return the 3D -> 2D projector when METHODE_2D is used

        Returns:
            *CELL_TO_POINT*: 3D -> 2D projector
        """

        ##Get 3D mesh (gauss to cells)
        timeStamps = mc.GetAllFieldIterations(self._medfilename_temp, "TMP")
        filefield_tmp = mc.MEDFileFieldMultiTS.New(self._medfilename_temp, "TMP")
        f_on_gauss = filefield_tmp.getFieldAtLevel(
            mc.ON_GAUSS_PT, timeStamps[-1][0], timeStamps[-1][1], 0
        )
        with disable_fpe():
            f_on_cells = f_on_gauss.voronoize(1e-12)
        self._mesh_3D_cells = f_on_cells.getMesh()
        if self._mesh_3D_cells.getMeshDimension() != 3:
            UTMESS("F", "RUPTURE4_14")

        ##Get 2D mesh
        ##TOFIX : découper le maillage 2D en cellules par points de Gauss
        self._mesh_proj_mc = self._mesh_proj.createMedCouplingMesh(spacedim=3)[0]
        if self._mesh_proj_mc.getMeshDimension() != 2:
            UTMESS("F", "RUPTURE4_15")

        ##Build projector 3D -> points (points = barycenters of 2D mesh cells)
        self._proj_3D_2D = CELL_TO_POINT(
            self._mesh_3D_cells, self._mesh_proj_mc, prec_rel=self._prec_proj
        )

    def compute_sig1_2D(self, sig1, pow_m, sigma_refe):
        """Projection of sig1 from 3D cells to barycenter of 2D cells

        Args:
            sig1 (*ComponentOnCells*): Major stress on 3D cells

        """

        ##TOFIX : aller récupérer les modèles affectés à _zone plutôt que de les inventer
        modele_tmp = AFFE_MODELE(
            MAILLAGE=self._result.getMesh().restrict(self._zone),
            AFFE=(
                _F(GROUP_MA=(self._zone), PHENOMENE="MECANIQUE", MODELISATION="3D"),
                _F(GROUP_MA=(self._zone), PHENOMENE="MECANIQUE", MODELISATION="3D_SI"),
            ),
        )

        ##3D ComponentOnCells to 3D medcoupling array
        sigmax = FieldOnCellsReal(modele_tmp, "ELGA", "SIEF_R")
        sigmax.setValues(0.0)
        sfield = sigmax.toSimpleFieldOnCells()
        sixx = sfield.SIXX
        sig1.restrict(self._zone_ids)
        sixx += sig1
        sfield.setComponentValues("SIXX", sixx.expand())
        sfield_mc = sfield.toMEDCouplingField(self._mesh_3D_cells)
        sfield_ar = sfield_mc.getArray()[:, 0]

        ##Projection 3D -> points
        sfield_ar_points = self._proj_3D_2D.Eval(sfield_ar)
        sfield_ar_points.setName("SIXX")

        ##Create 2D medcoupling field
        medc_cell_field = mc.MEDCouplingFieldDouble(mc.ON_CELLS, mc.ONE_TIME)
        medc_cell_field.setName("SIXX")
        if len(sfield_ar_points) == self._mesh_proj_mc.getNumberOfCells():
            medc_cell_field.setMesh(self._mesh_proj_mc)
        else:
            UTMESS("F", "RUPTURE4_16")
        medc_cell_field.setArray(sfield_ar_points**pow_m)
        medc_cell_field.setNature(mc.IntensiveConservation)
        medc_cell_field.setName("SIEF_ELGA")
        medc_cell_field.checkConsistencyLight()

        ##Compute table values
        coef_volu = self._coef_mult / self._weib_params["VOLU_REFE"]
        inv_m = 1 / pow_m
        sigma_u = self._weib_params["SIGM_CNV"] if self._use_function else sigma_refe

        ##INTESIXX
        intsig1pm = medc_cell_field.integral(0, True)

        # intfin(INTE_SIXX) = (COEF_MULT/VOLU_REFE*INTE_SIXX)**(1/bere_m)
        sigma_weibull = (coef_volu * intsig1pm) ** inv_m

        # sigwm(SIGMA_WEIBULL) = SIGMA_WEIBULL**bere_m
        sigma_weibullpm = sigma_weibull**pow_m

        # probaw(SIGMA_WEIBULL) = 1-exp(-SIGMA_WEIBULL**bere_m/sigma_u**bere_m)
        #  avec sigma_u=SIMG_CNV ou SIGM_REFE
        proba_weibull = 1.0 - exp(-((sigma_weibull / sigma_u) ** pow_m))

        return sigma_weibull, sigma_weibullpm, proba_weibull

    def main(self):
        """Compute Weibull stress and failure probability."""
        result = self._reswb
        model = result.getModel()
        params = result.getAccessParameters()
        logger.info("extracting integration scheme...")
        coor_elga = CALC_CHAM_ELEM(MODELE=model, OPTION="COOR_ELGA")
        coor_elga = coor_elga.toSimpleFieldOnCells()

        logger.info("starting computation...")

        table = TableBeremin(self._use_function)

        for sigma_thr in self._weib_params["SIGM_SEUIL"]:
            for sigma_refe in self._weib_params["SIGM_REFE"]:
                for pow_m in self._weib_params["M"]:

                    id_store = 0
                    previous = None

                    for idx, time in zip(params["NUME_ORDRE"], params["INST"]):

                        if self._use_indiplas:
                            indip = self.get_internal_variables(idx)

                        sieq = self.get_major_stress(idx)
                        sig1 = sieq.PRIN_3

                        if self._method_2D and not self._proj_3D_2D:
                            logger.info("building projector ...")
                            self.build_projector()

                        if self._use_indiplas:
                            sig1 = self.apply_threshold(sig1, indip)

                        sig1 = self.apply_stress_correction(sig1, idx, sigma_thr, sigma_refe)

                        if self._use_hist:
                            if previous is None:
                                previous = sig1.copy()
                            else:
                                sig1.maximum(previous)
                                previous = sig1.copy()

                        if DEBUG:
                            print("NEW: sigmax:", id_store, idx, sig1.sum())

                        if self._method_2D:
                            strwb, strwb_pm, proba = self.compute_sig1_2D(sig1, pow_m, sigma_refe)

                        else:
                            sig1 **= pow_m
                            self.store_sigm_maxi(id_store, time, sig1, model)

                            strwb, strwb_pm, proba = self.compute_table_values(
                                sig1, coor_elga.W, pow_m, sigma_refe
                            )

                        if self._use_function:
                            table.append(
                                id_store,
                                time,
                                self._zone,
                                pow_m,
                                "FONCTION",
                                sigma_thr,
                                strwb,
                                strwb_pm,
                                proba,
                            )
                        else:
                            table.append(
                                id_store,
                                time,
                                self._zone,
                                pow_m,
                                sigma_refe,
                                sigma_thr,
                                strwb,
                                strwb_pm,
                                proba,
                            )

                        id_store += 1

        return table

    def compute_table_values(self, sig1, weight, pow_m, sigma_refe):
        """Compute the values to be added into the result table.

        Args:
            sig1 (*ComponentOnCells*): Major stress ^ M.
            weight (*ComponentOnCells*): Weight of each integration point,
                **changed in place**.
            pow_m (float): M weibull parameter
            sigma_refe (float): sigm_refe parameter

        Returns:
            tuple: Values for each column.
        """
        coef_volu = self._coef_mult / self._weib_params["VOLU_REFE"]
        inv_m = 1 / pow_m
        sigma_u = self._weib_params["SIGM_CNV"] if self._use_function else sigma_refe

        logger.info("computing table values...")
        if DEBUG:
            print("NEW: sigpm:", sig1.sum())

        weight = weight.onSupportOf(sig1)
        # colonne INTE_SIXX
        # Restrict to GROUP_MA
        integ = sig1 * weight
        integ.restrict(self._zone_ids)
        intsig1pm = integ.sum()

        # intfin(INTE_SIXX) = (COEF_MULT/VOLU_REFE*INTE_SIXX)**(1/bere_m)
        sigma_weibull = (coef_volu * intsig1pm) ** inv_m

        # sigwm(SIGMA_WEIBULL) = SIGMA_WEIBULL**bere_m
        sigma_weibullpm = sigma_weibull**pow_m

        # probaw(SIGMA_WEIBULL) = 1-exp(-SIGMA_WEIBULL**bere_m/sigma_u**bere_m)
        #  avec sigma_u=SIMG_CNV ou SIGM_REFE
        proba_weibull = 1.0 - exp(-((sigma_weibull / sigma_u) ** pow_m))

        return sigma_weibull, sigma_weibullpm, proba_weibull

    def store_sigm_maxi(self, idx, time, sig1, model):
        """Store the major stress ^M into the output result.

        Args:
            idx (int): Storage index.
            time (float): Time value.
            sig1 (*ComponentOnCells*): Major stress ^ M.
            model (*Model*): Model object.
            cells_ids (list[int]): Cells list.
        """
        if not self._rout:
            return
        if idx == 0:
            self._rout.allocate(self._rsieq.getNumberOfIndexes())
        sigmax = FieldOnCellsReal(model, "ELGA", "SIEF_R")
        sigmax.setValues(0.0)
        sfield = sigmax.toSimpleFieldOnCells()
        sixx = sfield.SIXX
        sixx += sig1
        sfield.setComponentValues("SIXX", sixx.expand())
        fed = model.getFiniteElementDescriptor()
        sigmax = sfield.toFieldOnCells(fed, "RAPH_MECA", "PCONTPR")
        self._rout.setField(sigmax, "SIEF_ELGA", idx)
        self._rout.setModel(model, idx)
        self._rout.setTime(time, idx)


def post_beremin_ops(self, RESULTAT, GROUP_MA, DEFORMATION, FILTRE_SIGM, **args):
    """Main of POST_BEREMIN command."""
    post = PostBeremin(result=RESULTAT, strain_type=DEFORMATION, stress_option=FILTRE_SIGM)
    post.set_zone(GROUP_MA)
    post.set_indexes(intvar_idx=args["NUME_VARI"])
    if DEFORMATION == "GDEF_LOG":
        post.set_indexes(stress_idx=args["LIST_NUME_SIEF"])
    post.use_history(args["HIST_MAXI"] == "OUI")
    post.set_coef_mult(args["COEF_MULT"])
    post.set_projection_parameters(args)

    if args.get("SIGM_MAXI"):
        post._rout = NonLinearResult()
        self.register_result(post._rout, args["SIGM_MAXI"])

    post.set_weibull_parameters(args)
    post.setup_calc_result()

    table = post.main()
    sdtable = table.create_table()

    sdtable = CALC_TABLE(
        reuse=sdtable,
        TABLE=sdtable,
        ACTION=_F(
            OPERATION="TRI",
            NOM_PARA=("NUME_ORDRE", "M", "SIGM_REFE", "SIGM_SEUIL"),
            ORDRE="CROISSANT",
        ),
    )

    return sdtable


class TableBeremin(Table):
    """Helper object to build the result table."""

    types = None
    cols = [
        "NUME_ORDRE",
        "INST",
        "GROUP_MA",
        "M",
        "SIGM_REFE",
        "SIGM_SEUIL",
        "SIGMA_WEIBULL",
        "SIGMA_WEIBULL**M",
        "PROBA_WEIBULL",
    ]
    data = None
    __setattr__ = no_new_attributes(object.__setattr__)

    def __init__(self, use_function):
        Table.__init__(self)
        self.data = {}
        for col in self.cols:
            self.data.setdefault(col, [])
        if use_function:
            self.types = "IRKRKRRRR"
        else:
            self.types = "IRKRRRRRR"

    def append(self, *values):
        """Append the values of a row (columns order must match the definition).

        Args:
            values (list[misc]): Values of the row for each parameter.
        """
        assert len(values) == len(self.cols)
        for col, value in zip(self.cols, values):
            self.data[col].append(value)

    def create_table(self):
        """Create and return the Table datastructure.

        Returns:
            Table/table_sdaster: Table datastructure.
        """
        list_kws = []
        for col, typ in zip(self.cols, self.types):
            list_kws.append({"PARA": col, "LISTE_" + typ: self.data[col]})
        tab = CREA_TABLE(LISTE=list_kws)
        return tab
