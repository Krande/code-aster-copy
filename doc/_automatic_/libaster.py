# class AsterError in libaster

class AsterError(Exception):
    """Common base class for all non-exit exceptions.
    """
    
    # Method resolution order:
    #     AsterError
    #     builtins.Exception
    #     builtins.BaseException
    #     builtins.object
    
    # Data descriptors defined here:

# class ConvergenceError in libaster

class ConvergenceError(AsterError):
    """Common base class for all non-exit exceptions.
    """
    
    # Method resolution order:
    #     ConvergenceError
    #     AsterError
    #     builtins.Exception
    #     builtins.BaseException
    #     builtins.object

# class IntegrationError in libaster

class IntegrationError(AsterError):
    """Common base class for all non-exit exceptions.
    """
    
    # Method resolution order:
    #     IntegrationError
    #     AsterError
    #     builtins.Exception
    #     builtins.BaseException
    #     builtins.object

# class SolverError in libaster

class SolverError(AsterError):
    """Common base class for all non-exit exceptions.
    """
    
    # Method resolution order:
    #     SolverError
    #     AsterError
    #     builtins.Exception
    #     builtins.BaseException
    #     builtins.object

# class ContactError in libaster

class ContactError(AsterError):
    """Common base class for all non-exit exceptions.
    """
    
    # Method resolution order:
    #     ContactError
    #     AsterError
    #     builtins.Exception
    #     builtins.BaseException
    #     builtins.object

# class TimeLimitError in libaster

class TimeLimitError(AsterError):
    """Common base class for all non-exit exceptions.
    """
    
    # Method resolution order:
    #     TimeLimitError
    #     AsterError
    #     builtins.Exception
    #     builtins.BaseException
    #     builtins.object

# built-in function raiseAsterError in libaster

def raiseAsterError(idmess= 'VIDE_1'):
    pass

# class PythonBool in libaster

class PythonBool:
    """Members:
    
    NONE
    
    TRUE
    
    FALSE
    """
    
    # Method resolution order:
    #     PythonBool
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __eq__(self, other):
        pass
    
    def __getstate__(self):
        pass
    
    def __hash__(self):
        pass
    
    def __index__(self):
        pass
    
    def __init__(self, value):
        pass
    
    def __int__(self):
        pass
    
    def __ne__(self, other):
        pass
    
    def __repr__(self):
        pass
    
    def __setstate__(self, state):
        pass
    
    def name(self):
        pass
    
    #----------------------------------------------------------------------
    # Data descriptors defined here:
    
    @property
    def __members__(self):
        pass
    
    @property
    def name(self):
        """name(self: handle) -> str
        """
    
    @property
    def value(self):
        pass
    
    #----------------------------------------------------------------------
    # Data and other attributes defined here:
    
    FALSE = 0
    
    NONE = -1
    
    TRUE = 1

# class DataStructure in libaster

class DataStructure:
    pass
    
    # Method resolution order:
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __init__(self,  *args, **kwargs):
        """Initialize self.  See help(type(self)) for accurate signature.
        """
    
    def addDependency(self, ds):
        """Add a dependency to a *DataStructure*.
        
        Arguments:
            ds (*DataStructure*): Parent *DataStructure* to depend on.
        """
    
    def build(self):
        """Update the *DataStructure* attributes from the *Jeveux* objects.
        *Only use internally after calling fortran subroutines*.
        
        Returns:
            bool: *True* if all went ok, *False* otherwise.
        """
    
    def debugPrint(self, unit= 6):
        """Print the raw content of a *DataStructure* on the selected file.
        
        Args:
            unit (int): File number (default: 6, means stdout).
        """
    
    def getDependencies(self):
        """Return the explicit dependencies.
        
        Returns:
            list[*DataStructure*]: List of parents (dependencies) *DataStructure*.
        """
    
    def getName(self):
        """Return the internal (*Jeveux*) name of the *DataStructure*.
        
        Returns:
            str: Internal/*Jeveux* name.
        """
    
    def getTitle(self):
        """Return the tile of the *DataStructure* .
        
        Returns:
            str: Title of the *DataStructure*.
        """
    
    def getType(self):
        """Return the name of the *DataStructure* type.
        
        Returns:
            str: Name of the *DataStructure* type.
        """
    
    def id(self):
        """Return the identity of the object.
        
        Returns:
            int: Identifier (address as int).
        """
    
    def removeDependency(self, ds):
        """Remove a dependency to a *DataStructure*.
        
        Arguments:
            ds (*DataStructure*): Parent *DataStructure* to be removed from
                dependencies.
        """
    
    def resetDependencies(self):
        """Clear the list of explicit dependencies.
        """
    
    def setTitle(self, title):
        """Set the tile of the *DataStructure* .
        
        Arguments:
            title [str]: Title of the *DataStructure*.
        """
    
    #----------------------------------------------------------------------
    # Data descriptors defined here:
    
    @property
    def ptr_sdj(self):
        pass
    
    @property
    def userName(self):
        """str: Name of the user variable that holds this object.
        """

# built-in function debugJeveuxContent in libaster

def debugJeveuxContent(arg0):
    pass

# built-in function debugJeveuxExists in libaster

def debugJeveuxExists(arg0):
    pass

# built-in function use_count in libaster

def use_count(*args, **kwargs):
    """Overloaded function.
    
    1. use_count(arg0: Mesh) -> int
    
    2. use_count(arg0: Model) -> int
    
    3. use_count(arg0: DOFNumbering) -> int
    
    4. use_count(arg0: ElementaryMatrix<double, (PhysicalQuantityEnum)4>) -> int
    
    5. use_count(arg0: ElementaryMatrix<std::complex<double>, (PhysicalQuantityEnum)4>) -> int
    
    6. use_count(arg0: ElementaryMatrix<double, (PhysicalQuantityEnum)6>) -> int
    
    7. use_count(arg0: ElementaryMatrix<std::complex<double>, (PhysicalQuantityEnum)5>) -> int
    
    8. use_count(arg0: AssemblyMatrix<double, (PhysicalQuantityEnum)4>) -> int
    
    9. use_count(arg0: AssemblyMatrix<std::complex<double>, (PhysicalQuantityEnum)4>) -> int
    
    10. use_count(arg0: AssemblyMatrix<double, (PhysicalQuantityEnum)6>) -> int
    
    11. use_count(arg0: AssemblyMatrix<std::complex<double>, (PhysicalQuantityEnum)6>) -> int
    
    12. use_count(arg0: AssemblyMatrix<double, (PhysicalQuantityEnum)5>) -> int
    
    13. use_count(arg0: AssemblyMatrix<std::complex<double>, (PhysicalQuantityEnum)5>) -> int
    """

# class EntityType in libaster

class EntityType:
    """Members:
    
    GroupOfNodesType
    
    GroupOfCellsType
    
    AllMeshEntitiesType
    
    CellType
    
    NodeType
    
    NoType
    """
    
    # Method resolution order:
    #     EntityType
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __eq__(self, other):
        pass
    
    def __getstate__(self):
        pass
    
    def __hash__(self):
        pass
    
    def __index__(self):
        pass
    
    def __init__(self, value):
        pass
    
    def __int__(self):
        pass
    
    def __ne__(self, other):
        pass
    
    def __repr__(self):
        pass
    
    def __setstate__(self, state):
        pass
    
    def name(self):
        pass
    
    #----------------------------------------------------------------------
    # Data descriptors defined here:
    
    @property
    def __members__(self):
        pass
    
    @property
    def name(self):
        """name(self: handle) -> str
        """
    
    @property
    def value(self):
        pass
    
    #----------------------------------------------------------------------
    # Data and other attributes defined here:
    
    AllMeshEntitiesType = 2
    
    CellType = 3
    
    GroupOfCellsType = 1
    
    GroupOfNodesType = 0
    
    NoType = 5
    
    NodeType = 4

# class MeshEntity in libaster

class MeshEntity:
    pass
    
    # Method resolution order:
    #     MeshEntity
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __init__(self, arg0, arg1):
        pass
    
    def getNames(self):
        pass
    
    def getType(self):
        pass

# class AllMeshEntities in libaster

class AllMeshEntities(MeshEntity):
    pass
    
    # Method resolution order:
    #     AllMeshEntities
    #     MeshEntity
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __init__(self,  *args, **kwargs):
        """Initialize self.  See help(type(self)) for accurate signature.
        """

# class BaseMesh in libaster

class BaseMesh(DataStructure):
    pass
    
    # Method resolution order:
    #     BaseMesh
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __init__(self,  *args, **kwargs):
        """Initialize self.  See help(type(self)) for accurate signature.
        """
    
    def build(self):
        """Build list of Tables based on the mesh
        
        Returns:
            bool: true if building is ok
        """
    
    def getCellName(self, index):
        """Return the name of the given cell
        
        Arguments:
            index (int) : index of the cell
        
        Returns:
            str : name of the cell (stripped)
        """
    
    def getConnectivity(self):
        """Return the connectivity of the mesh as Python lists.
        
        Returns:
            list[list[int]]: List of, for each cell, a list of the nodes indexes.
        """
    
    def getCoordinates(self):
        """Return the coordinates of the mesh.
        
        Returns:
            MeshCoordinatesField: Field of the coordinates.
        """
    
    def getDimension(self):
        """Return the dimension of the mesh.
        
        Returns:
            int: 2 or 3
        """
    
    def getMedCellsTypes(self):
        """Return the Med type of each cell.
        
        Returns:
            list[int]: List of Med types.
        """
    
    def getMedConnectivity(self):
        """Return the connectivity of the mesh as Python lists following the Med numbering.
        
        Returns:
            list[list[int]]: List of, for each cell, a list of the nodes indexes.
        """
    
    def getNodeName(self, index):
        """Return the name of the given node
        
        Arguments:
            index (int) : index of the node
        
        Returns:
            str : name of the node (stripped)
        """
    
    def getNumberOfCells(self):
        """Return the number of cells of the mesh.
        
        Returns:
            int: Number of cells.
        """
    
    def getNumberOfNodes(self):
        """Return the number of nodes of the mesh.
        
        Returns:
            int: Number of nodes.
        """
    
    def getTable(self, identifier):
        """Extract a Table from the datastructure.
        
        Arguments:
            identifier (str): Table identifier.
        
        Returns:
            Table: Table stored with the given identifier.
        """
    
    def isParallel(self):
        """Tell if the mesh is distributed on parallel instances.
        
        Returns:
            bool: *False* for a centralized mesh, *True* for a parallel mesh.
        """
    
    def printMedFile(self, fileName, local= True):
        """Print the mesh in the MED format
        
        Arguments:
            filename (str): Name of the file
            local (bool=True) : print local values only (relevent for ParallelMesh only)
        
        Returns:
            Bool: True if of
        """
    
    def update(self):
        """Update the internal state of the datastructure.
        
        Returns:
            bool: *True* in case of success, *False* otherwise.
        """

# class Mesh in libaster

class Mesh(BaseMesh):
    pass
    
    # Method resolution order:
    #     Mesh
    #     BaseMesh
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __init__(self, *args, **kwargs):
        """Overloaded function.
        
        1. __init__(self: libaster.Mesh) -> None
        
        2. __init__(self: libaster.Mesh, arg0: str) -> None
        """
    
    def getCells(self, group_name= ''):
        """Return the list of the indexes of the cells that belong to a group of cells.
        
        Arguments:
            group_name (str): Name of the local group.
        
        Returns:
            list[int]: Indexes of the cells of the local group.
        """
    
    def getGroupsOfCells(self, local= False):
        """Return the list of the existing groups of cells.
        
        Returns:
            list[str]: List of groups names (stripped).
        """
    
    def getGroupsOfNodes(self, local= False):
        """Return the list of the existing groups of nodes.
        
        Arguments:
            local=false (bool): not used (for compatibilty with ParallelMesh)
        
        Returns:
            list[str]: List of groups names (stripped).
        """
    
    def getInnerNodes(self):
        """Return the list of the indexes of the nodes in the mesh
        
        Returns:
            list[int]: Indexes of the nodes.
        """
    
    def hasGroupOfCells(self, group_name, local= False):
        """The group exists in the mesh
        
        Arguments:
            group_name (str): Name of the group.
            local=false (bool): not used (for compatibilty with ParallelMesh)
        
        Returns:
            bool: *True* if exists, *False* otherwise.
        """
    
    def hasGroupOfNodes(self, group_name, local= False):
        """The group exists in the mesh
        
        Arguments:
            group_name (str): Name of the group.
            local=false (bool): not used (for compatibilty with ParallelMesh)
        
        Returns:
            bool: *True* if exists, *False* otherwise.
        """
    
    def isQuadratic(self):
        """To know if the mesh contains quadratic cells
        
        Returns:
            bool: *True* if the mesh contains quadratic cells, *False* otherwise.
        """
    
    def readAsterFile(self, filename):
        """Read a mesh file from ASTER format.
        
        Arguments:
            filename (str): Path to the file to be read.
        
        Returns:
            bool: *True* if succeeds, *False* otherwise.
        """
    
    def readGibiFile(self, filename):
        """Read a mesh file from GIBI format.
        
        Arguments:
            filename (str): Path to the file to be read.
        
        Returns:
            bool: *True* if succeeds, *False* otherwise.
        """
    
    def readGmshFile(self, filename):
        """Read a mesh file from GMSH format.
        
        Arguments:
            filename (str): Path to the file to be read.
        
        Returns:
            bool: *True* if succeeds, *False* otherwise.
        """
    
    def readMedFile(self, filename):
        """Read a mesh file from MED format.
        
        Arguments:
            filename (str): Path to the file to be read.
        
        Returns:
            bool: *True* if succeeds, *False* otherwise.
        """

# built-in function getMedCouplingConversionData in libaster

def getMedCouplingConversionData(mesh):
    """Return three dictionnaries containing data to create an equivalent MedCoupling unstructured mesh.
    
    MedCoupling needs a mesh splitted by dimension for what concerns cells and groups of cells.
    The group of nodes all belongs to an unique level so there is no need to split them.
    
     - The first dictionnary (cells) contains for each dimension (the keys) :
       1. The connectivity
       2. The connectivity index
     - The second dictionnary (groups_c) contains for each dimension (the keys) a dictionnary
       which keys are the groups names at the items their cells.
     - The third dictionnary (groups_n) contains for each group of nodes (the keys)
       the nodes composing the group.
    
    Arguments:
        mesh (BaseMeshPtr): The aster mesh to be processed.
    
    Returns:
        tuple (cells, groups_c, groups_n) : The data to create the equivalent MedCoupling mesh.
    """

# class DiscreteComputation in libaster

class DiscreteComputation:
    pass
    
    # Method resolution order:
    #     DiscreteComputation
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __init__(self, arg0):
        pass
    
    def getCompressibilityMatrix(self, groupOfCells= []):
        """Return the elementary matrices for compressibility acoustic matrix.
        Option MASS_ACOU.
        
        Arguments:
            groupOfCells (list[str]): compute matrices on given groups of cells.
                If it empty, the full model is used
        Returns:
            ElementaryMatrix: elementary mass matrix
        """
    
    def getContactForces(self, geom, displ, displ_step, time_prev, time_step, data, coef_cont, coef_frot):
        """Compute contact and friction forces
        
        Arguments:
            geom (MeshCoordinatesField): coordinates of mesh used to compute normal
            displ (FieldOnNodes): displacement field at begin of current time
            displ_step (FieldOnNodes): field of increment of displacement
            time_prev (float): time at begin of the step
            time_curr (float): delta time between begin and end of the step
            data (FieldOnCellsReal): contact data
            coef_cont (FeildOnNodesReal) : contact coefficient
            coef_frot (FeildOnNodesReal) : friction coefficient
        
        Returns:
            FieldOnNodesReal: contact and friction forces
        """
    
    def getContactMatrix(self, geom, displ, displ_step, time_prev, time_step, data, coef_cont, coef_frot):
        """Compute contact matrix
        
        Arguments:
            geom (MeshCoordinatesField): coordinates of mesh used to compute normal
            displ (FieldOnNodes): displacement field at begin of current time
            displ_step (FieldOnNodes): field of increment of displacement
            time_prev (float): time at begin of the step
            time_curr (float): delta time between begin and end of the step
            data (FieldOnCellsReal): contact data
            coef_cont (FeildOnNodesReal) : contact coefficient
            coef_frot (FeildOnNodesReal) : friction coefficient
        
        Returns:
            ElementaryMatrixDisplacementReal: contact and friction elementary matrix
        """
    
    def getDirichletBC(self, time):
        """Return the imposed displacement vector used to remove imposed DDL
        
        Arguments:
              time (float): Current time
        
        Returns:
              FieldOnNodes: imposed displacement vector
        """
    
    def getDualDisplacement(self, disp_curr, scaling= 1.0):
        """Return the Dirichlet load vector
        
        Arguments:
              disp_curr (FieldOnNodes): current displacement vector
        
        Returns:
              FieldOnNodes: Dirichlet load vector
        """
    
    def getDualElasticStiffnessMatrix(self):
        """Return elementary matrices for dual mechanical BC
        
        Returns:
            ElementaryMatrix: elementary matrices
        """
    
    def getDualForces(self, disp_curr):
        """Return the imposed displacement assembled vector
        
        Arguments:
              disp_curr (FieldOnNodes): current displacement vector
        
        Returns:
              FieldOnNodes: dual reaction vector (B^T*lambda)
        """
    
    def getDualLinearConductivityMatrix(self):
        """Return elementary matrices for dual thermal BC
        
        Returns:
            ElementaryMatrix: elementary matrices
        """
    
    def getDualLinearMobilityMatrix(self):
        """Return elementary matrices for dual acoustic BC
        
        Returns:
            ElementaryMatrix: elementary matrices
        """
    
    def getElasticStiffnessMatrix(self, time= 0.0, fourierMode= -1, groupOfCells= [], with_dual= True):
        """Return the elementary matrices for elastic Stiffness matrix.
        Option RIGI_MECA.
        
        Arguments:
              time (float): Current time for external state variable evaluation (default: 0.0)
              fourierMode (int): Fourier mode (default: -1)
              groupOfCells (list[str]): compute matrices on given groups of cells.
                  If it empty, the full model is used
              with_dual (bool): compute dual terms or not (default: True)
        Returns:
              ElementaryMatrix: elementary elastic Stiffness matrix
        """
    
    def getExchangeThermalMatrix(self, time):
        """Return the elementary matices for exhange thermal matrix.
        
        Arguments:
            time (float): Current time
        Returns:
            ElementaryMatrix: elementary exchange thermal matrices
        """
    
    def getExternalStateVariablesForces(self, time):
        """Compute load from external state variables
        
        Arguments:
              time (float): Current time
        
        Returns:
              FieldOnNodes: load from external state variables
        """
    
    def getFluidStructureMassMatrix(self, time= 0.0, groupOfCells= []):
        """Return the elementary matrices for fluid-structure mass matrix.
        Option MASS_FLUI_STRUC.
        
        Arguments:
              time (float): Current time for external state variable evaluation (default: 0.0)
              groupOfCells (list[str]): compute matrices on given groups of cells.
                  If it empty, the full model is used
        Returns:
              ElementaryMatrixReal: elementary fluid-structure mass matrix
        """
    
    def getFluidStructureStiffnessMatrix(self, time= 0.0, fourierMode= -1, groupOfCells= []):
        """Return the elementary matrices for fluid-structure stiffness matrix.
        Option RIGI_FLUI_STRUC.
        
        Arguments:
              time (float): Current time for external state variable evaluation (default: 0.0)
              fourierMode (int): Fourier mode (default: -1)
              groupOfCells (list[str]): compute matrices on given groups of cells.
                  If it empty, the full model is used
        Returns:
              ElementaryMatrixReal: elementary fluid-structure Stiffness matrix
        """
    
    def getGeometricStiffnessMatrix(self, sief_elga, strx_elga= None, displ= None, modeFourier= -1, groupOfCells= []):
        """Return the elementary matrices for geometric Stiffness matrix.
        Option RIGI_MECA_HYST.
        
        Arguments:
            sief_elga (FieldOnCellsReal) : stress at Gauss points
            strx_elga (FieldOnCellsReal) : stress at Gauss points for structural element
            displ (FieldOnNodesReal) : displacement field
            groupOfCells (list[str]): compute matrices on given groups of cells.
                If it empty, the full model is used
        Returns:
            ElementaryMatrixComplex: elementary geometric rigidity matrix
        """
    
    def getGyroscopicDampingMatrix(self, groupOfCells= []):
        """Return the elementary matrices for gyroscopic damping matrix.
        Option MECA_GYRO.
        
        Arguments:
            groupOfCells (list[str]): compute matrices on given groups of cells.
                If it empty, the full model is used
        Returns:
            ElementaryMatrixReal: elementary gyroscopic damping matrix
        """
    
    def getGyroscopicStiffnessMatrix(self, groupOfCells= []):
        """Return the elementary matrices for gyroscopic Stiffness matrix.
        Option RIGI_GYRO.
        
        Arguments:
            groupOfCells (list[str]): compute matrices on given groups of cells.
                If it empty, the full model is used
        Returns:
            ElementaryMatrixReal: elementary gyroscopic rigidity matrix
        """
    
    def getHystereticStiffnessMatrix(self, stiffnessMatrix, time= 0.0, groupOfCells= []):
        """Return the elementary matrices for viscoelastic Stiffness matrix.
        Option RIGI_MECA_HYST.
        
        Arguments:
            stiffnessMatrix : elementary stiffness matrix
            time (float): Current time for external state variable evaluation (default: 0.0)
            groupOfCells (list[str]): compute matrices on given groups of cells.
                If it empty, the full model is used
        Returns:
            ElementaryMatrixComplex: elementary viscoelastic rigidity matrix
        """
    
    def getImpedanceBoundaryMatrix(self, groupOfCells= []):
        """Return the elementary matrices for impedance (mechanical) matrix.
        Option IMPE_MECA.
        
        Returns:
            ElementaryMatrixReal: impedance mechanical matrix
        """
    
    def getImpedanceMatrix(self):
        """Return the elementary matrices for impedance (acoustic) damping matrix.
        Option AMOR_ACOU.
        
        Returns:
            ElementaryMatrixReal: elementary damping matrix
        """
    
    def getImpedanceWaveMatrix(self, groupOfCells= []):
        """Return the elementary matrices for impedance (mechanical) matrix
        from an harmonic wave.
        Option ONDE_FLUI.
        
        Returns:
            ElementaryMatrixReal: impedance wave matrix
        """
    
    def getImposedDualBC(self, *args, **kwargs):
        """Overloaded function.
        
        1. getImposedDualBC(self: libaster.DiscreteComputation, time: float, time_step: float, theta: float) -> FieldOnNodes<double>
        
        
              Return the imposed nodal BC assembled vector
        
              Arguments:
                    time (float): Current time
                    time_step (float): Time increment
                    theta (float): Theta parameter for integration
        
              Returns:
                    FieldOnNodes: imposed dual field
                
        
        2. getImposedDualBC(self: libaster.DiscreteComputation, time: float) -> FieldOnNodes<double>
        
        
              Return the imposed nodal BC assembled vector
        
              Arguments:
                    time (float): Current time
        
              Returns:
                    FieldOnNodes: imposed dual field
        """
    
    def getIncrementalDirichletBC(self, time, disp):
        """Return the incremental imposed displacement vector used to remove imposed DDL
        for incremental resolution.
        
        incr_disp = getDirichletBC(time) - disp, with 0.0 for DDL not imposed
        
        Arguments:
              time (float): Current time
              disp (FieldOnNodes): displacement field at current time
        
        Returns:
              FieldOnNodes: incremental imposed displacement vector
        """
    
    def getInternalForces(self, displ, displ_step, stress, internVar, time_prev, time_step, groupOfCells= []):
        """Compute internal forces (integration of behaviour)
        
        Arguments:
            displ (FieldOnNodes): displacement field at begin of current time
            displ_step (FieldOnNodes): field of increment of displacement
            stress (FieldOnCells): field of stress at begin of current time
            internVar (FieldOnCells): field of internal state variables at begin of current time
            time_prev (float): time at begin of the step
            time_curr (float): delta time between begin and end of the step
            groupOfCells (list[str]): compute matrices on given groups of cells.
        
        Returns:
            tuple (tuple): return code error (FieldOnCells),
            error code flag (integer),
            internal state variables VARI_ELGA (FieldOnCells),
            Cauchy stress SIEF_ELGA (FieldOnCells),
            field of internal forces (FieldOnNodesReal),
        """
    
    def getLinearCapacityMatrix(self, time, groupOfCells= []):
        """Return the elementary matrices for linear Capacity matrix in thermal computation.
        Option MASS_THER.
        
        Arguments:
            time (float): current time to evaluate rho_cp
            groupOfCells (list[str]): compute matrices on given groups of cells.
                If it empty, the full model is used
        Returns:
            ElementaryMatrix: elementary mass matrix
        """
    
    def getLinearConductivityMatrix(self, time, fourierMode= 0, groupOfCells= [], with_dual= True):
        """Return the elementary matices for linear thermal matrix.
        Option RIGI_THER.
        
        Arguments:
              time (float): Current time
              fourierMode (int): Fourier mode (default: -1)
              groupOfCells (list[str]): compute matrices on given groups of cells.
                If it empty, the full model is used
              with_dual (bool): compute dual terms or not (default: True)
        Returns:
              ElementaryMatrix: elementary linear thermal matrices
        """
    
    def getLinearMobilityMatrix(self, groupOfCells= [], with_dual= True):
        """Return the elementary matices for linear mobility acoustic matrix
        Option RIGI_ACOU.
        
        Arguments:
            groupOfCells (list[str]): compute matrices on given groups of cells.
            with_dual (bool): compute dual terms or not (default: True)
        
        Returns:
            ElementaryMatrix: elementary linear acoustic matrices
        """
    
    def getMechanicalDampingMatrix(self, getMechanicalMassMatrix= None, stiffnessMatrix= None, time= 0.0, groupOfCells= []):
        """Return the elementary matrices for damping matrix.
        Option AMOR_MECA.
        
        Arguments:
            getMechanicalMassMatrix : elementary mass matrix
            stiffnessMatrix : elementary stiffness matrix
            time (float): Current time for external state variable evaluation (default: 0.0)
            groupOfCells (list[str]): compute matrices on given groups of cells.
                If it empty, the full model is used
        Returns:
            ElementaryMatrixReal: elementary damping matrix
        """
    
    def getMechanicalMassMatrix(self, diagonal, time= 0.0, groupOfCells= []):
        """Return the elementary matrices for mechanical mass matrix
        Option MASS_MECA.
        
        Arguments:
            diagonal (bool) : True for diagonal mass matrix else False.
            time (float): Current time for external state variable evaluation (default: 0.0)
            groupOfCells (list[str]): compute matrices on given groups of cells.
                If it empty, the full model is used
        Returns:
            ElementaryMatrix: elementary mass matrix
        """
    
    def getNeumannForces(self, time= 0.0, time_step= 0.0, theta= 1.0, previousPrimalField= None):
        """Return the Neumann forces vector
        
        Arguments:
              time (float): Current time
              time_step (float): Time increment
              theta (float): Theta parameter for time-integration
              previousPrimalField (fieldOnNodesReal): solution field at previous time
        
        Returns:
              FieldOnNodes: Neumann forces vector
        """
    
    def getPhysicalProblem(self):
        """Get physical probelm
        
        Returns:
              PhysicalProblem: physical problem
        """
    
    def getPredictionTangentStiffnessMatrix(self, displ, displ_step, stress, internVar, time_prev, time_step, groupOfCells= []):
        """Compute jacobian matrix for Newton algorithm, Euler prediction
        
        Arguments:
            displ (FieldOnNodes): displacement field at begin of current time
            displ_step (FieldOnNodes): field of increment of displacement
            stress (FieldOnCells): field of stress at begin of current time
            internVar (FieldOnCells): field of internal state variables at begin of current time
            time_prev (float): time at begin of the step
            time_curr (float): delta time between begin and end of the step
            groupOfCells (list[str]): compute matrices on given groups of cells.
        
        Returns:
            tuple (tuple): return code error (FieldOnCellsLong),
            error code flag (int),
            elementary tangent matrix (ElementaryMatrixDisplacementReal),
        """
    
    def getRotationalStiffnessMatrix(self, groupOfCells= []):
        """Return the elementary matrices for rotational Stiffness matrix.
        Option RIGI_ROTA.
        
        Arguments:
            groupOfCells (list[str]): compute matrices on given groups of cells.
                If it empty, the full model is used
        Returns:
            ElementaryMatrixReal: elementary rotational rigidity matrix
        """
    
    def getTangentStiffnessMatrix(self, displ, displ_step, stress, internVar, time_prev, time_step, groupOfCells= []):
        """Compute jacobian matrix for Newton algorithm
        
        Arguments:
            displ (FieldOnNodes): displacement field at begin of current time
            displ_step (FieldOnNodes): field of increment of displacement
            stress (FieldOnCells): field of stress at begin of current time
            internVar (FieldOnCells): field of internal state variables at begin of current time
            time_prev (float): time at begin of the step
            time_curr (float): delta time between begin and end of the step
            groupOfCells (list[str]): compute matrices on given groups of cells.
        
        Returns:
            tuple (tuple): return code error (FieldOnCellsLong),
            error code flag (int),
            elementary tangent matrix (ElementaryMatrixDisplacementReal)
        """
    
    def getTransientThermalForces(self, time, time_step, theta, previousPrimalField= None):
        """Compute Transient Thermal Load
        
        Arguments:
              time (float): Current time
              time_step (float): Time increment
              theta (float): Theta parameter for integration
              previousPrimalField (fieldOnNodesReal): solution field at previous time
        
        Returns:
              FieldOnNodes: load from external state variables
        """

# class FieldOnNodesDescription in libaster

class FieldOnNodesDescription(DataStructure):
    pass
    
    # Method resolution order:
    #     FieldOnNodesDescription
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __init__(self, *args, **kwargs):
        """Overloaded function.
        
        1. __init__(self: libaster.FieldOnNodesDescription) -> None
        
        2. __init__(self: libaster.FieldOnNodesDescription, arg0: str) -> None
        """

# class BaseDOFNumbering in libaster

class BaseDOFNumbering(DataStructure):
    pass
    
    # Method resolution order:
    #     BaseDOFNumbering
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __init__(self,  *args, **kwargs):
        """Initialize self.  See help(type(self)) for accurate signature.
        """
    
    def addDirichletBC(self, *args, **kwargs):
        """Overloaded function.
        
        1. addDirichletBC(self: libaster.BaseDOFNumbering, arg0: DirichletBC) -> None
        
        2. addDirichletBC(self: libaster.BaseDOFNumbering, arg0: DirichletBC, arg1: Function) -> None
        
        3. addDirichletBC(self: libaster.BaseDOFNumbering, arg0: DirichletBC, arg1: Formula) -> None
        
        4. addDirichletBC(self: libaster.BaseDOFNumbering, arg0: DirichletBC, arg1: Function2D) -> None
        """
    
    def addFiniteElementDescriptor(self, arg0):
        pass
    
    def addLoad(self, *args, **kwargs):
        """Overloaded function.
        
        1. addLoad(self: libaster.BaseDOFNumbering, arg0: MechanicalLoad<ConstantFieldOnCells<double> >) -> None
        
        2. addLoad(self: libaster.BaseDOFNumbering, arg0: MechanicalLoad<ConstantFieldOnCells<double> >, arg1: Function) -> None
        
        3. addLoad(self: libaster.BaseDOFNumbering, arg0: MechanicalLoad<ConstantFieldOnCells<double> >, arg1: Formula) -> None
        
        4. addLoad(self: libaster.BaseDOFNumbering, arg0: MechanicalLoad<ConstantFieldOnCells<double> >, arg1: Function2D) -> None
        
        5. addLoad(self: libaster.BaseDOFNumbering, arg0: MechanicalLoad<ConstantFieldOnCells<JeveuxString<24> > >) -> None
        
        6. addLoad(self: libaster.BaseDOFNumbering, arg0: MechanicalLoad<ConstantFieldOnCells<JeveuxString<24> > >, arg1: Function) -> None
        
        7. addLoad(self: libaster.BaseDOFNumbering, arg0: MechanicalLoad<ConstantFieldOnCells<JeveuxString<24> > >, arg1: Formula) -> None
        
        8. addLoad(self: libaster.BaseDOFNumbering, arg0: MechanicalLoad<ConstantFieldOnCells<JeveuxString<24> > >, arg1: Function2D) -> None
        
        9. addLoad(self: libaster.BaseDOFNumbering, arg0: MechanicalLoad<ConstantFieldOnCells<std::complex<double> > >) -> None
        
        10. addLoad(self: libaster.BaseDOFNumbering, arg0: MechanicalLoad<ConstantFieldOnCells<std::complex<double> > >, arg1: Function) -> None
        
        11. addLoad(self: libaster.BaseDOFNumbering, arg0: MechanicalLoad<ConstantFieldOnCells<std::complex<double> > >, arg1: Formula) -> None
        
        12. addLoad(self: libaster.BaseDOFNumbering, arg0: MechanicalLoad<ConstantFieldOnCells<std::complex<double> > >, arg1: Function2D) -> None
        
        13. addLoad(self: libaster.BaseDOFNumbering, arg0: ThermalLoad<ConstantFieldOnCells<double> >) -> None
        
        14. addLoad(self: libaster.BaseDOFNumbering, arg0: ThermalLoad<ConstantFieldOnCells<double> >, arg1: Function) -> None
        
        15. addLoad(self: libaster.BaseDOFNumbering, arg0: ThermalLoad<ConstantFieldOnCells<double> >, arg1: Formula) -> None
        
        16. addLoad(self: libaster.BaseDOFNumbering, arg0: ThermalLoad<ConstantFieldOnCells<double> >, arg1: Function2D) -> None
        
        17. addLoad(self: libaster.BaseDOFNumbering, arg0: ThermalLoad<ConstantFieldOnCells<JeveuxString<24> > >) -> None
        
        18. addLoad(self: libaster.BaseDOFNumbering, arg0: ThermalLoad<ConstantFieldOnCells<JeveuxString<24> > >, arg1: Function) -> None
        
        19. addLoad(self: libaster.BaseDOFNumbering, arg0: ThermalLoad<ConstantFieldOnCells<JeveuxString<24> > >, arg1: Formula) -> None
        
        20. addLoad(self: libaster.BaseDOFNumbering, arg0: ThermalLoad<ConstantFieldOnCells<JeveuxString<24> > >, arg1: Function2D) -> None
        
        21. addLoad(self: libaster.BaseDOFNumbering, arg0: AcousticLoad<ConstantFieldOnCells<std::complex<double> > >) -> None
        
        22. addLoad(self: libaster.BaseDOFNumbering, arg0: AcousticLoad<ConstantFieldOnCells<std::complex<double> > >, arg1: Function) -> None
        
        23. addLoad(self: libaster.BaseDOFNumbering, arg0: AcousticLoad<ConstantFieldOnCells<std::complex<double> > >, arg1: Formula) -> None
        
        24. addLoad(self: libaster.BaseDOFNumbering, arg0: AcousticLoad<ConstantFieldOnCells<std::complex<double> > >, arg1: Function2D) -> None
        
        25. addLoad(self: libaster.BaseDOFNumbering, arg0: ParallelMechanicalLoad<ConstantFieldOnCells<double> >) -> None
        
        26. addLoad(self: libaster.BaseDOFNumbering, arg0: ParallelMechanicalLoad<ConstantFieldOnCells<double> >, arg1: Function) -> None
        
        27. addLoad(self: libaster.BaseDOFNumbering, arg0: ParallelMechanicalLoad<ConstantFieldOnCells<double> >, arg1: Formula) -> None
        
        28. addLoad(self: libaster.BaseDOFNumbering, arg0: ParallelMechanicalLoad<ConstantFieldOnCells<double> >, arg1: Function2D) -> None
        
        29. addLoad(self: libaster.BaseDOFNumbering, arg0: ParallelMechanicalLoad<ConstantFieldOnCells<JeveuxString<24> > >) -> None
        
        30. addLoad(self: libaster.BaseDOFNumbering, arg0: ParallelMechanicalLoad<ConstantFieldOnCells<JeveuxString<24> > >, arg1: Function) -> None
        
        31. addLoad(self: libaster.BaseDOFNumbering, arg0: ParallelMechanicalLoad<ConstantFieldOnCells<JeveuxString<24> > >, arg1: Formula) -> None
        
        32. addLoad(self: libaster.BaseDOFNumbering, arg0: ParallelMechanicalLoad<ConstantFieldOnCells<JeveuxString<24> > >, arg1: Function2D) -> None
        
        33. addLoad(self: libaster.BaseDOFNumbering, arg0: ParallelThermalLoad<ConstantFieldOnCells<double> >) -> None
        
        34. addLoad(self: libaster.BaseDOFNumbering, arg0: ParallelThermalLoad<ConstantFieldOnCells<double> >, arg1: Function) -> None
        
        35. addLoad(self: libaster.BaseDOFNumbering, arg0: ParallelThermalLoad<ConstantFieldOnCells<double> >, arg1: Formula) -> None
        
        36. addLoad(self: libaster.BaseDOFNumbering, arg0: ParallelThermalLoad<ConstantFieldOnCells<double> >, arg1: Function2D) -> None
        
        37. addLoad(self: libaster.BaseDOFNumbering, arg0: ParallelThermalLoad<ConstantFieldOnCells<JeveuxString<24> > >) -> None
        
        38. addLoad(self: libaster.BaseDOFNumbering, arg0: ParallelThermalLoad<ConstantFieldOnCells<JeveuxString<24> > >, arg1: Function) -> None
        
        39. addLoad(self: libaster.BaseDOFNumbering, arg0: ParallelThermalLoad<ConstantFieldOnCells<JeveuxString<24> > >, arg1: Formula) -> None
        
        40. addLoad(self: libaster.BaseDOFNumbering, arg0: ParallelThermalLoad<ConstantFieldOnCells<JeveuxString<24> > >, arg1: Function2D) -> None
        """
    
    def computeNumbering(self):
        pass
    
    def computeRenumbering(self):
        pass
    
    def getDescription(self):
        pass
    
    def getDirichletBCDOFs(self):
        """Return a vector which describes DOFs that are imposed by Dirichlet BC.
        
        The vector has a size equals to the number of DOFs. For each dof, the value is equal to one
        if Dirichel BC is imposed to this DOF else zero
        
        Be carefull all Dirichlet BC have to be added before to call this function.
        
        Returns:
            tuple(int): a list with dirichlet imposition.
        """
    
    def getFiniteElementDescriptors(self):
        pass
    
    def getListOfLoads(self):
        pass
    
    def getMesh(self):
        """Return the mesh
        
        Returns:
            MeshPtr: a pointer to the mesh
        """
    
    def getModel(self):
        """Return the model
        
        Returns:
            ModelPtr: a pointer to the model
        """
    
    def getPhysicalQuantity(self):
        """Returns the name of the physical quantity that is numbered.
        
        Returns:
            str: physical quantity name.
        """
    
    def hasDirichletBC(self):
        """The list of loads used to build numbering contains Dirichlet BC.
        
        Returns:
            bool: *True* if Dirichlet BC are present, *False* otherwise.
        """
    
    def isParallel(self):
        """The numbering is distributed across MPI processes for High Performance Computing.
        
        Returns:
            bool: *True* if used, *False* otherwise.
        """
    
    def setDescription(self, arg0):
        pass
    
    def setElementaryMatrix(self, *args, **kwargs):
        """Overloaded function.
        
        1. setElementaryMatrix(self: libaster.BaseDOFNumbering, arg0: ElementaryMatrix<double, (PhysicalQuantityEnum)4>) -> None
        
        2. setElementaryMatrix(self: libaster.BaseDOFNumbering, arg0: ElementaryMatrix<std::complex<double>, (PhysicalQuantityEnum)4>) -> None
        
        3. setElementaryMatrix(self: libaster.BaseDOFNumbering, arg0: ElementaryMatrix<double, (PhysicalQuantityEnum)6>) -> None
        
        4. setElementaryMatrix(self: libaster.BaseDOFNumbering, arg0: ElementaryMatrix<std::complex<double>, (PhysicalQuantityEnum)5>) -> None
        """
    
    def setEmpty(self, arg0):
        pass
    
    def setModel(self, arg0):
        pass

# class DOFNumbering in libaster

class DOFNumbering(BaseDOFNumbering):
    pass
    
    # Method resolution order:
    #     DOFNumbering
    #     BaseDOFNumbering
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __init__(self, *args, **kwargs):
        """Overloaded function.
        
        1. __init__(self: libaster.DOFNumbering) -> None
        
        2. __init__(self: libaster.DOFNumbering, arg0: str) -> None
        
        3. __init__(self: libaster.DOFNumbering, arg0: str, arg1: Model, arg2: ListOfLoads, arg3: libaster.FieldOnNodesDescription) -> None
        """
    
    def getComponentAssociatedToRow(self, row, local= False):
        """Returns the component name associated to a dof index.
        
        - If the row is associated to a physical DOF, the name of the component is returned.
        
        - If the row is associated to a Lagrange multiplier DOF for a Dirichlet boundary
          condition, the name of the component which is constrained by the multiplier is
          returned, precedeed by 'LAGR:', e.g. 'LAGR:DX'.
        
        - If the row is associated to a Lagrange multiplier DOF for a multipoint-constraint
          (MPC) implying several DOF, 'LAGR:MPC' is returned (since no component can be
          identified).
        
        Arguments:
            row (int): Index of the dof.
            local (bool, optional): not used (default: false).
        
        Returns:
            str: component name.
        """
    
    def getComponents(self):
        """Returns all the component names assigned in the numbering.
        
        Returns:
            str: component names.
        """
    
    def getComponentsAssociatedToNode(self, node, local= False):
        """Returns the components name associated to a node index.
        
        Arguments:
            node (int): Index of the node.
            local (bool, optional): not used (default: false).
        
        Returns:
            str: component names.
        """
    
    def getNodeAssociatedToRow(self, row, local= False):
        """Returns the node index associated to a dof index.
        
        Arguments:
            row (int): Index of the dof.
            local (bool, optional): not used (default: false).
        
        Returns:
            int: index of the dof.
        """
    
    def getNumberOfDofs(self, local= False):
        """Returns the number of DOFs.
        
        Arguments:
            local (bool, optional): not used (default: false).
        
        Returns:
            int: number of DOFs.
        """
    
    def getRowAssociatedToNodeComponent(self, node, component, local= False):
        """Returns the index of the dof associated to a node.
        
        Arguments:
            node (int): Index of the node.
            component (str): name of the component
            local (bool, optional): not used (default: false).
        
        Returns:
            int: index of the dof.
        """
    
    def getRowsAssociatedToLagrangeMultipliers(self, local= False):
        """Returns the indexes of the Lagrange multipliers dof.
        
        Arguments:
            local (bool, optional): not used (default: false).
        
        Returns:
            [int]: indexes of the Lagrange multipliers dof.
        """
    
    def getRowsAssociatedToPhysicalDofs(self, local= False):
        """Returns the indexes of the physical dof.
        
        Arguments:
            local (bool, optional): not used (default: false).
        
        Returns:
            [int]: indexes of the physical dof.
        """
    
    def isRowAssociatedToPhysical(self, row, local= False):
        """If the row is associated to a physical DOF, return True
        
        If the row is associated to a Lagrange multiplier DOF for a Dirichlet boundary
          condition, return False
        
        Arguments:
            row (int): Index of the dof.
            local (bool, optional): not used (default: false).
        
        Returns:
            int: index of the dof.
        """
    
    def useLagrangeMultipliers(self):
        """Lagrange multipliers are used for BC or MPC.
        
        Returns:
            bool: *True* if used, *False* otherwise.
        """
    
    def useSingleLagrangeMultipliers(self):
        """Single Lagrange multipliers are used for BC or MPC.
        
        Returns:
            bool: *True* if used, *False* otherwise.
        """

# class ElementaryCharacteristics in libaster

class ElementaryCharacteristics(DataStructure):
    pass
    
    # Method resolution order:
    #     ElementaryCharacteristics
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __init__(self, *args, **kwargs):
        """Overloaded function.
        
        1. __init__(self: libaster.ElementaryCharacteristics, arg0: Model) -> None
        
        2. __init__(self: libaster.ElementaryCharacteristics, arg0: str, arg1: Model) -> None
        """
    
    def getMesh(self):
        pass
    
    def getModel(self):
        pass

# class FiniteElementDescriptor in libaster

class FiniteElementDescriptor(DataStructure):
    pass
    
    # Method resolution order:
    #     FiniteElementDescriptor
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __init__(self, *args, **kwargs):
        """Overloaded function.
        
        1. __init__(self: libaster.FiniteElementDescriptor, arg0: libaster.BaseMesh) -> None
        
        2. __init__(self: libaster.FiniteElementDescriptor, arg0: str, arg1: libaster.BaseMesh) -> None
        
        3. __init__(self: libaster.FiniteElementDescriptor, arg0: libaster.FiniteElementDescriptor, arg1: List[str]) -> None
        
        4. __init__(self: libaster.FiniteElementDescriptor, arg0: Model, arg1: List[str]) -> None
        """
    
    def getListOfGroupOfElements(self):
        pass
    
    def getMesh(self):
        pass
    
    def getModel(self):
        pass
    
    def getPhysics(self):
        pass
    
    def getVirtualCellsDescriptor(self):
        pass
    
    def setModel(self, arg0):
        pass
    
    def transferDofDescriptorFrom(self, arg0):
        pass

# class FiberGeometry in libaster

class FiberGeometry(DataStructure):
    pass
    
    # Method resolution order:
    #     FiberGeometry
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __init__(self, *args, **kwargs):
        """Overloaded function.
        
        1. __init__(self: libaster.FiberGeometry) -> None
        
        2. __init__(self: libaster.FiberGeometry, arg0: str) -> None
        """

# class DataField in libaster

class DataField(DataStructure):
    pass
    
    # Method resolution order:
    #     DataField
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __init__(self, *args, **kwargs):
        """Overloaded function.
        
        1. __init__(self: libaster.DataField) -> None
        
        2. __init__(self: libaster.DataField) -> None
        
        3. __init__(self: libaster.DataField, arg0: str) -> None
        
        4. __init__(self: libaster.DataField, arg0: str, arg1: str) -> None
        """
    
    def getFieldType(self):
        """Get field type between "ELEM", "ELGA", "ELNO", "NOEU", "CART".
        
        Returns:
            str: field type
        """

# class FieldOnCellsReal in libaster

class FieldOnCellsReal(DataField):
    pass
    
    # Method resolution order:
    #     FieldOnCellsReal
    #     DataField
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __add__(self, arg0):
        pass
    
    def __getitem__(self, arg0):
        pass
    
    def __iadd__(self, arg0):
        pass
    
    def __init__(self, *args, **kwargs):
        """Overloaded function.
        
        1. __init__(self: libaster.FieldOnCellsReal) -> None
        
        2. __init__(self: libaster.FieldOnCellsReal, arg0: str) -> None
        
        3. __init__(self: libaster.FieldOnCellsReal, arg0: Model, arg1: BehaviourProperty, arg2: str) -> None
        
        4. __init__(self: libaster.FieldOnCellsReal, arg0: Model, arg1: BehaviourProperty, arg2: str, arg3: libaster.ElementaryCharacteristics) -> None
        
        5. __init__(self: libaster.FieldOnCellsReal, arg0: libaster.FieldOnCellsReal) -> None
        """
    
    def __isub__(self, arg0):
        pass
    
    def __len__(self):
        pass
    
    def __mul__(self, arg0):
        pass
    
    def __neg__(self):
        pass
    
    def __rmul__(self, arg0):
        pass
    
    def __setitem__(self, arg0, arg1):
        pass
    
    def __sub__(self, arg0):
        pass
    
    def build(self, feds= []):
        pass
    
    def dot(self, other):
        """Return the dot product of two fields
        
        Arguments:
            field (FieldOnCells): other field
        
        Returns:
            float: dot product
        """
    
    def duplicate(self):
        """Return a duplicated FieldOnCellsReal as a copy
        
        Returns:
            FieldOnCellsReal
        """
    
    def exportToSimpleFieldOnCells(self):
        pass
    
    def getComponents(self):
        """Get list of components
        
        Returns:
            list[str]: list of components
        """
    
    def getDescription(self, *args, **kwargs):
        """Overloaded function.
        
        1. getDescription(self: libaster.FieldOnCellsReal) -> libaster.FiniteElementDescriptor
        
        
                    Return the descriptor associated with the FieldOnCellsReal object
        
                    Returns:
                        FiniteElementDescriptor: FiniteElementDescriptor Object
                    
        
        2. getDescription(self: libaster.FieldOnCellsReal) -> libaster.FiniteElementDescriptor
        """
    
    def getMesh(self):
        """Return the Mesh associated with the FieldOnCellsReal object
        
        Returns:
            BaseMesh: Mesh object
        """
    
    def getModel(self):
        """Return the Model associated with the FieldOnCellsReal object
        
        Returns:
            Model: Model object
        """
    
    def getNumberOfComponents(self):
        """Get number of components
        
        Returns:
            int: number of components
        """
    
    def getPhysicalQuantity(self):
        """Get physical quantity
        
        Returns:
            str: physical quantity
        """
    
    def getValues(self):
        """Return a list of values as (x1, y1, z1, x2, y2, z2...)
        
        Returns:
            list[float]: List of values.
        """
    
    def norm(self, arg0):
        """Return the euclidean norm of the field
        
        Arguments:
            normType (str): "NORM_1", "NORM_2", "NORM_INFINITY"
        
        Returns:
            float: euclidean norm
        """
    
    def printMedFile(self, filename, local= True):
        """Print the field in MED format.
        
        Arguments:
            filename (str): Path to the file to be printed.
            local (bool): Print local values only (relevant for ParallelMesh only,
                default: *True*)
        
        Returns:
            bool: *True* if succeeds, *False* otherwise.
        """
    
    def setDescription(self, arg0):
        pass
    
    def setValues(self, value):
        """Set values of the field
        
        Arguments:
            float: value to set
        """
    
    def size(self):
        """Return the size of the field
        
        Returns:
            int: number of element in the field
        """
    
    def transform(self, func):
        """Apply Function to each value of _ValuesList of the FieldOnCells object.
        
        Arguments:
            func (Function): Callable Python object
        
        Returns:
            FieldOnCellsReal: New FieldOnCells object with the trasformed values
        """

# class FieldOnCellsComplex in libaster

class FieldOnCellsComplex(DataField):
    pass
    
    # Method resolution order:
    #     FieldOnCellsComplex
    #     DataField
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __add__(self, arg0):
        pass
    
    def __getitem__(self, arg0):
        pass
    
    def __iadd__(self, arg0):
        pass
    
    def __init__(self, *args, **kwargs):
        """Overloaded function.
        
        1. __init__(self: libaster.FieldOnCellsComplex) -> None
        
        2. __init__(self: libaster.FieldOnCellsComplex, arg0: str) -> None
        
        3. __init__(self: libaster.FieldOnCellsComplex, arg0: libaster.FieldOnCellsComplex) -> None
        """
    
    def __isub__(self, arg0):
        pass
    
    def __len__(self):
        pass
    
    def __mul__(self, arg0):
        pass
    
    def __rmul__(self, arg0):
        pass
    
    def __setitem__(self, arg0, arg1):
        pass
    
    def __sub__(self, arg0):
        pass
    
    def build(self, feds= []):
        pass
    
    def duplicate(self):
        pass
    
    def getDescription(self):
        pass
    
    def getMesh(self):
        """Return the Mesh associated with the FieldOnCellsReal object
        
        Returns:
            BaseMesh: Mesh object
        """
    
    def getModel(self):
        """Return the Model associated with the FieldOnCellsReal object
        
        Returns:
            Model: Model object
        """
    
    def getValues(self):
        """Return a list of values as (x1, y1, z1, x2, y2, z2...)
        
        Returns:
            list[complex]: List of values.
        """
    
    def printMedFile(self, filename, local= True):
        """Print the field in MED format.
        
        Arguments:
            filename (str): Path to the file to be printed.
        
        Returns:
            bool: *True* if succeeds, *False* otherwise.
        """
    
    def setDescription(self, arg0):
        pass
    
    def size(self):
        """Return the size of the field
        
        Returns:
            int: number of element in the field
        """
    
    def transform(self, func):
        """Apply Function to each value of _ValuesList of the FieldOnCells object.
        
        Arguments:
            func (Function): Callable Python object
        
        Returns:
            bool: New FieldOnCells object with the trasformed values
        """

# class FieldOnCellsLong in libaster

class FieldOnCellsLong(DataField):
    pass
    
    # Method resolution order:
    #     FieldOnCellsLong
    #     DataField
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __add__(self, arg0):
        pass
    
    def __getitem__(self, arg0):
        pass
    
    def __iadd__(self, arg0):
        pass
    
    def __init__(self, *args, **kwargs):
        """Overloaded function.
        
        1. __init__(self: libaster.FieldOnCellsLong) -> None
        
        2. __init__(self: libaster.FieldOnCellsLong, arg0: str) -> None
        
        3. __init__(self: libaster.FieldOnCellsLong, arg0: libaster.FieldOnCellsLong) -> None
        """
    
    def __isub__(self, arg0):
        pass
    
    def __len__(self):
        pass
    
    def __mul__(self, arg0):
        pass
    
    def __rmul__(self, arg0):
        pass
    
    def __setitem__(self, arg0, arg1):
        pass
    
    def __sub__(self, arg0):
        pass
    
    def build(self, feds= []):
        pass
    
    def duplicate(self):
        pass
    
    def getDescription(self):
        pass
    
    def getMesh(self):
        """Return the Mesh associated with the FieldOnCellsReal object
        
        Returns:
            BaseMesh: Mesh object
        """
    
    def getModel(self):
        """Return the Model associated with the FieldOnCellsReal object
        
        Returns:
            Model: Model object
        """
    
    def getValues(self):
        """Return a list of values as (x1, y1, z1, x2, y2, z2...)
        
        Returns:
            list[int]: List of values.
        """
    
    def printMedFile(self, filename, local= True):
        """Print the field in MED format.
        
        Arguments:
            filename (str): Path to the file to be printed.
        
        Returns:
            bool: *True* if succeeds, *False* otherwise.
        """
    
    def setDescription(self, arg0):
        pass
    
    def size(self):
        """Return the size of the field
        
        Returns:
            int: number of element in the field
        """

# class FieldOnCellsChar8 in libaster

class FieldOnCellsChar8(DataField):
    pass
    
    # Method resolution order:
    #     FieldOnCellsChar8
    #     DataField
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __init__(self, *args, **kwargs):
        """Overloaded function.
        
        1. __init__(self: libaster.FieldOnCellsChar8) -> None
        
        2. __init__(self: libaster.FieldOnCellsChar8, arg0: str) -> None
        
        3. __init__(self: libaster.FieldOnCellsChar8, arg0: libaster.FieldOnCellsChar8) -> None
        """
    
    def build(self, feds= []):
        pass
    
    def getDescription(self):
        """Return the description associated with the FieldOnCellsChar8 object
        
        Returns:
            FiniteElementDescriptor: FiniteElementDescriptor Object
        """
    
    def getMesh(self):
        """Return the Mesh associated with the FieldOnCellsChar8 object
        
        Returns:
            BaseMesh: Mesh object
        """
    
    def getModel(self):
        """Return the Model associated with the FieldOnCellsReal object
        
        Returns:
            Model: Model object
        """
    
    def setDescription(self, arg0):
        pass

# class FieldOnNodesReal in libaster

class FieldOnNodesReal(DataField):
    pass
    
    # Method resolution order:
    #     FieldOnNodesReal
    #     DataField
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __add__(self, arg0):
        pass
    
    def __getitem__(self, arg0):
        pass
    
    def __iadd__(self, arg0):
        pass
    
    def __imul__(self, arg0):
        pass
    
    def __init__(self, *args, **kwargs):
        """Overloaded function.
        
        1. __init__(self: libaster.FieldOnNodesReal) -> None
        
        2. __init__(self: libaster.FieldOnNodesReal, arg0: str) -> None
        
        3. __init__(self: libaster.FieldOnNodesReal, arg0: libaster.FieldOnNodesReal) -> None
        
        4. __init__(self: libaster.FieldOnNodesReal, arg0: Model) -> None
        
        5. __init__(self: libaster.FieldOnNodesReal, arg0: libaster.BaseDOFNumbering) -> None
        
        6. __init__(self: libaster.FieldOnNodesReal, arg0: MeshCoordinatesField) -> None
        """
    
    def __isub__(self, arg0):
        pass
    
    def __mul__(self, arg0):
        pass
    
    def __neg__(self):
        pass
    
    def __rmul__(self, arg0):
        pass
    
    def __setitem__(self, arg0, arg1):
        pass
    
    def __sub__(self, arg0):
        pass
    
    def build(self):
        pass
    
    def dot(self, other):
        """Return the dot product of two fields
        
        Arguments:
            other (FieldOnNodes): other field
        
        Returns:
            float: dot product
        """
    
    def duplicate(self):
        pass
    
    def exportToSimpleFieldOnNodes(self):
        pass
    
    def getComponents(self):
        """Get list of components
        
        Returns:
            list[str]: list of components
        """
    
    def getDescription(self):
        pass
    
    def getMesh(self, *args, **kwargs):
        """Overloaded function.
        
        1. getMesh(self: libaster.FieldOnNodesReal) -> libaster.BaseMesh
        
        2. getMesh(self: libaster.FieldOnNodesReal) -> libaster.BaseMesh
        """
    
    def getNodesAndComponentsNumberFromDOF(self):
        """Return a list of values such that for each DOF, it gives the node id and component id
        as [dof1=[node_1, comp_1], dof2=[node_1, comp_2], ....]
        
        Returns:
            list[[int, int]]: List of values (node, component) for each DOF.
        """
    
    def getNumberOfComponents(self):
        """Get number of components
        
        Returns:
            int: number of components
        """
    
    def getPhysicalQuantity(self):
        pass
    
    def getValues(self):
        """Return a list of values as (x1, y1, z1, x2, y2, z2...)
        
        Returns:
            list[float]: List of values.
        """
    
    def norm(self, normType, list_cmp= []):
        """Return the euclidean norm of the field
        
        Arguments:
            normType (str): "NORM_1", "NORM_2", "NORM_INFINITY"
            list_cmp (list[str]) : list of components used to compute norm (default: all)
        
        Returns:
            float: euclidean norm
        """
    
    def printMedFile(self, fileName, local= True):
        pass
    
    def scale(self, vect):
        """Scale in-place the field by a diagonal matrix stored as an array
        
        Arguments:
            vect (float): diagonal matrix stored as an array
        """
    
    def setDescription(self, arg0):
        pass
    
    def setMesh(self, arg0):
        pass
    
    def setValues(self, *args, **kwargs):
        """Overloaded function.
        
        1. setValues(self: libaster.FieldOnNodesReal, value: float) -> None
        
        
                    Set values of the field
        
                    Arguments:
                        value (float): value to set
                    
        
        2. setValues(self: libaster.FieldOnNodesReal, value: List[float]) -> None
        
        
                    Set values of the field
        
                    Arguments:
                        value (list[float]): list of values to set
        """
    
    def size(self):
        """Return the size of the field
        
        Returns:
            int: number of element in the field
        """
    
    def updateValuePointers(self):
        pass

# class FieldOnNodesComplex in libaster

class FieldOnNodesComplex(DataField):
    pass
    
    # Method resolution order:
    #     FieldOnNodesComplex
    #     DataField
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __getitem__(self, arg0):
        pass
    
    def __init__(self, *args, **kwargs):
        """Overloaded function.
        
        1. __init__(self: libaster.FieldOnNodesComplex) -> None
        
        2. __init__(self: libaster.FieldOnNodesComplex, arg0: str) -> None
        
        3. __init__(self: libaster.FieldOnNodesComplex, arg0: libaster.FieldOnNodesComplex) -> None
        
        4. __init__(self: libaster.FieldOnNodesComplex, arg0: Model) -> None
        
        5. __init__(self: libaster.FieldOnNodesComplex, arg0: libaster.BaseDOFNumbering) -> None
        """
    
    def __setitem__(self, arg0, arg1):
        pass
    
    def build(self):
        pass
    
    def dot(self, other):
        """Return the dot product of two complex fields
        
        Arguments:
            other (FieldOnNodes): other field
        
        Returns:
            complex: dot product
        """
    
    def exportToSimpleFieldOnNodes(self):
        pass
    
    def getComponents(self):
        """Get list of components
        
        Returns:
            list[str]: list of components
        """
    
    def getDescription(self):
        pass
    
    def getMesh(self, *args, **kwargs):
        """Overloaded function.
        
        1. getMesh(self: libaster.FieldOnNodesComplex) -> libaster.BaseMesh
        
        2. getMesh(self: libaster.FieldOnNodesComplex) -> libaster.BaseMesh
        """
    
    def getNumberOfComponents(self):
        """Get number of components
        
        Returns:
            int: number of components
        """
    
    def getPhysicalQuantity(self):
        pass
    
    def getValues(self):
        """Return a list of complex values as [x11, x21, ..., xm1, x12, x22, ..., xm2...]
        (m is the total number of componenets)
        
        Returns:
            list[complex]: List of values.
        """
    
    def norm(self, normType, list_cmp= []):
        """Return the euclidean norm of the field
        
        Arguments:
            normType (str): "NORM_1", "NORM_2", "NORM_INFINITY"
            list_cmp (list[str]) : list of components used to compute norm (default: all)
        
        Returns:
            float: euclidean norm
        """
    
    def printMedFile(self, arg0, arg1):
        pass
    
    def scale(self, vect):
        """Scale in-place the field by a diagonal matrix stored as an array
        
        Arguments:
            vect (float): diagonal matrix stored as an array
        """
    
    def setDescription(self, arg0):
        pass
    
    def setMesh(self, arg0):
        pass
    
    def setValues(self, *args, **kwargs):
        """Overloaded function.
        
        1. setValues(self: libaster.FieldOnNodesComplex, value: complex) -> None
        
        
                    Set values of the field
        
                    Arguments:
                        value (complex): value to set
                    
        
        2. setValues(self: libaster.FieldOnNodesComplex, value: List[complex]) -> None
        
        
                    Set values of the field
        
                    Arguments:
                        value (list[complex]): list of values to set
        """
    
    def updateValuePointers(self):
        pass

# class FieldOnNodesLong in libaster

class FieldOnNodesLong(DataField):
    pass
    
    # Method resolution order:
    #     FieldOnNodesLong
    #     DataField
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __init__(self, *args, **kwargs):
        """Overloaded function.
        
        1. __init__(self: libaster.FieldOnNodesLong) -> None
        
        2. __init__(self: libaster.FieldOnNodesLong, arg0: str) -> None
        
        3. __init__(self: libaster.FieldOnNodesLong, arg0: libaster.FieldOnNodesLong) -> None
        
        4. __init__(self: libaster.FieldOnNodesLong, arg0: libaster.BaseDOFNumbering) -> None
        """
    
    def getDescription(self):
        pass
    
    def getMesh(self):
        pass
    
    def setDescription(self, arg0):
        pass
    
    def setMesh(self, arg0):
        pass

# class FieldOnNodesChar8 in libaster

class FieldOnNodesChar8(DataField):
    pass
    
    # Method resolution order:
    #     FieldOnNodesChar8
    #     DataField
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __init__(self, *args, **kwargs):
        """Overloaded function.
        
        1. __init__(self: libaster.FieldOnNodesChar8) -> None
        
        2. __init__(self: libaster.FieldOnNodesChar8, arg0: str) -> None
        
        3. __init__(self: libaster.FieldOnNodesChar8, arg0: libaster.FieldOnNodesChar8) -> None
        
        4. __init__(self: libaster.FieldOnNodesChar8, arg0: libaster.BaseDOFNumbering) -> None
        """
    
    def getDescription(self):
        pass
    
    def getMesh(self):
        pass
    
    def setDescription(self, arg0):
        pass
    
    def setMesh(self, arg0):
        pass

# class ConstantFieldOnCellsReal in libaster

class ConstantFieldOnCellsReal(DataField):
    pass
    
    # Method resolution order:
    #     ConstantFieldOnCellsReal
    #     DataField
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __init__(self, *args, **kwargs):
        """Overloaded function.
        
        1. __init__(self: libaster.ConstantFieldOnCellsReal, arg0: libaster.BaseMesh) -> None
        
        2. __init__(self: libaster.ConstantFieldOnCellsReal, arg0: str, arg1: libaster.BaseMesh) -> None
        """
    
    def getMesh(self):
        pass

# class ConstantFieldOnCellsChar16 in libaster

class ConstantFieldOnCellsChar16(DataField):
    pass
    
    # Method resolution order:
    #     ConstantFieldOnCellsChar16
    #     DataField
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __init__(self, *args, **kwargs):
        """Overloaded function.
        
        1. __init__(self: libaster.ConstantFieldOnCellsChar16, arg0: libaster.BaseMesh) -> None
        
        2. __init__(self: libaster.ConstantFieldOnCellsChar16, arg0: str, arg1: libaster.BaseMesh) -> None
        """
    
    def getMesh(self):
        pass

# class ConstantFieldOnCellsLong in libaster

class ConstantFieldOnCellsLong(DataField):
    pass
    
    # Method resolution order:
    #     ConstantFieldOnCellsLong
    #     DataField
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __init__(self, *args, **kwargs):
        """Overloaded function.
        
        1. __init__(self: libaster.ConstantFieldOnCellsLong, arg0: libaster.BaseMesh) -> None
        
        2. __init__(self: libaster.ConstantFieldOnCellsLong, arg0: str, arg1: libaster.BaseMesh) -> None
        """
    
    def getMesh(self):
        pass

# class SimpleFieldOnCellsReal in libaster

class SimpleFieldOnCellsReal(DataStructure):
    pass
    
    # Method resolution order:
    #     SimpleFieldOnCellsReal
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __init__(self, *args, **kwargs):
        """Overloaded function.
        
        1. __init__(self: libaster.SimpleFieldOnCellsReal) -> None
        
        2. __init__(self: libaster.SimpleFieldOnCellsReal, arg0: str) -> None
        """
    
    def getCellsWithComponents(self):
        """Returns the list of cells where the field is defined.
        
        Returns:
            tuple (int): Indexes of cells where the field is defined.
        """
    
    def getFieldLocation(self):
        pass
    
    def getMaxNumberOfPoints(self):
        pass
    
    def getMaxNumberOfSubPoints(self):
        pass
    
    def getNameOfComponent(self, arg0):
        pass
    
    def getNameOfComponents(self):
        pass
    
    def getNumberOfCells(self):
        pass
    
    def getNumberOfComponents(self):
        pass
    
    def getNumberOfComponentsForSubpointsOfCell(self, arg0):
        pass
    
    def getNumberOfPointsOfCell(self, arg0):
        pass
    
    def getNumberOfSubPointsOfCell(self, arg0):
        pass
    
    def getPhysicalQuantity(self):
        pass
    
    def getValue(self, ima, icmp, ipt, ispt):
        """Returns the value of the `icmp` component of the field on the `ima` cell,
        at the `ipt` point, at the `ispt` sub-point.
        
        Args:
            ima  (int): Index of cells.
            icmp (int): Index of component.
            ipt  (int): Index of point.
            ispt (int): Index of sub-point.
        
        Returns:
            float: Value of field at *ima*, of *icmp*, at *ipt*, at *ispt*;
                     NaN if the position is not allocated.
        """
    
    def getValues(self, copy= False):
        """Returns two numpy arrays with shape ( number_of_cells_with_components, number_of_components )
        The first array contains the field values while the second one is a mask
        which is `True` if the corresponding value exists, `False` otherwise.
        
        Where the mask is `False` the corresponding value is set to zero.
        
        Args:
                copy (bool): If True copy the data, default: *False*
        
        Returns:
            ndarray (float): Field values.
            ndarray (bool): Mask for the field values.
        """
    
    def updateValuePointers(self):
        pass

# class SimpleFieldOnNodesReal in libaster

class SimpleFieldOnNodesReal(DataStructure):
    pass
    
    # Method resolution order:
    #     SimpleFieldOnNodesReal
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __init__(self, *args, **kwargs):
        """Overloaded function.
        
        1. __init__(self: libaster.SimpleFieldOnNodesReal) -> None
        
        2. __init__(self: libaster.SimpleFieldOnNodesReal, arg0: str) -> None
        """
    
    def getNameOfComponent(self, arg0):
        pass
    
    def getNameOfComponents(self):
        pass
    
    def getNumberOfComponents(self):
        pass
    
    def getNumberOfNodes(self):
        pass
    
    def getPhysicalQuantity(self):
        pass
    
    def getValue(self, ino, icmp):
        """Returns the value of the `icmp` component of the field on the `ino` node.
        
        Arguments:
                ino (int): Index of node.
                icmp (int): Index of component.
        
        Returns:
            (float): The field value. NaN is returned if the position is not allocated.
        """
    
    def getValues(self, copy= False):
        """Returns two numpy arrays with shape ( number_of_components, space_dimension )
        The first array contains the field values while the second one is a mask
        which is `True` if the corresponding value exists, `False` otherwise.
        
        Where the mask is `False` the corresponding value is set to zero.
        
        Args:
                copy (bool): If True copy the data, default: *False*
        
        Returns:
            ndarray (float): Field values.
            ndarray (bool): Mask for the field values.
        """
    
    def updateValuePointers(self):
        pass

# class SimpleFieldOnNodesComplex in libaster

class SimpleFieldOnNodesComplex(DataStructure):
    pass
    
    # Method resolution order:
    #     SimpleFieldOnNodesComplex
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __init__(self, *args, **kwargs):
        """Overloaded function.
        
        1. __init__(self: libaster.SimpleFieldOnNodesComplex) -> None
        
        2. __init__(self: libaster.SimpleFieldOnNodesComplex, arg0: str) -> None
        """
    
    def getNameOfComponent(self, arg0):
        pass
    
    def getNameOfComponents(self):
        pass
    
    def getNumberOfComponents(self):
        pass
    
    def getNumberOfNodes(self):
        pass
    
    def getPhysicalQuantity(self):
        pass
    
    def getValue(self, ino, icmp):
        """Returns the value of the `icmp` component of the field on the `ino` node.
        
        Arguments:
            ino (int): Index of node.
            icmp (int): Index of component.
        
        Returns:
            complex: The field value. NaN is returned if the position is not allocated.
        """
    
    def getValues(self, copy= False):
        """Returns two numpy arrays with shape ( number_of_components, space_dimension )
        The first array contains the field values while the second one is a mask
        which is `True` if the corresponding value exists, `False` otherwise.
        
        Where the mask is `False` the corresponding value is set to zero.
        
        Args:
                copy (bool): If True copy the data, default: *False*
        
        Returns:
            ndarray (complex): Field values.
            ndarray (bool): Mask for the field values.
        """
    
    def updateValuePointers(self):
        pass

# class Table in libaster

class Table(DataStructure):
    pass
    
    # Method resolution order:
    #     Table
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __init__(self, *args, **kwargs):
        """Overloaded function.
        
        1. __init__(self: libaster.Table) -> None
        
        2. __init__(self: libaster.Table, arg0: str) -> None
        """

# class TableOfFunctions in libaster

class TableOfFunctions(Table):
    pass
    
    # Method resolution order:
    #     TableOfFunctions
    #     Table
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __init__(self, *args, **kwargs):
        """Overloaded function.
        
        1. __init__(self: libaster.TableOfFunctions) -> None
        
        2. __init__(self: libaster.TableOfFunctions, arg0: str) -> None
        """
    
    def addFunction(self, arg0):
        pass
    
    def getFunction(self, arg0):
        pass
    
    def getNumberOfFunctions(self):
        pass

# class TableContainer in libaster

class TableContainer(Table):
    pass
    
    # Method resolution order:
    #     TableContainer
    #     Table
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __init__(self, *args, **kwargs):
        """Overloaded function.
        
        1. __init__(self: libaster.TableContainer) -> None
        
        2. __init__(self: libaster.TableContainer, arg0: str) -> None
        """
    
    def addObject(self, *args, **kwargs):
        """Overloaded function.
        
        1. addObject(self: libaster.TableContainer, arg0: str, arg1: ElementaryMatrix<double, (PhysicalQuantityEnum)4>) -> None
        
        2. addObject(self: libaster.TableContainer, arg0: str, arg1: ElementaryMatrix<double, (PhysicalQuantityEnum)6>) -> None
        
        3. addObject(self: libaster.TableContainer, arg0: str, arg1: ElementaryVector<double, (PhysicalQuantityEnum)4>) -> None
        
        4. addObject(self: libaster.TableContainer, arg0: str, arg1: ElementaryVector<double, (PhysicalQuantityEnum)6>) -> None
        
        5. addObject(self: libaster.TableContainer, arg0: str, arg1: libaster.FieldOnCellsReal) -> None
        
        6. addObject(self: libaster.TableContainer, arg0: str, arg1: libaster.FieldOnNodesReal) -> None
        
        7. addObject(self: libaster.TableContainer, arg0: str, arg1: FunctionComplex) -> None
        
        8. addObject(self: libaster.TableContainer, arg0: str, arg1: GeneralizedAssemblyMatrix<double>) -> None
        
        9. addObject(self: libaster.TableContainer, arg0: str, arg1: libaster.DataField) -> None
        
        10. addObject(self: libaster.TableContainer, arg0: str, arg1: ModeResult) -> None
        
        11. addObject(self: libaster.TableContainer, arg0: str, arg1: libaster.ConstantFieldOnCellsReal) -> None
        
        12. addObject(self: libaster.TableContainer, arg0: str, arg1: Function2D) -> None
        
        13. addObject(self: libaster.TableContainer, arg0: str, arg1: libaster.Table) -> None
        
        14. addObject(self: libaster.TableContainer, name: str, object: Function) -> None
        
        
                    Store a *DataStructure* in the table.
        
                    Arguments:
                        name (str): String to identify the object in the table.
                        object (misc): Object to be inserted, can be a Function, Mesh, Fields...
        """
    
    def getConstantFieldOnCellsReal(self, arg0):
        pass
    
    def getDataField(self, arg0):
        pass
    
    def getElementaryMatrixDisplacementReal(self, arg0):
        pass
    
    def getElementaryMatrixTemperatureReal(self, arg0):
        pass
    
    def getElementaryVectorDisplacementReal(self, arg0):
        pass
    
    def getElementaryVectorTemperatureReal(self, arg0):
        pass
    
    def getFieldOnCellsReal(self, arg0):
        pass
    
    def getFieldOnNodesReal(self, arg0):
        pass
    
    def getFunction(self, arg0):
        pass
    
    def getFunction2D(self, arg0):
        pass
    
    def getFunctionComplex(self, arg0):
        pass
    
    def getGeneralizedAssemblyMatrix(self, arg0):
        pass
    
    def getModeResult(self, arg0):
        pass
    
    def getTable(self, arg0):
        pass

# class TimeStepper in libaster

class TimeStepper(DataStructure):
    pass
    
    # Method resolution order:
    #     TimeStepper
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __init__(self, *args, **kwargs):
        """Overloaded function.
        
        1. __init__(self: libaster.TimeStepper) -> None
        
        2. __init__(self: libaster.TimeStepper, arg0: str) -> None
        """
    
    def getValues(self):
        pass
    
    def setValues(self, arg0):
        pass

# class GeneralizedDOFNumbering in libaster

class GeneralizedDOFNumbering(DataStructure):
    pass
    
    # Method resolution order:
    #     GeneralizedDOFNumbering
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __init__(self, *args, **kwargs):
        """Overloaded function.
        
        1. __init__(self: libaster.GeneralizedDOFNumbering) -> None
        
        2. __init__(self: libaster.GeneralizedDOFNumbering, arg0: str) -> None
        """
    
    def getGeneralizedModel(self):
        pass
    
    def getModalBasis(self):
        pass
    
    def setGeneralizedModel(self, arg0):
        pass
    
    def setModalBasis(self, *args, **kwargs):
        """Overloaded function.
        
        1. setModalBasis(self: libaster.GeneralizedDOFNumbering, arg0: ModeResult) -> bool
        
        2. setModalBasis(self: libaster.GeneralizedDOFNumbering, arg0: GeneralizedModeResult) -> bool
        """

# class FluidStructureInteraction in libaster

class FluidStructureInteraction(DataStructure):
    pass
    
    # Method resolution order:
    #     FluidStructureInteraction
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __init__(self, *args, **kwargs):
        """Overloaded function.
        
        1. __init__(self: libaster.FluidStructureInteraction) -> None
        
        2. __init__(self: libaster.FluidStructureInteraction, arg0: str) -> None
        """

# class TurbulentSpectrum in libaster

class TurbulentSpectrum(DataStructure):
    pass
    
    # Method resolution order:
    #     TurbulentSpectrum
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __init__(self, *args, **kwargs):
        """Overloaded function.
        
        1. __init__(self: libaster.TurbulentSpectrum) -> None
        
        2. __init__(self: libaster.TurbulentSpectrum, arg0: str) -> None
        """

# class GenericFunction in libaster

class GenericFunction(DataStructure):
    pass
    
    # Method resolution order:
    #     GenericFunction
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __init__(self,  *args, **kwargs):
        """Initialize self.  See help(type(self)) for accurate signature.
        """
    
    def getProperties(self):
        pass
    
    def setExtrapolation(self, arg0):
        pass

# class ListOfLoads in libaster

class ListOfLoads(DataStructure):
    pass
    
    # Method resolution order:
    #     ListOfLoads
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __init__(self, *args, **kwargs):
        """Overloaded function.
        
        1. __init__(self: libaster.ListOfLoads) -> None
        
        2. __init__(self: libaster.ListOfLoads, arg0: str) -> None
        
        3. __init__(self: libaster.ListOfLoads, arg0: Model) -> None
        
        4. __init__(self: libaster.ListOfLoads, arg0: str, arg1: Model) -> None
        """
    
    def addContactLoadDescriptor(self, FED_Slave, FED_Pair):
        """Add contact load descriptor.
        
        Arguments:
            FED_Slave (FiniteElementDescriptor): Finite Element Descriptor defining
                slave cells (in DEFI_CONTACT)
            FED_Pair (FiniteElementDescriptor): Finite Element Descriptor defining
                list of contact pair
        """
    
    def getContactLoadDescriptor(self):
        """Get contact load descriptors.
        
        Returns:
            (FiniteElementDescriptor): Finite Element Descriptor defining
                slave cells (in DEFI_CONTACT)
            (FiniteElementDescriptor): Finite Element Descriptor defining
                list of contact pair
        """
    
    def getDirichletBCs(self):
        """Return list of DirichletBCs
        
        Returns:
            ListDiriBC: a list of DirichletBC
        """
    
    def getMechanicalLoadsFunction(self):
        """Return list of Function mechanical loads
        
        Returns:
            ListMecaLoadFunction: a list of Function mechanical loads
        """
    
    def getMechanicalLoadsReal(self):
        """Return list of real mechanical loads
        
        Returns:
            ListMecaLoadReal: a list of real mechanical loads
        """
    
    def getModel(self):
        """Return the model used
        
        Returns:
            Model: model used
        """
    
    def getParallelMechanicalLoadsFunction(self):
        """Return list of function parallel mechanical loads
        
        Returns:
            ListParaMecaLoadFunction: a list of function parallel mechanical loads
        """
    
    def getParallelMechanicalLoadsReal(self):
        """Return list of real parallel mechanical loads
        
        Returns:
            ListParaMecaLoadReal: a list of real parallel mechanical loads
        """
    
    def hasDirichletBC(self):
        """Dirichlet BCs have been added or not ?
        
        Returns:
            bool: True if Dirichlet BCs have been added
        """
    
    def hasExternalLoad(self):
        """External load (= not Dirichlet BCs) have been added or not ?
        
        Returns:
            bool: True if External load have been added
        """
    
    def isEmpty(self):
        """The list of loads is empty or not.
        
        Returns:
            bool: True if empty
        """

# class BaseFunction in libaster

class BaseFunction(GenericFunction):
    pass
    
    # Method resolution order:
    #     BaseFunction
    #     GenericFunction
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __init__(self,  *args, **kwargs):
        """Initialize self.  See help(type(self)) for accurate signature.
        """
    
    def getValues(self):
        """Return a list of the values of the function as (x1, x2, ..., y1, y2, ...)
        
        Returns:
            list[float]: List of values (size = 2 * *size()*).
        """
    
    def setInterpolation(self, arg0):
        pass
    
    def setParameterName(self, arg0):
        pass
    
    def setResultName(self, arg0):
        pass
    
    def setValues(self, arg0, arg1):
        pass

# class Function in libaster

class Function(BaseFunction):
    pass
    
    # Method resolution order:
    #     Function
    #     BaseFunction
    #     GenericFunction
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __init__(self, *args, **kwargs):
        """Overloaded function.
        
        1. __init__(self: libaster.Function) -> None
        
        2. __init__(self: libaster.Function, arg0: str) -> None
        """
    
    def setAsConstant(self):
        pass
    
    def setValues(self, arg0, arg1):
        pass
    
    def size(self):
        """Return the number of points of the function.
        
        Returns:
            int: Number of points.
        """

# class FunctionComplex in libaster

class FunctionComplex(BaseFunction):
    pass
    
    # Method resolution order:
    #     FunctionComplex
    #     BaseFunction
    #     GenericFunction
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __init__(self, *args, **kwargs):
        """Overloaded function.
        
        1. __init__(self: libaster.FunctionComplex) -> None
        
        2. __init__(self: libaster.FunctionComplex, arg0: str) -> None
        """
    
    def setValues(self, *args, **kwargs):
        """Overloaded function.
        
        1. setValues(self: libaster.FunctionComplex, arg0: List[float], arg1: List[float]) -> None
        
        2. setValues(self: libaster.FunctionComplex, arg0: List[float], arg1: List[complex]) -> None
        """
    
    def size(self):
        """Return the number of points of the function.
        
        Returns:
            int: Number of points.
        """

# class Formula in libaster

class Formula(GenericFunction):
    pass
    
    # Method resolution order:
    #     Formula
    #     GenericFunction
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __init__(self, *args, **kwargs):
        """Overloaded function.
        
        1. __init__(self: libaster.Formula) -> None
        
        2. __init__(self: libaster.Formula, arg0: str) -> None
        """
    
    def evaluate(self, *val):
        """Evaluate the formula with the given variables values.
        
        Arguments:
            val (list[float]): List of the values of the variables.
        
        Returns:
            float/complex: Value of the formula for these values.
        """
    
    def getContext(self):
        """Return the context used to evaluate the formula.
        
        Returns:
            dict: Context used for evaluation.
        """
    
    def getExpression(self):
        """Return expression of the formula.
        
        Returns:
            str: Expression of the formula.
        """
    
    def getProperties(self):
        """Return the properties of the formula (for compatibility with function objects).
        
        Returns:
            list[str]: List of 6 strings as function objects contain.
        """
    
    def getVariables(self):
        """Return the variables names.
        
        Returns:
            list[str]: List of the names of the variables.
        """
    
    def setComplex(self):
        """Set the type of the formula as complex.
        """
    
    def setContext(self, context):
        """Define the context holding objects required to evaluate the expression.
        
        Arguments:
            context (dict): Context for the evaluation.
        """
    
    def setExpression(self, expression):
        """Define the expression of the formula.
        
        If the expression uses non builtin objects, the evaluation context must be
        defined using `:func:setContext`.
        
        Arguments:
            expression (str): Expression of the formula.
        """
    
    def setVariables(self, varnames):
        """Define the variables names.
        
        Arguments:
            varnames (list[str]): List of variables names.
        """

# built-in function jeveux_init in libaster

def jeveux_init(mpi_comm):
    """Initialize the memory manager (Jeveux).
    
    Arguments:
        mpi_comm (int): Identifier of MPI communicator (from ``py2f()``).
    """

# built-in function jeveux_finalize in libaster

def jeveux_finalize(arg0):
    """Finalize the memory manager (Jeveux).
    """

# built-in function call_oper in libaster

def call_oper(syntax, jxveri):
    """Call a Fortran operator ('op' subroutine).
    
    Arguments:
        syntax (CommandSyntax): Object containing the user syntax.
        jxveri (int): If non null `JXVERI` is called after calling the operator.
    """

# built-in function call_oper_init in libaster

def call_oper_init():
    """Execute initializations before and after operator but without executing any
    operator.
    """

# built-in function call_ops in libaster

def call_ops(syntax, ops):
    """Call a Fortran 'ops' subroutine.
    
    Arguments:
        syntax (CommandSyntax): Object containing the user syntax.
        ops (int): Number of the `ops00x` subroutine.
    """

# built-in function call_debut in libaster

def call_debut(syntax):
    """Call a Fortran 'debut' subroutine.
    
    Arguments:
        syntax (CommandSyntax): Object containing the user syntax.
    """

# built-in function call_poursuite in libaster

def call_poursuite(syntax):
    """Call a Fortran 'poursuite' subroutine.
    
    Arguments:
        syntax (CommandSyntax): Object containing the user syntax.
    """

# built-in function cmd_ctxt_enter in libaster

def cmd_ctxt_enter():
    """Call Fortran 'cmd_ctxt_enter' subroutine.
    """

# built-in function cmd_ctxt_exit in libaster

def cmd_ctxt_exit():
    """Call Fortran 'cmd_ctxt_exit' subroutine.
    """

# built-in function write in libaster

def write(text):
    """Print a string using the fortran subroutine.
    
    Arguments:
        text (str): Text to be printed.
    """

# built-in function affich in libaster

def affich(code, text):
    """Print a string using the fortran subroutine on an internal file.
    
    Arguments:
        code (str): Code name of the internal file : 'MESSAGE' or 'CODE'.
        text (str): Text to be printed.
    """

# built-in function jeveux_status in libaster

def jeveux_status():
    """Return the status of the Jeveux memory manager.
    
    Returns:
        int: 0 after initialization and after shutdown, 1 during the execution.
    """

# built-in function jeveux_delete in libaster

def jeveux_delete(prefix):
    """Force manual deletion of Jeveux objects.
    
    Warning: Use only for debugging usages, it is dangerous for objects integrity
    and cpu consuming.
    
    Arguments:
        prefix (str): Root name of the Jeveux datastructure.
    """

# built-in function deleteTemporaryObjects in libaster

def deleteTemporaryObjects():
    """Delete temporary Jeveux objects
    """

# built-in function deleteCachedObjects in libaster

def deleteCachedObjects():
    """Delete temporary and cached Jeveux objects (temporary, matrix, base, ...)
    """

# built-in function onFatalError in libaster

def onFatalError(value= ''):
    """Get/set the behavior in case of error.
    
    Arguments:
        value (str, optional): Set the new behavior in case of error (one of "ABORT",
            "EXCEPTION", "EXCEPTION+VALID" or "INIT"). If `value` is not provided,
            the current behavior is returned.
    
    Returns:
        str: Current value
    """

# built-in function set_option in libaster

def set_option(arg0, arg1):
    """Set an option value to be used from Fortran operators.
    
    Arguments:
        option (str): Option name.
        value (float): Option value.
    """

# class Function2D in libaster

class Function2D(GenericFunction):
    pass
    
    # Method resolution order:
    #     Function2D
    #     GenericFunction
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __init__(self, *args, **kwargs):
        """Overloaded function.
        
        1. __init__(self: libaster.Function2D) -> None
        
        2. __init__(self: libaster.Function2D, arg0: str) -> None
        """
    
    def getParameters(self):
        """Return a list of the values of the parameter as (x1, x2, ...)
        
        Returns:
            list[float]: List of values (size = number of functions).
        """
    
    def getProperties(self):
        pass
    
    def getValues(self):
        """Return a list of the values of the functions as [F1, F2, ...]
        where Fi is (x1, x2, ..., y1, y2, ...).
        
        Returns:
            list[list [float]]: List of values (size = number of functions).
        """

# class Contact in libaster

class Contact(DataStructure):
    pass
    
    # Method resolution order:
    #     Contact
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __init__(self, *args, **kwargs):
        """Overloaded function.
        
        1. __init__(self: libaster.Contact, arg0: str, arg1: Model) -> None
        
        2. __init__(self: libaster.Contact, arg0: Model) -> None
        """
    
    def getFiniteElementDescriptor(self):
        pass
    
    def getModel(self):
        pass

# class ContactAlgo in libaster

class ContactAlgo:
    """Members:
    
    Lagrangian
    
    Nitsche
    
    Penalization
    """
    
    # Method resolution order:
    #     ContactAlgo
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __eq__(self, other):
        pass
    
    def __getstate__(self):
        pass
    
    def __hash__(self):
        pass
    
    def __index__(self):
        pass
    
    def __init__(self, value):
        pass
    
    def __int__(self):
        pass
    
    def __ne__(self, other):
        pass
    
    def __repr__(self):
        pass
    
    def __setstate__(self, state):
        pass
    
    def name(self):
        pass
    
    #----------------------------------------------------------------------
    # Data descriptors defined here:
    
    @property
    def __members__(self):
        pass
    
    @property
    def name(self):
        """name(self: handle) -> str
        """
    
    @property
    def value(self):
        pass
    
    #----------------------------------------------------------------------
    # Data and other attributes defined here:
    
    Lagrangian = 0
    
    Nitsche = 1
    
    Penalization = 2

# class ContactVariant in libaster

class ContactVariant:
    """Members:
    
    Empty
    
    Rapide
    
    Robust
    
    Symetric
    """
    
    # Method resolution order:
    #     ContactVariant
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __eq__(self, other):
        pass
    
    def __getstate__(self):
        pass
    
    def __hash__(self):
        pass
    
    def __index__(self):
        pass
    
    def __init__(self, value):
        pass
    
    def __int__(self):
        pass
    
    def __ne__(self, other):
        pass
    
    def __repr__(self):
        pass
    
    def __setstate__(self, state):
        pass
    
    def name(self):
        pass
    
    #----------------------------------------------------------------------
    # Data descriptors defined here:
    
    @property
    def __members__(self):
        pass
    
    @property
    def name(self):
        """name(self: handle) -> str
        """
    
    @property
    def value(self):
        pass
    
    #----------------------------------------------------------------------
    # Data and other attributes defined here:
    
    Empty = 0
    
    Rapide = 1
    
    Robust = 2
    
    Symetric = 3

# class ContactType in libaster

class ContactType:
    """Members:
    
    Unilateral
    
    Bilateral
    """
    
    # Method resolution order:
    #     ContactType
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __eq__(self, other):
        pass
    
    def __getstate__(self):
        pass
    
    def __hash__(self):
        pass
    
    def __index__(self):
        pass
    
    def __init__(self, value):
        pass
    
    def __int__(self):
        pass
    
    def __ne__(self, other):
        pass
    
    def __repr__(self):
        pass
    
    def __setstate__(self, state):
        pass
    
    def name(self):
        pass
    
    #----------------------------------------------------------------------
    # Data descriptors defined here:
    
    @property
    def __members__(self):
        pass
    
    @property
    def name(self):
        """name(self: handle) -> str
        """
    
    @property
    def value(self):
        pass
    
    #----------------------------------------------------------------------
    # Data and other attributes defined here:
    
    Bilateral = 1
    
    Unilateral = 0

# class FrictionAlgo in libaster

class FrictionAlgo:
    """Members:
    
    Lagrangian
    
    Nitsche
    
    Penalization
    """
    
    # Method resolution order:
    #     FrictionAlgo
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __eq__(self, other):
        pass
    
    def __getstate__(self):
        pass
    
    def __hash__(self):
        pass
    
    def __index__(self):
        pass
    
    def __init__(self, value):
        pass
    
    def __int__(self):
        pass
    
    def __ne__(self, other):
        pass
    
    def __repr__(self):
        pass
    
    def __setstate__(self, state):
        pass
    
    def name(self):
        pass
    
    #----------------------------------------------------------------------
    # Data descriptors defined here:
    
    @property
    def __members__(self):
        pass
    
    @property
    def name(self):
        """name(self: handle) -> str
        """
    
    @property
    def value(self):
        pass
    
    #----------------------------------------------------------------------
    # Data and other attributes defined here:
    
    Lagrangian = 0
    
    Nitsche = 1
    
    Penalization = 2

# class FrictionType in libaster

class FrictionType:
    """Members:
    
    Without
    
    Tresca
    
    Coulomb
    
    Stick
    """
    
    # Method resolution order:
    #     FrictionType
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __eq__(self, other):
        pass
    
    def __getstate__(self):
        pass
    
    def __hash__(self):
        pass
    
    def __index__(self):
        pass
    
    def __init__(self, value):
        pass
    
    def __int__(self):
        pass
    
    def __ne__(self, other):
        pass
    
    def __repr__(self):
        pass
    
    def __setstate__(self, state):
        pass
    
    def name(self):
        pass
    
    #----------------------------------------------------------------------
    # Data descriptors defined here:
    
    @property
    def __members__(self):
        pass
    
    @property
    def name(self):
        """name(self: handle) -> str
        """
    
    @property
    def value(self):
        pass
    
    #----------------------------------------------------------------------
    # Data and other attributes defined here:
    
    Coulomb = 2
    
    Stick = 3
    
    Tresca = 1
    
    Without = 0

# class PairingAlgo in libaster

class PairingAlgo:
    """Members:
    
    Mortar
    """
    
    # Method resolution order:
    #     PairingAlgo
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __eq__(self, other):
        pass
    
    def __getstate__(self):
        pass
    
    def __hash__(self):
        pass
    
    def __index__(self):
        pass
    
    def __init__(self, value):
        pass
    
    def __int__(self):
        pass
    
    def __ne__(self, other):
        pass
    
    def __repr__(self):
        pass
    
    def __setstate__(self, state):
        pass
    
    def name(self):
        pass
    
    #----------------------------------------------------------------------
    # Data descriptors defined here:
    
    @property
    def __members__(self):
        pass
    
    @property
    def name(self):
        """name(self: handle) -> str
        """
    
    @property
    def value(self):
        pass
    
    #----------------------------------------------------------------------
    # Data and other attributes defined here:
    
    Mortar = 0

# class InitialState in libaster

class InitialState:
    """Members:
    
    Interpenetrated
    
    Yes
    
    No
    """
    
    # Method resolution order:
    #     InitialState
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __eq__(self, other):
        pass
    
    def __getstate__(self):
        pass
    
    def __hash__(self):
        pass
    
    def __index__(self):
        pass
    
    def __init__(self, value):
        pass
    
    def __int__(self):
        pass
    
    def __ne__(self, other):
        pass
    
    def __repr__(self):
        pass
    
    def __setstate__(self, state):
        pass
    
    def name(self):
        pass
    
    #----------------------------------------------------------------------
    # Data descriptors defined here:
    
    @property
    def __members__(self):
        pass
    
    @property
    def name(self):
        """name(self: handle) -> str
        """
    
    @property
    def value(self):
        pass
    
    #----------------------------------------------------------------------
    # Data and other attributes defined here:
    
    Interpenetrated = 0
    
    No = 1
    
    Yes = 2

# class ContactParameter in libaster

class ContactParameter:
    pass
    
    # Method resolution order:
    #     ContactParameter
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __init__(self):
        pass
    
    def getAlgorithm(self):
        """Return the contact algorithm used. It is a value of an enum
        
        Returns:
            ContactAlgo: contact algorithm.
        """
    
    def getCoefficient(self):
        """Return the contact coefficient used. It is a value of a float
        
        Returns:
            float: contact coefficient.
        """
    
    def getType(self):
        """Return the contact type used. It is a value of an enum
        
        Returns:
            ContactType: contact type.
        """
    
    def getVariant(self):
        """Return the contact variant used. It is a value of an enum
        
        Returns:
            ContactVariant: contact variant.
        """
    
    def setAlgorithm(self, algo):
        """Set the contact algorithm used. It is a value of an enum
        
        Arguments:
            ContactAlgo: contact algorithm.
        """
    
    def setCoefficient(self, coeff):
        """Set the contact coefficient used. It is a value of a float
        
        Arguments:
            float: contact coefficient.
        """
    
    def setType(self, type):
        """Set the contact type used. It is a value of an enum
        
        Arguments:
            ContactType: contact type.
        """
    
    def setVariant(self, variant):
        """Set the contact variant used. It is a value of an enum
        
        Arguments:
            ContactVariant: contact variant.
        """

# class FrictionParameter in libaster

class FrictionParameter:
    pass
    
    # Method resolution order:
    #     FrictionParameter
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __init__(self):
        pass
    
    def getAlgorithm(self):
        """Return the Friction algorithm used. It is a value of an enum
        
        Returns:
            FrictionAlgo: Friction algorithm.
        """
    
    def getCoefficient(self):
        """Return the Friction coefficient used. It is a value of a float
        
        Returns:
            float: Friction coefficient.
        """
    
    def getCoulomb(self):
        """Return the Coulomb coefficient used. It is a value of a float
        
        Returns:
            float: Coulomb coefficient.
        """
    
    def getTresca(self):
        """Return the Tresca coefficient used. It is a value of a float
        
        Returns:
            float: Tresca coefficient.
        """
    
    def getType(self):
        """Return the Friction type used. It is a value of an enum
        
        Returns:
            FrictionType: Friction type.
        """
    
    def setAlgorithm(self, algo):
        """Set the Friction algorithm used. It is a value of an enum
        
        Arguments:
            FrictionAlgo: Friction algorithm.
        """
    
    def setCoefficient(self, coeff):
        """Set the Friction coefficient used. It is a value of a float
        
        Arguments:
            float: Friction coefficient.
        """
    
    def setCoulomb(self, coulomb):
        """Set the Coulomb coefficient used. It is a value of a float
        
        Arguments:
            float: Coulomb coefficient.
        """
    
    def setTresca(self, tresca):
        """Set the Tresca coefficient used. It is a value of a float
        
        Arguments:
            float: Tresca coefficient.
        """
    
    def setType(self, type):
        """Set the Friction type used. It is a value of an enum
        
        Arguments:
            FrictionType: Friction type.
        """
    
    #----------------------------------------------------------------------
    # Data descriptors defined here:
    
    @property
    def hasFriction(self):
        """bool: enable or disable the use of friction.
        """

# class PairingParameter in libaster

class PairingParameter:
    pass
    
    # Method resolution order:
    #     PairingParameter
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __init__(self):
        pass
    
    def getAlgorithm(self):
        """Return the Pairing algorithm used. It is a value of an enum
        
        Returns:
            PairingAlgo: Pairing algorithm.
        """
    
    def getDistanceFunction(self):
        """Return the fictive distance function. It is a value of a pointer
        
        Returns:
            GenericFunction: FunctionPtr/ FormulaPtr/ Function2DPtr.
        """
    
    def getDistanceRatio(self):
        """Return the pairing distance ratio used. It is a value of a float
        
        Returns:
            float: pairing distance.
        """
    
    def getElementaryCharacteristics(self):
        """Return the elementary characteristics. It is a value of a pointer
        
        Returns:
            ElementaryCharacteristicsPtr: cara_elel pointer.
        """
    
    def getInitialState(self):
        """Return the initial contact state. It is a value of an enum
        
        Returns:
            InitialState: Initial contact state.
        """
    
    def getThreshold(self):
        """Return the distance threshold. It is a value of a float
        
        Returns:
            float: distance threshold.
        """
    
    def setAlgorithm(self, algo):
        """Set the Pairing algorithm used. It is a value of an enum
        
        Arguments:
            PairingAlgo: Pairing algorithm.
        """
    
    def setDistanceFunction(self, dist_supp):
        """Set the fictive distance function. It is a value of a pointer
        
        Arguments:
            GenericFunction: FunctionPtr/ FormulaPtr/ Function2DPtr.
        """
    
    def setDistanceRatio(self, dist_ratio):
        """Set the pairing distance ratio used. It is a value of a float
        
        Arguments:
            float: pairing distance ratio.
        """
    
    def setElementaryCharacteristics(self, cara):
        """Set the elementary characteristics. It is a value of a pointer
        
        Arguments:
            ElementaryCharacteristicsPtr: cara_elel pointer.
        """
    
    def setInitialState(self, cont_init):
        """Set the initial contact state. It is a value of an enum
        
        Arguments:
            InitialState: Initial contact state.
        """
    
    def setThreshold(self, seuil):
        """Set the distance threshold used. It is a value of a float
        
        Arguments:
            float: distance threshold.
        """
    
    #----------------------------------------------------------------------
    # Data descriptors defined here:
    
    @property
    def hasBeamDistance(self):
        """bool: enable or disable the use of a fictive distance for beam.
        """
    
    @property
    def hasShellDistance(self):
        """bool: enable or disable the use of a fictive distance for shell.
        """

# class ContactNew in libaster

class ContactNew(DataStructure):
    pass
    
    # Method resolution order:
    #     ContactNew
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __init__(self, *args, **kwargs):
        """Overloaded function.
        
        1. __init__(self: libaster.ContactNew, arg0: str, arg1: Model) -> None
        
        2. __init__(self: libaster.ContactNew, arg0: Model) -> None
        """
    
    def appendContactZone(self, zone):
        """Append a new contact zone to the contact definition
        
        Arguments:
            zone (ContactZone): contact zone to append
        """
    
    def build(self):
        """Build and check internal objects
        """
    
    def getContactZone(self, zone_id):
        """Return the specified contact zone
        
        Arguments:
            zone_id (int): index of the contact zone (0-based)
        
        Returns:
            ContactZone: contact zone.
        """
    
    def getContactZones(self):
        """Return the list of contact zones
        
        Returns:
            List[ContactZone]: List of contact zones.
        """
    
    def getFiniteElementDescriptor(self):
        """Return the finite element descriptor to define virtual cells for Lagrange multipliers
        
        Returns:
            FiniteElementDescriptor: fed.
        """
    
    def getMesh(self):
        """Return the mesh used in the contact definition
        
        Returns:
            Mesh: mesh.
        """
    
    def getModel(self):
        """Return the model used in the contact definition
        
        Returns:
            Model: model.
        """
    
    def getNumberOfContactZones(self):
        """Return the number of contact zones used
        
        Returns:
            inter: number of contact zones.
        """
    
    def getVerbosity(self):
        """Get level of verbosity:*
              0- without
              1- normal
              2- detailled
        
        Returns:
            integer: level of verbosity
        """
    
    def setVerbosity(self, level):
        """Set level of verbosity:
              0- without
              1- normal (default)
              2- detailled
        
        Arguments:
            level (int): level of verbosity
        """
    
    #----------------------------------------------------------------------
    # Data descriptors defined here:
    
    @property
    def hasFriction(self):
        """bool: enable or disable the use of friction.
        """
    
    @property
    def hasSmoothing(self):
        """bool: enable or disable  the use of smoothing.
        """

# class FrictionNew in libaster

class FrictionNew(ContactNew):
    pass
    
    # Method resolution order:
    #     FrictionNew
    #     ContactNew
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __init__(self, *args, **kwargs):
        """Overloaded function.
        
        1. __init__(self: libaster.FrictionNew, arg0: str, arg1: Model) -> None
        
        2. __init__(self: libaster.FrictionNew, arg0: Model) -> None
        """

# class ContactZone in libaster

class ContactZone(DataStructure):
    pass
    
    # Method resolution order:
    #     ContactZone
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __init__(self, *args, **kwargs):
        """Overloaded function.
        
        1. __init__(self: libaster.ContactZone, arg0: str, arg1: Model) -> None
        
        2. __init__(self: libaster.ContactZone, arg0: Model) -> None
        """
    
    def build(self):
        """Build and check internal objects
        """
    
    def getContactParameter(self):
        """Get contact parameters defining method, coefficient...
        
        Returns:
            ContactParameter: contact parameters
        """
    
    def getExcludedSlaveCells(self):
        """Get excluded groups of cells on slave side
        
        Returns:
            str: excluded groups' names
        """
    
    def getFrictionParameter(self):
        """Get friction parameters defining method, coefficient...
        
        Returns:
            FrictionParameter: friction parameters
        """
    
    def getMasterCellNeighbors(self, cell_number):
        """Get the master cells in the neighbor of a given master cell number
        
        Arguments:
            int: master cell number
        """
    
    def getMasterCellsFromNode(self, node_number):
        """Get the master cells associtaed with a node number
        
        Arguments:
            int: node number
        """
    
    def getMesh(self):
        """Return the mesh used in the contact zone definition
        
        Returns:
            BaseMesh: mesh.
        """
    
    def getModel(self):
        """Return the model used in the contact zone definition
        
        Returns:
            Model: model.
        """
    
    def getPairingParameter(self):
        """Get pairing parameters defining algorithm, distance...
        
        Returns:
            PairingParameter: pairing parameters
        """
    
    def getSlaveCellNeighbors(self, cell_number):
        """Get the slave cells in the neighbor of a given slave cell number
        
        Arguments:
            int: slave cell number
        """
    
    def getSlaveCells(self):
        """Get slave's cells index
        
        Returns:
            list[int]: slave's cells index
        """
    
    def getSlaveCellsFromNode(self, node_number):
        """Get the slave cells associtaed with a node number
        
        Arguments:
            int: node number
        """
    
    def getSlaveNodes(self):
        """Get slave's nodes index
        
        Returns:
            list[int]: slave's nodes index
        """
    
    def getVerbosity(self):
        """Get level of verbosity:
              0- without
              1- normal
              2- detailled
        
        Returns:
            integer: level of verbosity
        """
    
    def setContactParameter(self, contact):
        """Set contact parameters defining method, coefficient...
        
        Arguments:
            ContactParameter: contact parameters
        """
    
    def setExcludedSlaveGroupOfCells(self, groups):
        """Set excluded groups of cells on slave side
        
        Arguments:
            str: excluded groups' names
        """
    
    def setExcludedSlaveGroupOfNodes(self, groups):
        """Set excluded groups of nodes on slave side
        
        Arguments:
            str: excluded groups' names
        """
    
    def setFrictionParameter(self, friction):
        """Set friction parameters defining method, coefficient...
        
        Arguments:
            FrictionParameter: friction parameters
        """
    
    def setMasterGroupOfCells(self, master_name):
        """Set master's name of group of cells
        
        Arguments:
            str: master's name
        """
    
    def setPairingParameter(self, pairing):
        """Set pairing parameters defining algorithm, distance...
        
        Arguments:
            PairingParameter: pairing parameters
        """
    
    def setSlaveGroupOfCells(self, slave_name):
        """Set slave's name of group of cells
        
        Arguments:
            str: slave's name
        """
    
    def setVerbosity(self, level):
        """Set level of verbosity:
              0- without
              1- normal (default)
              2- detailled
        
        Arguments:
            integer: level of verbosity
        """
    
    #----------------------------------------------------------------------
    # Data descriptors defined here:
    
    @property
    def checkNormals(self):
        """bool: Attribute that holds the checking of outwards normals.
        """
    
    @property
    def hasFriction(self):
        """bool: enable or disable the use of friction.
        """
    
    @property
    def hasSmoothing(self):
        """bool: enable or disable  the use of smoothing.
        """

# class ContactPairing in libaster

class ContactPairing(DataStructure):
    pass
    
    # Method resolution order:
    #     ContactPairing
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __init__(self, *args, **kwargs):
        """Overloaded function.
        
        1. __init__(self: libaster.ContactPairing, arg0: str, arg1: libaster.ContactNew) -> None
        
        2. __init__(self: libaster.ContactPairing, arg0: libaster.ContactNew) -> None
        """
    
    def clear(self):
        """clean all the paring quantities of all zones
        Returns:
            bool: true if the pairing quantities are cleared
        """
    
    def clearZone(self, zone_index):
        """clean all the paring quantities of zone zonde_index
        Arguments:
            zone_index(int)
        Returns:
            bool: true if the pairing quantities are cleared
        """
    
    def compute(self):
        """Compute the pairing quantities associated with the zones
        
        Returns:
            bool: True if the pairing quantities are updated appropriately
        """
    
    def computeZone(self, zone_index):
        """Compute the pairing quantities associated with the zone zone_index
        Arguments:
            zone_index(int)
        Returns:
            bool: True if the pairing quantities are updated appropriately
        """
    
    def getCoordinates(self):
        """Coordinates of nodes used for pairing (almost always different from the intial mesh).
        
        Returns:
            MeshCoordinatesFieldPtr: the coordinates field
        """
    
    def getFiniteElementDescriptor(self):
        """Return Finite Element Descriptor for virtual cells from pairing.
        
        Returns:
            FiniteElementDescriptor: finite element for virtual cells
        """
    
    def getListOfPairsOfZone(self, zone_index):
        """return list of pairs associated with the zone izone
        Arguments:
            zone_index(int)
        Returns:
            List[List[int]]: List of pairs
        """
    
    def getNumberOfPairs(self):
        """return the total number of pairs
        Returns:
            int: Total number of pairs
        """
    
    def getNumberOfPairsOfZone(self, zone_index):
        """return number of  pairs associated with the zone zone_index
        Arguments:
            zone_index(int)
        Returns:
            int: number of pairs
        """
    
    def getSlaveIntersectionPoints(self, zone_index):
        """Get the intersection points beetween a master and slave cells in the parametric
        slave space. The maximum number of points is 8.
        
        Arguments:
            zone_index(int) : index of zone
        
        Returns:
            list[list]: list of list of intersection points (each intersection is of size 16)
        """
    
    def setCoordinates(self, coordinates):
        """Set the mesh coordinates field
        
        Arguments:
            coordinates (MeshCoordinatesField) : coordinates to use for pairing
        """
    
    def updateCoordinates(self, disp):
        """Update the mesh coordinates given a displacement field
        """

# class ContactComputation in libaster

class ContactComputation:
    pass
    
    # Method resolution order:
    #     ContactComputation
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __init__(self, arg0):
        pass
    
    def contactCoefficient(self):
        """Compute contact coefficient at the nodes of the slave surface based on values of COEF_CONT
        and COEF_FROT
        
        Returns:
            FieldOnNodesReal: contact coefficient (= COEF_CONT)
            FieldOnNodesReal: friction coefficient (= COEF_FROT)
        """
    
    def contactData(self, pairing, material, initial_contact):
        """Compute contact data (cf. MMCHML) as input to compute contact forces and matrices.
        
        Arguments:
            pairing (ContactPairing): pairing object
            material (MaterialField): material field
            initial_contact (bool): True to use value in contact definition (CONTACT_INIT).
        
        Returns:
            FieldOnCellsReal: contact data
        """
    
    def geometricGap(self, coordinates):
        """Compute geometric gap and indicator using projection. The indicator is equal to 0 for
        a node with no projection (gap value is Nan) found else 1.
        
        Arguments:
            coordinates (MeshCoordinatesField): (current) coordinates of mesh
        
        Returns:
            FieldOnNodesReal: gap field.
            FieldOnNodesReal: gap indicator.
        """

# class BaseAssemblyMatrix in libaster

class BaseAssemblyMatrix(DataStructure):
    pass
    
    # Method resolution order:
    #     BaseAssemblyMatrix
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __init__(self, *args, **kwargs):
        """Overloaded function.
        
        1. __init__(self: libaster.BaseAssemblyMatrix, arg0: str) -> None
        
        2. __init__(self: libaster.BaseAssemblyMatrix, arg0: str, arg1: str) -> None
        
        3. __init__(self: libaster.BaseAssemblyMatrix, arg0: PhysicalProblem, arg1: str) -> None
        """
    
    def addDirichletBC(self, currentLoad, func):
        pass
    
    def assemble(self, clean= True):
        """Assembly matrix from elementar matrices added.
        
        Arguments:
            clean (bool) : Clean elementary matrices after building (default = true)
        """
    
    def getDOFNumbering(self):
        pass
    
    def getDirichletBCDOFs(self):
        """Return a vector with DOFs eliminated by Dirichlet boundaries conditions (if it exists)
        
        Returns:
            tuple(int): a vector with DOFs eliminated by Dirichlet boundaries conditions of
                size = neq + 1,
                tuple(ieq = 0, neq - 1) = 1 then DOF eliminated else 0,
                tuple(neq) = number of DOFs eliminated.
        """
    
    def getLagrangeScaling(self):
        """Return the scaling used for Lagrange multipliers. It returns 1 if no Lagrange.
        
        Returns:
            float: scaling used for Lagrange multipliers. It returns 1 if no Lagrange
            are present.
        """
    
    def getListOfLoads(self):
        """Return the list of loads.
        
        Returns:
            ListOfLoads: a pointer to the list of loads
        """
    
    def getMesh(self):
        """Return the mesh.
        
        Returns:
            Mesh: a pointer to the mesh
        """
    
    def getModel(self):
        """Return the model.
        
        Returns:
            Model: a pointer to the model
        """
    
    def hasDirichletEliminationDOFs(self):
        """Tell if matrix has some DOFs eliminated by Dirichlet boundaries conditions.
        
        Returns:
            bool: *True* if matrix has some DOFs eliminated by Dirichlet boundaries conditions else *False*
        """
    
    def isEmpty(self):
        """Tell if the matrix is empty.
        
        Returns:
            bool: *True* if the matrix is empty.
        """
    
    def isFactorized(self):
        """Tell if the matrix is factorized.
        
        Returns:
            bool: *True* if the matrix is factorized, *False* otherwise.
        """
    
    def print(self, format= 'ASTER', unit= 6):
        """Print the matrix in code_aster or matlab format (with information on the DOF).
        
        Arguments:
            format (str): 'ASTER' (default) or 'MATLAB'
        """
    
    def setDOFNumbering(self, arg0):
        pass
    
    def setListOfLoads(self, load):
        """Set the list of loads.
        
        Arguments:
            ListOfLoads: a pointer to the list of loads to set
        """
    
    def setSolverName(self, arg0):
        pass
    
    def symmetrize(self):
        """Make the assembly matrix symmetric in place
        """
    
    def transpose(self):
        pass
    
    def updateDOFNumbering(self):
        pass

# class AssemblyMatrixDisplacementReal in libaster

class AssemblyMatrixDisplacementReal(BaseAssemblyMatrix):
    pass
    
    # Method resolution order:
    #     AssemblyMatrixDisplacementReal
    #     BaseAssemblyMatrix
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __add__(self, arg0):
        pass
    
    def __iadd__(self, arg0):
        pass
    
    def __imul__(self, arg0):
        pass
    
    def __init__(self, *args, **kwargs):
        """Overloaded function.
        
        1. __init__(self: libaster.AssemblyMatrixDisplacementReal) -> None
        
        2. __init__(self: libaster.AssemblyMatrixDisplacementReal, arg0: str) -> None
        
        3. __init__(self: libaster.AssemblyMatrixDisplacementReal, arg0: PhysicalProblem) -> None
        
        4. __init__(self: libaster.AssemblyMatrixDisplacementReal, arg0: libaster.AssemblyMatrixDisplacementReal) -> None
        """
    
    def __isub__(self, arg0):
        pass
    
    def __mul__(self, arg0):
        pass
    
    def __neg__(self):
        pass
    
    def __rmul__(self, arg0):
        pass
    
    def __sub__(self, arg0):
        pass
    
    def addElementaryMatrix(self, matr_elem, coeff= 1.0):
        """Add elementary matrix to assemble such that during assembling Mat += coeff * matr_elem
        
        Arguments:
            matr_elem [ElementaryMatrixDisplacementReal]: elementary matrix to add
            coeff [float]: assembling factor (default = 1.0)
        """
    
    def clearElementaryMatrix(self):
        pass
    
    def defineSolver(self):
        pass
    
    def duplicate(self):
        pass
    
    def getMaterialField(self):
        pass
    
    def getNumberOfElementaryMatrix(self):
        pass
    
    def scale(self, arg0, arg1):
        """Scale the matrix in place using right and left multiplication by diagonal matrices stored as vectors
        
        Arguments:
            lvect (list[float]): List of the values.
            rvect (list[float]): List of the values.
        """
    
    def setValues(self, arg0, arg1, arg2):
        """Erase the assembly matrix and set new values in it.
        
        The new values are in coordinate format (i, j, aij). The matrix  must be stored in CSR format.
        There is no rule for the indices - they can be in arbitrary order and can be repeated. Repeated
        indices are sumed according to an assembly process.
        
        Arguments:
            idx (list[int]): List of the row indices.
            jdx (list[int]): List of the column indices.
            values (list[float]): List of the values.
        """
    
    def size(self, local= True):
        """Get the size of the matrix
        
        Arguments:
            local (bool) local or global size
        """

# class AssemblyMatrixDisplacementComplex in libaster

class AssemblyMatrixDisplacementComplex(BaseAssemblyMatrix):
    pass
    
    # Method resolution order:
    #     AssemblyMatrixDisplacementComplex
    #     BaseAssemblyMatrix
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __imul__(self, arg0):
        pass
    
    def __init__(self, *args, **kwargs):
        """Overloaded function.
        
        1. __init__(self: libaster.AssemblyMatrixDisplacementComplex) -> None
        
        2. __init__(self: libaster.AssemblyMatrixDisplacementComplex, arg0: str) -> None
        """
    
    def __mul__(self, arg0):
        pass
    
    def __neg__(self):
        pass
    
    def __rmul__(self, arg0):
        pass
    
    def addElementaryMatrix(self, matr_elem, coeff= 1.0):
        """Add elementary matrix to assemble such that during assembling Mat += coeff * matr_elem
        
        Arguments:
            matr_elem [ElementaryMatrixDisplacementComplex]: elementary matrix to add
            coeff [float]: assembling factor (default = 1.0)
        """
    
    def clearElementaryMatrix(self):
        pass
    
    def defineSolver(self):
        pass
    
    def duplicate(self):
        pass
    
    def getMaterialField(self):
        pass
    
    def getNumberOfElementaryMatrix(self):
        pass
    
    def setValues(self, arg0, arg1, arg2):
        """Erase the assembly matrix and set new values in it.
        
        The new values are in coordinate format (i, j, aij). The matrix  must be stored in CSR format.
        There is no rule for the indices - they can be in arbitrary order and can be repeated. Repeated
        indices are sumed according to an assembly process.
        
        Arguments:
            idx (list[int]): List of the row indices.
            jdx (list[int]): List of the column indices.
            values (list[float]): List of the values.
        """
    
    def transposeConjugate(self):
        pass

# class AssemblyMatrixTemperatureReal in libaster

class AssemblyMatrixTemperatureReal(BaseAssemblyMatrix):
    pass
    
    # Method resolution order:
    #     AssemblyMatrixTemperatureReal
    #     BaseAssemblyMatrix
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __add__(self, arg0):
        pass
    
    def __iadd__(self, arg0):
        pass
    
    def __imul__(self, arg0):
        pass
    
    def __init__(self, *args, **kwargs):
        """Overloaded function.
        
        1. __init__(self: libaster.AssemblyMatrixTemperatureReal) -> None
        
        2. __init__(self: libaster.AssemblyMatrixTemperatureReal, arg0: str) -> None
        
        3. __init__(self: libaster.AssemblyMatrixTemperatureReal, arg0: PhysicalProblem) -> None
        """
    
    def __isub__(self, arg0):
        pass
    
    def __mul__(self, arg0):
        pass
    
    def __neg__(self):
        pass
    
    def __rmul__(self, arg0):
        pass
    
    def __sub__(self, arg0):
        pass
    
    def addElementaryMatrix(self, matr_elem, coeff= 1.0):
        """Add elementary matrix to assemble such that during assembling Mat += coeff * matr_elem
        
        Arguments:
            matr_elem [ElementaryMatrixTemperatureReal]: elementary matrix to add
            coeff [float]: assembling factor (default = 1.0)
        """
    
    def clearElementaryMatrix(self):
        pass
    
    def defineSolver(self):
        pass
    
    def duplicate(self):
        pass
    
    def getMaterialField(self):
        pass
    
    def getNumberOfElementaryMatrix(self):
        pass
    
    def scale(self, arg0, arg1):
        """Scale the matrix in place using right and left multiplication by diagonal matrices stored as vectors
        
        Arguments:
            lvect (list[float]): List of the values.
            rvect (list[float]): List of the values.
        """
    
    def setValues(self, arg0, arg1, arg2):
        """Erase the assembly matrix and set new values in it.
        
        The new values are in coordinate format (i, j, aij). The matrix  must be stored in CSR format.
        There is no rule for the indices - they can be in arbitrary order and can be repeated. Repeated
        indices are sumed according to an assembly process.
        
        Arguments:
            idx (list[int]): List of the row indices.
            jdx (list[int]): List of the column indices.
            values (list[float]): List of the values.
        """
    
    def size(self, local= True):
        """Get the size of the matrix
        
        Arguments:
            local (bool) local or global size
        """

# class AssemblyMatrixTemperatureComplex in libaster

class AssemblyMatrixTemperatureComplex(BaseAssemblyMatrix):
    pass
    
    # Method resolution order:
    #     AssemblyMatrixTemperatureComplex
    #     BaseAssemblyMatrix
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __init__(self, *args, **kwargs):
        """Overloaded function.
        
        1. __init__(self: libaster.AssemblyMatrixTemperatureComplex) -> None
        
        2. __init__(self: libaster.AssemblyMatrixTemperatureComplex, arg0: str) -> None
        """
    
    def addElementaryMatrix(self, matr_elem, coeff= 1.0):
        """Add elementary matrix to assemble such that during assembling Mat += coeff * matr_elem
        
        Arguments:
            matr_elem [ElementaryMatrixDisplacementReal]: elementary matrix to add
            coeff [float]: assembling factor (default = 1.0)
        """
    
    def clearElementaryMatrix(self):
        pass
    
    def getMaterialField(self):
        pass
    
    def getNumberOfElementaryMatrix(self):
        pass
    
    def transposeConjugate(self):
        pass

# class AssemblyMatrixPressureReal in libaster

class AssemblyMatrixPressureReal(BaseAssemblyMatrix):
    pass
    
    # Method resolution order:
    #     AssemblyMatrixPressureReal
    #     BaseAssemblyMatrix
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __add__(self, arg0):
        pass
    
    def __iadd__(self, arg0):
        pass
    
    def __imul__(self, arg0):
        pass
    
    def __init__(self, *args, **kwargs):
        """Overloaded function.
        
        1. __init__(self: libaster.AssemblyMatrixPressureReal) -> None
        
        2. __init__(self: libaster.AssemblyMatrixPressureReal, arg0: str) -> None
        """
    
    def __isub__(self, arg0):
        pass
    
    def __mul__(self, arg0):
        pass
    
    def __neg__(self):
        pass
    
    def __rmul__(self, arg0):
        pass
    
    def __sub__(self, arg0):
        pass
    
    def addElementaryMatrix(self, matr_elem, coeff= 1.0):
        """Add elementary matrix to assemble such that during assembling Mat += coeff * matr_elem
        
        Arguments:
            matr_elem [ElementaryMatrixDisplacementReal]: elementary matrix to add
            coeff [float]: assembling factor (default = 1.0)
        """
    
    def clearElementaryMatrix(self):
        pass
    
    def duplicate(self):
        pass
    
    def getMaterialField(self):
        pass
    
    def getNumberOfElementaryMatrix(self):
        pass
    
    def setValues(self, arg0, arg1, arg2):
        """Erase the assembly matrix and set new values in it.
        
        The new values are in coordinate format (i, j, aij). The matrix  must be stored in CSR format.
        There is no rule for the indices - they can be in arbitrary order and can be repeated. Repeated
        indices are sumed according to an assembly process.
        
        Arguments:
            idx (list[int]): List of the row indices.
            jdx (list[int]): List of the column indices.
            values (list[float]): List of the values.
        """

# class AssemblyMatrixPressureComplex in libaster

class AssemblyMatrixPressureComplex(BaseAssemblyMatrix):
    pass
    
    # Method resolution order:
    #     AssemblyMatrixPressureComplex
    #     BaseAssemblyMatrix
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __add__(self, arg0):
        pass
    
    def __iadd__(self, arg0):
        pass
    
    def __imul__(self, arg0):
        pass
    
    def __init__(self, *args, **kwargs):
        """Overloaded function.
        
        1. __init__(self: libaster.AssemblyMatrixPressureComplex) -> None
        
        2. __init__(self: libaster.AssemblyMatrixPressureComplex, arg0: str) -> None
        """
    
    def __isub__(self, arg0):
        pass
    
    def __mul__(self, arg0):
        pass
    
    def __neg__(self):
        pass
    
    def __rmul__(self, arg0):
        pass
    
    def __sub__(self, arg0):
        pass
    
    def addElementaryMatrix(self, matr_elem, coeff= 1.0):
        """Add elementary matrix to assemble such that during assembling Mat += coeff * matr_elem
        
        Arguments:
            matr_elem [ElementaryMatrixPressureComplex]: elementary matrix to add
            coeff [float]: assembling factor (default = 1.0)
        """
    
    def clearElementaryMatrix(self):
        pass
    
    def defineSolver(self):
        pass
    
    def duplicate(self):
        pass
    
    def getMaterialField(self):
        pass
    
    def getNumberOfElementaryMatrix(self):
        pass
    
    def transposeConjugate(self):
        pass

# class ElementaryTermReal in libaster

class ElementaryTermReal(DataField):
    pass
    
    # Method resolution order:
    #     ElementaryTermReal
    #     DataField
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __init__(self, *args, **kwargs):
        """Overloaded function.
        
        1. __init__(self: libaster.ElementaryTermReal) -> None
        
        2. __init__(self: libaster.ElementaryTermReal, arg0: str) -> None
        """
    
    def getFiniteElementDescriptor(self):
        """Return the finite element descriptor
        
        Returns:
            FiniteElementDescriptor: finite element descriptor
        """
    
    def getLocalMode(self):
        """Return the local mode.
        
        Returns:
            str: the local mode
        """
    
    def getMesh(self):
        """Return the mesh
        
        Returns:
            BaseMesh: a pointer to the mesh
        """
    
    def getOption(self):
        """Return the optior used to compute it
        
        Returns:
            str: name of the option
        """
    
    def getPhysicalQuantity(self):
        """Return the physical quantity
        
        Returns:
            str: name of the physical quantity
        """

# class ElementaryTermComplex in libaster

class ElementaryTermComplex(DataField):
    pass
    
    # Method resolution order:
    #     ElementaryTermComplex
    #     DataField
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __init__(self, *args, **kwargs):
        """Overloaded function.
        
        1. __init__(self: libaster.ElementaryTermComplex) -> None
        
        2. __init__(self: libaster.ElementaryTermComplex, arg0: str) -> None
        """
    
    def getFiniteElementDescriptor(self):
        """Return the finite element descriptor
        
        Returns:
            FiniteElementDescriptor: finite element descriptor
        """
    
    def getLocalMode(self):
        """Return the local mode.
        
        Returns:
            str: the local mode
        """
    
    def getMesh(self):
        """Return the mesh
        
        Returns:
            BaseMesh: a pointer to the mesh
        """
    
    def getOption(self):
        """Return the optior used to compute it
        
        Returns:
            str: name of the option
        """
    
    def getPhysicalQuantity(self):
        """Return the physical quantity
        
        Returns:
            str: name of the physical quantity
        """

# class BaseElementaryMatrix in libaster

class BaseElementaryMatrix(DataStructure):
    pass
    
    # Method resolution order:
    #     BaseElementaryMatrix
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __init__(self,  *args, **kwargs):
        """Initialize self.  See help(type(self)) for accurate signature.
        """
    
    def getElementaryCharacteristics(self):
        pass
    
    def getMaterialField(self):
        pass
    
    def getMesh(self):
        pass
    
    def getModel(self):
        pass
    
    def setElementaryCharacteristics(self, arg0):
        pass
    
    def setMaterialField(self, arg0):
        pass
    
    def setModel(self, arg0):
        pass
    
    def setPhysicalProblem(self, arg0):
        pass

# class ElementaryMatrixDisplacementReal in libaster

class ElementaryMatrixDisplacementReal(BaseElementaryMatrix):
    pass
    
    # Method resolution order:
    #     ElementaryMatrixDisplacementReal
    #     BaseElementaryMatrix
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __imul__(self, arg0):
        pass
    
    def __init__(self, *args, **kwargs):
        """Overloaded function.
        
        1. __init__(self: libaster.ElementaryMatrixDisplacementReal) -> None
        
        2. __init__(self: libaster.ElementaryMatrixDisplacementReal, arg0: str) -> None
        """
    
    def addElementaryTerm(self, *args, **kwargs):
        """Overloaded function.
        
        1. addElementaryTerm(self: libaster.ElementaryMatrixDisplacementReal, arg0: libaster.ElementaryTermReal) -> None
        
        2. addElementaryTerm(self: libaster.ElementaryMatrixDisplacementReal, arg0: List[libaster.ElementaryTermReal]) -> None
        """
    
    def build(self):
        pass
    
    def getElementaryTerms(self):
        pass
    
    def getFiniteElementDescriptors(self):
        pass
    
    def hasElementaryTerms(self):
        pass

# class ElementaryMatrixDisplacementComplex in libaster

class ElementaryMatrixDisplacementComplex(BaseElementaryMatrix):
    pass
    
    # Method resolution order:
    #     ElementaryMatrixDisplacementComplex
    #     BaseElementaryMatrix
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __init__(self, *args, **kwargs):
        """Overloaded function.
        
        1. __init__(self: libaster.ElementaryMatrixDisplacementComplex) -> None
        
        2. __init__(self: libaster.ElementaryMatrixDisplacementComplex, arg0: str) -> None
        """
    
    def addElementaryTerm(self, *args, **kwargs):
        """Overloaded function.
        
        1. addElementaryTerm(self: libaster.ElementaryMatrixDisplacementComplex, arg0: libaster.ElementaryTermComplex) -> None
        
        2. addElementaryTerm(self: libaster.ElementaryMatrixDisplacementComplex, arg0: List[libaster.ElementaryTermComplex]) -> None
        """
    
    def build(self):
        pass
    
    def getElementaryTerms(self):
        pass
    
    def getFiniteElementDescriptors(self):
        pass
    
    def hasElementaryTerms(self):
        pass

# class ElementaryMatrixTemperatureReal in libaster

class ElementaryMatrixTemperatureReal(BaseElementaryMatrix):
    pass
    
    # Method resolution order:
    #     ElementaryMatrixTemperatureReal
    #     BaseElementaryMatrix
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __imul__(self, arg0):
        pass
    
    def __init__(self, *args, **kwargs):
        """Overloaded function.
        
        1. __init__(self: libaster.ElementaryMatrixTemperatureReal) -> None
        
        2. __init__(self: libaster.ElementaryMatrixTemperatureReal, arg0: str) -> None
        """
    
    def addElementaryTerm(self, *args, **kwargs):
        """Overloaded function.
        
        1. addElementaryTerm(self: libaster.ElementaryMatrixTemperatureReal, arg0: libaster.ElementaryTermReal) -> None
        
        2. addElementaryTerm(self: libaster.ElementaryMatrixTemperatureReal, arg0: List[libaster.ElementaryTermReal]) -> None
        """
    
    def build(self):
        pass
    
    def getElementaryTerms(self):
        pass
    
    def getFiniteElementDescriptors(self):
        pass
    
    def hasElementaryTerms(self):
        pass

# class ElementaryMatrixPressureComplex in libaster

class ElementaryMatrixPressureComplex(BaseElementaryMatrix):
    pass
    
    # Method resolution order:
    #     ElementaryMatrixPressureComplex
    #     BaseElementaryMatrix
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __init__(self, *args, **kwargs):
        """Overloaded function.
        
        1. __init__(self: libaster.ElementaryMatrixPressureComplex) -> None
        
        2. __init__(self: libaster.ElementaryMatrixPressureComplex, arg0: str) -> None
        """
    
    def addElementaryTerm(self, *args, **kwargs):
        """Overloaded function.
        
        1. addElementaryTerm(self: libaster.ElementaryMatrixPressureComplex, arg0: libaster.ElementaryTermComplex) -> None
        
        2. addElementaryTerm(self: libaster.ElementaryMatrixPressureComplex, arg0: List[libaster.ElementaryTermComplex]) -> None
        """
    
    def build(self):
        pass
    
    def getElementaryTerms(self):
        pass
    
    def getFiniteElementDescriptors(self):
        pass
    
    def hasElementaryTerms(self):
        pass

# class BaseElementaryVector in libaster

class BaseElementaryVector(DataStructure):
    pass
    
    # Method resolution order:
    #     BaseElementaryVector
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __init__(self, *args, **kwargs):
        """Overloaded function.
        
        1. __init__(self: libaster.BaseElementaryVector) -> None
        
        2. __init__(self: libaster.BaseElementaryVector, arg0: str) -> None
        
        3. __init__(self: libaster.BaseElementaryVector, arg0: PhysicalProblem) -> None
        """
    
    def addLoad(self, arg0):
        pass
    
    def assembleWithLoadFunctions(self, dofNume, time= 0.0):
        pass
    
    def assembleWithMask(self, arg0, arg1, arg2):
        pass
    
    def build(self):
        pass
    
    def setElementaryCharacteristics(self, arg0):
        pass
    
    def setListOfLoads(self, arg0):
        pass
    
    def setMaterialField(self, arg0):
        pass
    
    def setModel(self, arg0):
        pass
    
    def setPhysicalProblem(self, phys_pb):
        """Set the physical problem
        
        Arguments:
            phys_pb (PhysicalProblem): the physical problem.
        """
    
    def setType(self, arg0):
        pass

# class ElementaryVectorReal in libaster

class ElementaryVectorReal(BaseElementaryVector):
    pass
    
    # Method resolution order:
    #     ElementaryVectorReal
    #     BaseElementaryVector
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __init__(self, *args, **kwargs):
        """Overloaded function.
        
        1. __init__(self: libaster.ElementaryVectorReal) -> None
        
        2. __init__(self: libaster.ElementaryVectorReal, arg0: str) -> None
        
        3. __init__(self: libaster.ElementaryVectorReal, arg0: PhysicalProblem) -> None
        """
    
    def assemble(self, arg0):
        pass
    
    def getVeass(self):
        pass

# class ElementaryVectorComplex in libaster

class ElementaryVectorComplex(BaseElementaryVector):
    pass
    
    # Method resolution order:
    #     ElementaryVectorComplex
    #     BaseElementaryVector
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __init__(self, *args, **kwargs):
        """Overloaded function.
        
        1. __init__(self: libaster.ElementaryVectorComplex) -> None
        
        2. __init__(self: libaster.ElementaryVectorComplex, arg0: str) -> None
        
        3. __init__(self: libaster.ElementaryVectorComplex, arg0: PhysicalProblem) -> None
        """
    
    def assemble(self, arg0):
        pass
    
    def getVeass(self):
        pass

# class ElementaryVectorDisplacementReal in libaster

class ElementaryVectorDisplacementReal(ElementaryVectorReal):
    pass
    
    # Method resolution order:
    #     ElementaryVectorDisplacementReal
    #     ElementaryVectorReal
    #     BaseElementaryVector
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __init__(self, *args, **kwargs):
        """Overloaded function.
        
        1. __init__(self: libaster.ElementaryVectorDisplacementReal) -> None
        
        2. __init__(self: libaster.ElementaryVectorDisplacementReal, arg0: str) -> None
        
        3. __init__(self: libaster.ElementaryVectorDisplacementReal, arg0: PhysicalProblem) -> None
        """

# class ElementaryVectorTemperatureReal in libaster

class ElementaryVectorTemperatureReal(ElementaryVectorReal):
    pass
    
    # Method resolution order:
    #     ElementaryVectorTemperatureReal
    #     ElementaryVectorReal
    #     BaseElementaryVector
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __init__(self, *args, **kwargs):
        """Overloaded function.
        
        1. __init__(self: libaster.ElementaryVectorTemperatureReal) -> None
        
        2. __init__(self: libaster.ElementaryVectorTemperatureReal, arg0: str) -> None
        
        3. __init__(self: libaster.ElementaryVectorTemperatureReal, arg0: PhysicalProblem) -> None
        """

# class ElementaryVectorPressureComplex in libaster

class ElementaryVectorPressureComplex(ElementaryVectorComplex):
    pass
    
    # Method resolution order:
    #     ElementaryVectorPressureComplex
    #     ElementaryVectorComplex
    #     BaseElementaryVector
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __init__(self, *args, **kwargs):
        """Overloaded function.
        
        1. __init__(self: libaster.ElementaryVectorPressureComplex) -> None
        
        2. __init__(self: libaster.ElementaryVectorPressureComplex, arg0: str) -> None
        
        3. __init__(self: libaster.ElementaryVectorPressureComplex, arg0: PhysicalProblem) -> None
        """

# class GeneralizedAssemblyMatrix in libaster

class GeneralizedAssemblyMatrix(DataStructure):
    pass
    
    # Method resolution order:
    #     GeneralizedAssemblyMatrix
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __init__(self,  *args, **kwargs):
        """Initialize self.  See help(type(self)) for accurate signature.
        """
    
    def getGeneralizedDOFNumbering(self):
        pass
    
    def getModalBasis(self):
        pass
    
    def setGeneralizedDOFNumbering(self, arg0):
        pass
    
    def setModalBasis(self, *args, **kwargs):
        """Overloaded function.
        
        1. setModalBasis(self: libaster.GeneralizedAssemblyMatrix, arg0: GeneralizedModeResult) -> bool
        
        2. setModalBasis(self: libaster.GeneralizedAssemblyMatrix, arg0: ModeResult) -> bool
        """

# class GeneralizedAssemblyMatrixReal in libaster

class GeneralizedAssemblyMatrixReal(GeneralizedAssemblyMatrix):
    pass
    
    # Method resolution order:
    #     GeneralizedAssemblyMatrixReal
    #     GeneralizedAssemblyMatrix
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __init__(self, *args, **kwargs):
        """Overloaded function.
        
        1. __init__(self: libaster.GeneralizedAssemblyMatrixReal) -> None
        
        2. __init__(self: libaster.GeneralizedAssemblyMatrixReal, arg0: str) -> None
        """

# class GeneralizedAssemblyMatrixComplex in libaster

class GeneralizedAssemblyMatrixComplex(GeneralizedAssemblyMatrix):
    pass
    
    # Method resolution order:
    #     GeneralizedAssemblyMatrixComplex
    #     GeneralizedAssemblyMatrix
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __init__(self, *args, **kwargs):
        """Overloaded function.
        
        1. __init__(self: libaster.GeneralizedAssemblyMatrixComplex) -> None
        
        2. __init__(self: libaster.GeneralizedAssemblyMatrixComplex, arg0: str) -> None
        """

# class GeneralizedAssemblyVector in libaster

class GeneralizedAssemblyVector(DataStructure):
    pass
    
    # Method resolution order:
    #     GeneralizedAssemblyVector
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __init__(self,  *args, **kwargs):
        """Initialize self.  See help(type(self)) for accurate signature.
        """

# class GeneralizedAssemblyVectorReal in libaster

class GeneralizedAssemblyVectorReal(GeneralizedAssemblyVector):
    pass
    
    # Method resolution order:
    #     GeneralizedAssemblyVectorReal
    #     GeneralizedAssemblyVector
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __init__(self, *args, **kwargs):
        """Overloaded function.
        
        1. __init__(self: libaster.GeneralizedAssemblyVectorReal) -> None
        
        2. __init__(self: libaster.GeneralizedAssemblyVectorReal, arg0: str) -> None
        """

# class GeneralizedAssemblyVectorComplex in libaster

class GeneralizedAssemblyVectorComplex(GeneralizedAssemblyVector):
    pass
    
    # Method resolution order:
    #     GeneralizedAssemblyVectorComplex
    #     GeneralizedAssemblyVector
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __init__(self, *args, **kwargs):
        """Overloaded function.
        
        1. __init__(self: libaster.GeneralizedAssemblyVectorComplex) -> None
        
        2. __init__(self: libaster.GeneralizedAssemblyVectorComplex, arg0: str) -> None
        """

# class InterspectralMatrix in libaster

class InterspectralMatrix(DataStructure):
    pass
    
    # Method resolution order:
    #     InterspectralMatrix
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __init__(self, *args, **kwargs):
        """Overloaded function.
        
        1. __init__(self: libaster.InterspectralMatrix) -> None
        
        2. __init__(self: libaster.InterspectralMatrix, arg0: str) -> None
        """

# class LinearSolver in libaster

class LinearSolver(DataStructure):
    pass
    
    # Method resolution order:
    #     LinearSolver
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __init__(self,  *args, **kwargs):
        """Initialize self.  See help(type(self)) for accurate signature.
        """
    
    def build(self):
        """build internal objects of the solver
        
        Returns:
             bool: True if the building is a success, else False
        """
    
    def deleteFactorizedMatrix(self):
        """delete the factorized matrix and its preconditionner if created.
        This is the case for Mumps and Petsc.
        
        Returns:
             bool: True if success, else False
        """
    
    def enableXfem(self):
        """Enable preconditionning for XFEM modeling.
        """
    
    def factorize(self, matrix):
        """Factorize the matrix.
        
        Arguments:
            matrix [BaseAssemblyMatrix] : matrix to factorize
        
        Returns:
            bool: True if factorization is a success, else False
        """
    
    def getMatrix(self):
        """return the factorized matrix
        
        Returns:
            BaseAssemblyMatrix: factorized matrix
        """
    
    def getPrecondMatrix(self):
        """return the preconditionning matrix
        
        Returns:
            BaseAssemblyMatrix: preconditionning matrix
        """
    
    def getSolverName(self):
        """Get the name of the solver used between 'MUMPS', 'PETSC', 'MULT_FRONT' and 'PETSC'
        
        Returns:
             str: name of the solver used
        """
    
    def setCataPath(self, path):
        """Set the path of the catalog that defines the solver keywords.
        It can be command name or a path as *code_aster.Cata.Commons.xxxx*.
        
        Arguments:
             path (str): command name or path of the catalog.
        """
    
    def setKeywords(self, arg0):
        pass
    
    def solve(self, *args, **kwargs):
        """Overloaded function.
        
        1. solve(self: libaster.LinearSolver, rhs: libaster.FieldOnNodesReal, dirichletBC: libaster.FieldOnNodesReal = None) -> libaster.FieldOnNodesReal
        
        2. solve(self: libaster.LinearSolver, rhs: libaster.FieldOnNodesComplex, dirichletBC: libaster.FieldOnNodesComplex = None) -> libaster.FieldOnNodesComplex
        """
    
    def supportParallelMesh(self):
        """tell if the solver is enable in HPC
        
        Returns:
             bool: True if the solver support ParallelMesh, else False
        """

# class MultFrontSolver in libaster

class MultFrontSolver(LinearSolver):
    pass
    
    # Method resolution order:
    #     MultFrontSolver
    #     LinearSolver
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __init__(self, *args, **kwargs):
        """Overloaded function.
        
        1. __init__(self: libaster.MultFrontSolver) -> None
        
        2. __init__(self: libaster.MultFrontSolver, arg0: str) -> None
        """

# class LdltSolver in libaster

class LdltSolver(LinearSolver):
    pass
    
    # Method resolution order:
    #     LdltSolver
    #     LinearSolver
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __init__(self, *args, **kwargs):
        """Overloaded function.
        
        1. __init__(self: libaster.LdltSolver) -> None
        
        2. __init__(self: libaster.LdltSolver, arg0: str) -> None
        """

# class MumpsSolver in libaster

class MumpsSolver(LinearSolver):
    pass
    
    # Method resolution order:
    #     MumpsSolver
    #     LinearSolver
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __init__(self, *args, **kwargs):
        """Overloaded function.
        
        1. __init__(self: libaster.MumpsSolver) -> None
        
        2. __init__(self: libaster.MumpsSolver, arg0: str) -> None
        """

# class PetscSolver in libaster

class PetscSolver(LinearSolver):
    pass
    
    # Method resolution order:
    #     PetscSolver
    #     LinearSolver
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __init__(self, *args, **kwargs):
        """Overloaded function.
        
        1. __init__(self: libaster.PetscSolver) -> None
        
        2. __init__(self: libaster.PetscSolver, arg0: str) -> None
        """

# class GcpcSolver in libaster

class GcpcSolver(LinearSolver):
    pass
    
    # Method resolution order:
    #     GcpcSolver
    #     LinearSolver
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __init__(self, *args, **kwargs):
        """Overloaded function.
        
        1. __init__(self: libaster.GcpcSolver) -> None
        
        2. __init__(self: libaster.GcpcSolver, arg0: str) -> None
        """

# class GenericModalBasis in libaster

class GenericModalBasis(DataStructure):
    pass
    
    # Method resolution order:
    #     GenericModalBasis
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __init__(self,  *args, **kwargs):
        """Initialize self.  See help(type(self)) for accurate signature.
        """

# class StandardModalBasis in libaster

class StandardModalBasis(GenericModalBasis):
    pass
    
    # Method resolution order:
    #     StandardModalBasis
    #     GenericModalBasis
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __init__(self, *args, **kwargs):
        """Overloaded function.
        
        1. __init__(self: libaster.StandardModalBasis) -> None
        
        2. __init__(self: libaster.StandardModalBasis, arg0: str) -> None
        """

# class RitzBasis in libaster

class RitzBasis(GenericModalBasis):
    pass
    
    # Method resolution order:
    #     RitzBasis
    #     GenericModalBasis
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __init__(self, *args, **kwargs):
        """Overloaded function.
        
        1. __init__(self: libaster.RitzBasis) -> None
        
        2. __init__(self: libaster.RitzBasis, arg0: str) -> None
        """

# class InterfaceType in libaster

class InterfaceType:
    """Members:
    
    MacNeal
    
    CraigBampton
    
    HarmonicCraigBampton
    
    None
    """
    
    # Method resolution order:
    #     InterfaceType
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __eq__(self, other):
        pass
    
    def __getstate__(self):
        pass
    
    def __hash__(self):
        pass
    
    def __index__(self):
        pass
    
    def __init__(self, value):
        pass
    
    def __int__(self):
        pass
    
    def __ne__(self, other):
        pass
    
    def __repr__(self):
        pass
    
    def __setstate__(self, state):
        pass
    
    def name(self):
        pass
    
    #----------------------------------------------------------------------
    # Data descriptors defined here:
    
    @property
    def __members__(self):
        pass
    
    @property
    def name(self):
        """name(self: handle) -> str
        """
    
    @property
    def value(self):
        pass
    
    #----------------------------------------------------------------------
    # Data and other attributes defined here:
    
    CraigBampton = 1
    
    HarmonicCraigBampton = 2
    
    MacNeal = 0

# class StructureInterface in libaster

class StructureInterface(DataStructure):
    pass
    
    # Method resolution order:
    #     StructureInterface
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __init__(self, *args, **kwargs):
        """Overloaded function.
        
        1. __init__(self: libaster.StructureInterface) -> None
        
        2. __init__(self: libaster.StructureInterface, arg0: str) -> None
        
        3. __init__(self: libaster.StructureInterface, arg0: libaster.DOFNumbering) -> None
        
        4. __init__(self: libaster.StructureInterface, arg0: str, arg1: libaster.DOFNumbering) -> None
        """

# class AcousticLoadComplex in libaster

class AcousticLoadComplex(DataStructure):
    pass
    
    # Method resolution order:
    #     AcousticLoadComplex
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __init__(self, *args, **kwargs):
        """Overloaded function.
        
        1. __init__(self: libaster.AcousticLoadComplex, arg0: Model) -> None
        
        2. __init__(self: libaster.AcousticLoadComplex, arg0: str, arg1: Model) -> None
        """
    
    def getFiniteElementDescriptor(self):
        pass
    
    def getMesh(self):
        pass
    
    def getModel(self):
        pass

# class DirichletBC in libaster

class DirichletBC(DataStructure):
    pass
    
    # Method resolution order:
    #     DirichletBC
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __init__(self,  *args, **kwargs):
        """Initialize self.  See help(type(self)) for accurate signature.
        """
    
    def build(self):
        pass
    
    def getModel(self):
        """Return the model
        
        Returns:
            ModelPtr: a pointer to the model
        """
    
    def getPhysics(self):
        """To know the physics supported by the model
        
        Returns:
            str: Mechanics or Thermal or Acoustic
        """

# class MechanicalDirichletBC in libaster

class MechanicalDirichletBC(DirichletBC):
    pass
    
    # Method resolution order:
    #     MechanicalDirichletBC
    #     DirichletBC
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __init__(self, *args, **kwargs):
        """Overloaded function.
        
        1. __init__(self: libaster.MechanicalDirichletBC, arg0: Model) -> None
        
        2. __init__(self: libaster.MechanicalDirichletBC, arg0: str, arg1: Model) -> None
        """
    
    def addBCOnCells(self, *args, **kwargs):
        """Overloaded function.
        
        1. addBCOnCells(self: libaster.MechanicalDirichletBC, arg0: PhysicalQuantityComponent, arg1: float, arg2: str) -> bool
        
        2. addBCOnCells(self: libaster.MechanicalDirichletBC, arg0: PhysicalQuantityComponent, arg1: float, arg2: List[str]) -> bool
        """
    
    def addBCOnNodes(self, *args, **kwargs):
        """Overloaded function.
        
        1. addBCOnNodes(self: libaster.MechanicalDirichletBC, arg0: PhysicalQuantityComponent, arg1: float, arg2: str) -> bool
        
        2. addBCOnNodes(self: libaster.MechanicalDirichletBC, arg0: PhysicalQuantityComponent, arg1: float, arg2: List[str]) -> bool
        """

# class ThermalDirichletBC in libaster

class ThermalDirichletBC(DirichletBC):
    pass
    
    # Method resolution order:
    #     ThermalDirichletBC
    #     DirichletBC
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __init__(self, *args, **kwargs):
        """Overloaded function.
        
        1. __init__(self: libaster.ThermalDirichletBC, arg0: Model) -> None
        
        2. __init__(self: libaster.ThermalDirichletBC, arg0: str, arg1: Model) -> None
        """
    
    def addBCOnCells(self, *args, **kwargs):
        """Overloaded function.
        
        1. addBCOnCells(self: libaster.ThermalDirichletBC, arg0: PhysicalQuantityComponent, arg1: float, arg2: str) -> bool
        
        2. addBCOnCells(self: libaster.ThermalDirichletBC, arg0: PhysicalQuantityComponent, arg1: float, arg2: List[str]) -> bool
        """
    
    def addBCOnNodes(self, *args, **kwargs):
        """Overloaded function.
        
        1. addBCOnNodes(self: libaster.ThermalDirichletBC, arg0: PhysicalQuantityComponent, arg1: float, arg2: str) -> bool
        
        2. addBCOnNodes(self: libaster.ThermalDirichletBC, arg0: PhysicalQuantityComponent, arg1: float, arg2: List[str]) -> bool
        
        3. addBCOnNodes(self: libaster.ThermalDirichletBC, arg0: PhysicalQuantityComponent, arg1: libaster.Function, arg2: List[str]) -> bool
        """

# class AcousticDirichletBC in libaster

class AcousticDirichletBC(DirichletBC):
    pass
    
    # Method resolution order:
    #     AcousticDirichletBC
    #     DirichletBC
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __init__(self, *args, **kwargs):
        """Overloaded function.
        
        1. __init__(self: libaster.AcousticDirichletBC, arg0: Model) -> None
        
        2. __init__(self: libaster.AcousticDirichletBC, arg0: str, arg1: Model) -> None
        """

# class MechanicalLoadReal in libaster

class MechanicalLoadReal(DataStructure):
    pass
    
    # Method resolution order:
    #     MechanicalLoadReal
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __init__(self, *args, **kwargs):
        """Overloaded function.
        
        1. __init__(self: libaster.MechanicalLoadReal, arg0: Model) -> None
        
        2. __init__(self: libaster.MechanicalLoadReal, arg0: str, arg1: Model) -> None
        """
    
    def getFiniteElementDescriptor(self):
        pass
    
    def getMesh(self):
        pass
    
    def getModel(self):
        pass
    
    def getTable(self, identifier):
        """Extract a Table from the datastructure.
        
        Arguments:
            identifier (str): Table identifier.
        
        Returns:
            Table: Table stored with the given identifier.
        """
    
    def hasLoadField(self, arg0):
        pass
    
    def updateValuePointers(self):
        pass

# class MechanicalLoadFunction in libaster

class MechanicalLoadFunction(DataStructure):
    pass
    
    # Method resolution order:
    #     MechanicalLoadFunction
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __init__(self, *args, **kwargs):
        """Overloaded function.
        
        1. __init__(self: libaster.MechanicalLoadFunction, arg0: Model) -> None
        
        2. __init__(self: libaster.MechanicalLoadFunction, arg0: str, arg1: Model) -> None
        """
    
    def getFiniteElementDescriptor(self):
        pass
    
    def getMesh(self):
        pass
    
    def getModel(self):
        pass
    
    def getTable(self, identifier):
        """Extract a Table from the datastructure.
        
        Arguments:
            identifier (str): Table identifier.
        
        Returns:
            Table: Table stored with the given identifier.
        """
    
    def hasLoadField(self, arg0):
        pass
    
    def updateValuePointers(self):
        pass

# class MechanicalLoadComplex in libaster

class MechanicalLoadComplex(DataStructure):
    pass
    
    # Method resolution order:
    #     MechanicalLoadComplex
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __init__(self, *args, **kwargs):
        """Overloaded function.
        
        1. __init__(self: libaster.MechanicalLoadComplex, arg0: Model) -> None
        
        2. __init__(self: libaster.MechanicalLoadComplex, arg0: str, arg1: Model) -> None
        """
    
    def getFiniteElementDescriptor(self):
        pass
    
    def getMesh(self):
        pass
    
    def getModel(self):
        pass
    
    def getTable(self, identifier):
        """Extract a Table from the datastructure.
        
        Arguments:
            identifier (str): Table identifier.
        
        Returns:
            Table: Table stored with the given identifier.
        """
    
    def hasLoadField(self, arg0):
        pass
    
    def updateValuePointers(self):
        pass

# class Loads in libaster

class Loads:
    """Members:
    
    NodalForce
    
    ForceOnEdge
    
    ForceOnFace
    
    LineicForce
    
    InternalForce
    
    ForceOnBeam
    
    ForceOnShell
    
    PressureOnPipe
    
    ImposedDoF
    
    DistributedPressure
    
    ImpedanceOnFace
    
    NormalSpeedOnFace
    
    WavePressureOnFace
    
    THMFlux
    """
    
    # Method resolution order:
    #     Loads
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __eq__(self, other):
        pass
    
    def __getstate__(self):
        pass
    
    def __hash__(self):
        pass
    
    def __index__(self):
        pass
    
    def __init__(self, value):
        pass
    
    def __int__(self):
        pass
    
    def __ne__(self, other):
        pass
    
    def __repr__(self):
        pass
    
    def __setstate__(self, state):
        pass
    
    def name(self):
        pass
    
    #----------------------------------------------------------------------
    # Data descriptors defined here:
    
    @property
    def __members__(self):
        pass
    
    @property
    def name(self):
        """name(self: handle) -> str
        """
    
    @property
    def value(self):
        pass
    
    #----------------------------------------------------------------------
    # Data and other attributes defined here:
    
    DistributedPressure = 9
    
    ForceOnBeam = 5
    
    ForceOnEdge = 1
    
    ForceOnFace = 2
    
    ForceOnShell = 6
    
    ImpedanceOnFace = 10
    
    ImposedDoF = 8
    
    InternalForce = 4
    
    LineicForce = 3
    
    NodalForce = 0
    
    NormalSpeedOnFace = 11
    
    PressureOnPipe = 7
    
    THMFlux = 13
    
    WavePressureOnFace = 12

# class NodalForceReal in libaster

class NodalForceReal(MechanicalLoadReal):
    pass
    
    # Method resolution order:
    #     NodalForceReal
    #     MechanicalLoadReal
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __init__(self, *args, **kwargs):
        """Overloaded function.
        
        1. __init__(self: libaster.NodalForceReal, arg0: Model) -> None
        
        2. __init__(self: libaster.NodalForceReal, arg0: str, arg1: Model) -> None
        """
    
    def build(self):
        pass
    
    def setValue(self, arg0, arg1):
        pass

# class NodalStructuralForceReal in libaster

class NodalStructuralForceReal(MechanicalLoadReal):
    pass
    
    # Method resolution order:
    #     NodalStructuralForceReal
    #     MechanicalLoadReal
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __init__(self, *args, **kwargs):
        """Overloaded function.
        
        1. __init__(self: libaster.NodalStructuralForceReal, arg0: Model) -> None
        
        2. __init__(self: libaster.NodalStructuralForceReal, arg0: str, arg1: Model) -> None
        """
    
    def build(self):
        pass
    
    def setValue(self, arg0, arg1):
        pass

# class ForceOnFaceReal in libaster

class ForceOnFaceReal(MechanicalLoadReal):
    pass
    
    # Method resolution order:
    #     ForceOnFaceReal
    #     MechanicalLoadReal
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __init__(self, *args, **kwargs):
        """Overloaded function.
        
        1. __init__(self: libaster.ForceOnFaceReal, arg0: Model) -> None
        
        2. __init__(self: libaster.ForceOnFaceReal, arg0: str, arg1: Model) -> None
        """
    
    def build(self):
        pass
    
    def setValue(self, arg0, arg1):
        pass

# class ForceOnEdgeReal in libaster

class ForceOnEdgeReal(MechanicalLoadReal):
    pass
    
    # Method resolution order:
    #     ForceOnEdgeReal
    #     MechanicalLoadReal
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __init__(self, *args, **kwargs):
        """Overloaded function.
        
        1. __init__(self: libaster.ForceOnEdgeReal, arg0: Model) -> None
        
        2. __init__(self: libaster.ForceOnEdgeReal, arg0: str, arg1: Model) -> None
        """
    
    def build(self):
        pass
    
    def setValue(self, arg0, arg1):
        pass

# class StructuralForceOnEdgeReal in libaster

class StructuralForceOnEdgeReal(MechanicalLoadReal):
    pass
    
    # Method resolution order:
    #     StructuralForceOnEdgeReal
    #     MechanicalLoadReal
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __init__(self, *args, **kwargs):
        """Overloaded function.
        
        1. __init__(self: libaster.StructuralForceOnEdgeReal, arg0: Model) -> None
        
        2. __init__(self: libaster.StructuralForceOnEdgeReal, arg0: str, arg1: Model) -> None
        """
    
    def build(self):
        pass
    
    def setValue(self, arg0, arg1):
        pass

# class LineicForceReal in libaster

class LineicForceReal(MechanicalLoadReal):
    pass
    
    # Method resolution order:
    #     LineicForceReal
    #     MechanicalLoadReal
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __init__(self, *args, **kwargs):
        """Overloaded function.
        
        1. __init__(self: libaster.LineicForceReal, arg0: Model) -> None
        
        2. __init__(self: libaster.LineicForceReal, arg0: str, arg1: Model) -> None
        """
    
    def build(self):
        pass
    
    def setValue(self, arg0, arg1):
        pass

# class InternalForceReal in libaster

class InternalForceReal(MechanicalLoadReal):
    pass
    
    # Method resolution order:
    #     InternalForceReal
    #     MechanicalLoadReal
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __init__(self, *args, **kwargs):
        """Overloaded function.
        
        1. __init__(self: libaster.InternalForceReal, arg0: Model) -> None
        
        2. __init__(self: libaster.InternalForceReal, arg0: str, arg1: Model) -> None
        """
    
    def build(self):
        pass
    
    def setValue(self, arg0, arg1):
        pass

# class StructuralForceOnBeamReal in libaster

class StructuralForceOnBeamReal(MechanicalLoadReal):
    pass
    
    # Method resolution order:
    #     StructuralForceOnBeamReal
    #     MechanicalLoadReal
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __init__(self, *args, **kwargs):
        """Overloaded function.
        
        1. __init__(self: libaster.StructuralForceOnBeamReal, arg0: Model) -> None
        
        2. __init__(self: libaster.StructuralForceOnBeamReal, arg0: str, arg1: Model) -> None
        """
    
    def build(self):
        pass
    
    def setValue(self, arg0, arg1):
        pass

# class LocalForceOnBeamReal in libaster

class LocalForceOnBeamReal(MechanicalLoadReal):
    pass
    
    # Method resolution order:
    #     LocalForceOnBeamReal
    #     MechanicalLoadReal
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __init__(self, *args, **kwargs):
        """Overloaded function.
        
        1. __init__(self: libaster.LocalForceOnBeamReal, arg0: Model) -> None
        
        2. __init__(self: libaster.LocalForceOnBeamReal, arg0: str, arg1: Model) -> None
        """
    
    def build(self):
        pass
    
    def setValue(self, arg0, arg1):
        pass

# class StructuralForceOnShellReal in libaster

class StructuralForceOnShellReal(MechanicalLoadReal):
    pass
    
    # Method resolution order:
    #     StructuralForceOnShellReal
    #     MechanicalLoadReal
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __init__(self, *args, **kwargs):
        """Overloaded function.
        
        1. __init__(self: libaster.StructuralForceOnShellReal, arg0: Model) -> None
        
        2. __init__(self: libaster.StructuralForceOnShellReal, arg0: str, arg1: Model) -> None
        """
    
    def build(self):
        pass
    
    def setValue(self, arg0, arg1):
        pass

# class LocalForceOnShellReal in libaster

class LocalForceOnShellReal(MechanicalLoadReal):
    pass
    
    # Method resolution order:
    #     LocalForceOnShellReal
    #     MechanicalLoadReal
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __init__(self, *args, **kwargs):
        """Overloaded function.
        
        1. __init__(self: libaster.LocalForceOnShellReal, arg0: Model) -> None
        
        2. __init__(self: libaster.LocalForceOnShellReal, arg0: str, arg1: Model) -> None
        """
    
    def build(self):
        pass
    
    def setValue(self, arg0, arg1):
        pass

# class PressureOnShellReal in libaster

class PressureOnShellReal(MechanicalLoadReal):
    pass
    
    # Method resolution order:
    #     PressureOnShellReal
    #     MechanicalLoadReal
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __init__(self, *args, **kwargs):
        """Overloaded function.
        
        1. __init__(self: libaster.PressureOnShellReal, arg0: Model) -> None
        
        2. __init__(self: libaster.PressureOnShellReal, arg0: str, arg1: Model) -> None
        """
    
    def build(self):
        pass
    
    def setValue(self, arg0, arg1):
        pass

# class PressureOnPipeReal in libaster

class PressureOnPipeReal(MechanicalLoadReal):
    pass
    
    # Method resolution order:
    #     PressureOnPipeReal
    #     MechanicalLoadReal
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __init__(self, *args, **kwargs):
        """Overloaded function.
        
        1. __init__(self: libaster.PressureOnPipeReal, arg0: Model) -> None
        
        2. __init__(self: libaster.PressureOnPipeReal, arg0: str, arg1: Model) -> None
        """
    
    def build(self):
        pass
    
    def setValue(self, arg0, arg1):
        pass

# class ImposedDisplacementReal in libaster

class ImposedDisplacementReal(MechanicalLoadReal):
    pass
    
    # Method resolution order:
    #     ImposedDisplacementReal
    #     MechanicalLoadReal
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __init__(self, *args, **kwargs):
        """Overloaded function.
        
        1. __init__(self: libaster.ImposedDisplacementReal, arg0: Model) -> None
        
        2. __init__(self: libaster.ImposedDisplacementReal, arg0: str, arg1: Model) -> None
        """
    
    def build(self):
        pass
    
    def setValue(self, arg0, arg1):
        pass

# class ImposedPressureReal in libaster

class ImposedPressureReal(MechanicalLoadReal):
    pass
    
    # Method resolution order:
    #     ImposedPressureReal
    #     MechanicalLoadReal
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __init__(self, *args, **kwargs):
        """Overloaded function.
        
        1. __init__(self: libaster.ImposedPressureReal, arg0: Model) -> None
        
        2. __init__(self: libaster.ImposedPressureReal, arg0: str, arg1: Model) -> None
        """
    
    def build(self):
        pass
    
    def setValue(self, arg0, arg1):
        pass

# class DistributedPressureReal in libaster

class DistributedPressureReal(MechanicalLoadReal):
    pass
    
    # Method resolution order:
    #     DistributedPressureReal
    #     MechanicalLoadReal
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __init__(self, *args, **kwargs):
        """Overloaded function.
        
        1. __init__(self: libaster.DistributedPressureReal, arg0: Model) -> None
        
        2. __init__(self: libaster.DistributedPressureReal, arg0: str, arg1: Model) -> None
        """
    
    def build(self):
        pass
    
    def setValue(self, arg0, arg1):
        pass

# class ImpedanceOnFaceReal in libaster

class ImpedanceOnFaceReal(MechanicalLoadReal):
    pass
    
    # Method resolution order:
    #     ImpedanceOnFaceReal
    #     MechanicalLoadReal
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __init__(self, *args, **kwargs):
        """Overloaded function.
        
        1. __init__(self: libaster.ImpedanceOnFaceReal, arg0: Model) -> None
        
        2. __init__(self: libaster.ImpedanceOnFaceReal, arg0: str, arg1: Model) -> None
        """
    
    def build(self):
        pass
    
    def setValue(self, arg0, arg1):
        pass

# class NormalSpeedOnFaceReal in libaster

class NormalSpeedOnFaceReal(MechanicalLoadReal):
    pass
    
    # Method resolution order:
    #     NormalSpeedOnFaceReal
    #     MechanicalLoadReal
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __init__(self, *args, **kwargs):
        """Overloaded function.
        
        1. __init__(self: libaster.NormalSpeedOnFaceReal, arg0: Model) -> None
        
        2. __init__(self: libaster.NormalSpeedOnFaceReal, arg0: str, arg1: Model) -> None
        """
    
    def build(self):
        pass
    
    def setValue(self, arg0, arg1):
        pass

# class WavePressureOnFaceReal in libaster

class WavePressureOnFaceReal(MechanicalLoadReal):
    pass
    
    # Method resolution order:
    #     WavePressureOnFaceReal
    #     MechanicalLoadReal
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __init__(self, *args, **kwargs):
        """Overloaded function.
        
        1. __init__(self: libaster.WavePressureOnFaceReal, arg0: Model) -> None
        
        2. __init__(self: libaster.WavePressureOnFaceReal, arg0: str, arg1: Model) -> None
        """
    
    def build(self):
        pass
    
    def setValue(self, arg0, arg1):
        pass

# class DistributedHeatFluxReal in libaster

class DistributedHeatFluxReal(MechanicalLoadReal):
    pass
    
    # Method resolution order:
    #     DistributedHeatFluxReal
    #     MechanicalLoadReal
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __init__(self, *args, **kwargs):
        """Overloaded function.
        
        1. __init__(self: libaster.DistributedHeatFluxReal, arg0: Model) -> None
        
        2. __init__(self: libaster.DistributedHeatFluxReal, arg0: str, arg1: Model) -> None
        """
    
    def build(self):
        pass
    
    def setValue(self, arg0, arg1):
        pass

# class DistributedHydraulicFluxReal in libaster

class DistributedHydraulicFluxReal(MechanicalLoadReal):
    pass
    
    # Method resolution order:
    #     DistributedHydraulicFluxReal
    #     MechanicalLoadReal
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __init__(self, *args, **kwargs):
        """Overloaded function.
        
        1. __init__(self: libaster.DistributedHydraulicFluxReal, arg0: Model) -> None
        
        2. __init__(self: libaster.DistributedHydraulicFluxReal, arg0: str, arg1: Model) -> None
        """
    
    def build(self):
        pass
    
    def setValue(self, arg0, arg1):
        pass

# class PhysicalQuantityComponent in libaster

class PhysicalQuantityComponent:
    """Members:
    
    Dx
    
    Dy
    
    Dz
    
    Drx
    
    Dry
    
    Drz
    
    Temp
    
    MiddleTemp
    
    Pres
    
    Fx
    
    Fy
    
    Fz
    
    Mx
    
    My
    
    Mz
    
    N
    
    Vy
    
    Vz
    
    Mt
    
    Mfy
    
    Mfz
    
    F1
    
    F2
    
    F3
    
    Mf1
    
    Mf2
    
    Impe
    
    Vnor
    
    Flun
    
    FlunHydr1
    
    FlunHydr2
    """
    
    # Method resolution order:
    #     PhysicalQuantityComponent
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __eq__(self, other):
        pass
    
    def __getstate__(self):
        pass
    
    def __hash__(self):
        pass
    
    def __index__(self):
        pass
    
    def __init__(self, value):
        pass
    
    def __int__(self):
        pass
    
    def __ne__(self, other):
        pass
    
    def __repr__(self):
        pass
    
    def __setstate__(self, state):
        pass
    
    def name(self):
        pass
    
    #----------------------------------------------------------------------
    # Data descriptors defined here:
    
    @property
    def __members__(self):
        pass
    
    @property
    def name(self):
        """name(self: handle) -> str
        """
    
    @property
    def value(self):
        pass
    
    #----------------------------------------------------------------------
    # Data and other attributes defined here:
    
    Drx = 3
    
    Dry = 4
    
    Drz = 5
    
    Dx = 0
    
    Dy = 1
    
    Dz = 2
    
    F1 = 21
    
    F2 = 22
    
    F3 = 23
    
    Flun = 28
    
    FlunHydr1 = 29
    
    FlunHydr2 = 30
    
    Fx = 9
    
    Fy = 10
    
    Fz = 11
    
    Impe = 26
    
    Mf1 = 24
    
    Mf2 = 25
    
    Mfy = 19
    
    Mfz = 20
    
    MiddleTemp = 7
    
    Mt = 18
    
    Mx = 12
    
    My = 13
    
    Mz = 14
    
    N = 15
    
    Pres = 8
    
    Temp = 6
    
    Vnor = 27
    
    Vy = 16
    
    Vz = 17

# class ForceReal in libaster

class ForceReal:
    pass
    
    # Method resolution order:
    #     ForceReal
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __init__(self):
        pass
    
    def debugPrint(self):
        pass
    
    def setValue(self, arg0, arg1):
        pass

# class StructuralForceReal in libaster

class StructuralForceReal:
    pass
    
    # Method resolution order:
    #     StructuralForceReal
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __init__(self):
        pass
    
    def debugPrint(self):
        pass
    
    def setValue(self, arg0, arg1):
        pass

# class LocalBeamForceReal in libaster

class LocalBeamForceReal:
    pass
    
    # Method resolution order:
    #     LocalBeamForceReal
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __init__(self):
        pass
    
    def debugPrint(self):
        pass
    
    def setValue(self, arg0, arg1):
        pass

# class LocalShellForceReal in libaster

class LocalShellForceReal:
    pass
    
    # Method resolution order:
    #     LocalShellForceReal
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __init__(self):
        pass
    
    def debugPrint(self):
        pass
    
    def setValue(self, arg0, arg1):
        pass

# class DisplacementReal in libaster

class DisplacementReal:
    pass
    
    # Method resolution order:
    #     DisplacementReal
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __init__(self):
        pass
    
    def debugPrint(self):
        pass
    
    def setValue(self, arg0, arg1):
        pass

# class PressureReal in libaster

class PressureReal:
    pass
    
    # Method resolution order:
    #     PressureReal
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __init__(self):
        pass
    
    def debugPrint(self):
        pass
    
    def setValue(self, arg0, arg1):
        pass

# class ImpedanceReal in libaster

class ImpedanceReal:
    pass
    
    # Method resolution order:
    #     ImpedanceReal
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __init__(self):
        pass
    
    def debugPrint(self):
        pass
    
    def setValue(self, arg0, arg1):
        pass

# class NormalSpeedReal in libaster

class NormalSpeedReal:
    pass
    
    # Method resolution order:
    #     NormalSpeedReal
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __init__(self):
        pass
    
    def debugPrint(self):
        pass
    
    def setValue(self, arg0, arg1):
        pass

# class HeatFluxReal in libaster

class HeatFluxReal:
    pass
    
    # Method resolution order:
    #     HeatFluxReal
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __init__(self):
        pass
    
    def debugPrint(self):
        pass
    
    def setValue(self, arg0, arg1):
        pass

# class HydraulicFluxReal in libaster

class HydraulicFluxReal:
    pass
    
    # Method resolution order:
    #     HydraulicFluxReal
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __init__(self):
        pass
    
    def debugPrint(self):
        pass
    
    def setValue(self, arg0, arg1):
        pass

# class ThermalLoadReal in libaster

class ThermalLoadReal(DataStructure):
    pass
    
    # Method resolution order:
    #     ThermalLoadReal
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __init__(self, *args, **kwargs):
        """Overloaded function.
        
        1. __init__(self: libaster.ThermalLoadReal, arg0: Model) -> None
        
        2. __init__(self: libaster.ThermalLoadReal, arg0: str, arg1: Model) -> None
        """
    
    def getFiniteElementDescriptor(self):
        pass
    
    def getMesh(self):
        pass
    
    def getModel(self):
        pass
    
    def hasLoadField(self, arg0):
        """Return true if the wanted field exists
        
        Arguments:
            str: name of the load field
        
        Returns:
            bool: field exists
        """
    
    def hasLoadResult(self):
        """Return true if the LoadResult structure exists
        
        Returns:
            bool: field exists
        """

# class ThermalLoadFunction in libaster

class ThermalLoadFunction(DataStructure):
    pass
    
    # Method resolution order:
    #     ThermalLoadFunction
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __init__(self, *args, **kwargs):
        """Overloaded function.
        
        1. __init__(self: libaster.ThermalLoadFunction, arg0: Model) -> None
        
        2. __init__(self: libaster.ThermalLoadFunction, arg0: str, arg1: Model) -> None
        """
    
    def getFiniteElementDescriptor(self):
        pass
    
    def getMesh(self):
        pass
    
    def getModel(self):
        pass
    
    def hasLoadField(self, arg0):
        """Return true if the wanted field exists
        
        Arguments:
            str: name of the load field
        
        Returns:
            bool: field exists
        """
    
    def hasLoadResult(self):
        """Return true if the LoadResult structure exists
        
        Returns:
            bool: field exists
        """

# class BehaviourDefinition in libaster

class BehaviourDefinition(DataStructure):
    pass
    
    # Method resolution order:
    #     BehaviourDefinition
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __init__(self, *args, **kwargs):
        """Overloaded function.
        
        1. __init__(self: libaster.BehaviourDefinition) -> None
        
        2. __init__(self: libaster.BehaviourDefinition, arg0: str) -> None
        """

# class Material in libaster

class Material(DataStructure):
    pass
    
    # Method resolution order:
    #     Material
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __init__(self, *args, **kwargs):
        """Overloaded function.
        
        1. __init__(self: libaster.Material) -> None
        
        2. __init__(self: libaster.Material, arg0: str) -> None
        
        3. __init__(self: libaster.Material, arg0: libaster.Material) -> None
        """
    
    def getFunction(self, materialName, propertyName):
        """Return the value of a property stored as a function.
        
        Raise an exception if the property does not exist.
        
        Arguments:
            materialName (str): Material name (without "_FO").
            propertyName (str): Property name.
        
        Returns:
            *Function*: Function object, *None* if the property does not exist or is not a function.
        """
    
    def getMaterialNames(self):
        """Return the list of the material names.
        
        Returns:
            list[str]: List of material names (without "_FO")
        """
    
    def getValueComplex(self, materialName, propertyName):
        """Return the value of a property stored as a complex.
        
        Raise an exception if the property does not exist.
        
        Arguments:
            materialName (str): Material name (without "_FO").
            propertyName (str): Property name.
        
        Returns:
            complex: Property value.
        """
    
    def getValueReal(self, materialName, propertyName):
        """Return the value of a property stored as a real.
        
        Raise an exception if the property does not exist.
        
        Arguments:
            materialName (str): Material name (without "_FO").
            propertyName (str): Property name.
        
        Returns:
            float: Property value.
        """
    
    def size(self):
        """Return the number of material names.
        
        Returns:
            int: Number of material names.
        """

# class PartOfMaterialField in libaster

class PartOfMaterialField:
    pass
    
    # Method resolution order:
    #     PartOfMaterialField
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __init__(self, *args, **kwargs):
        """Overloaded function.
        
        1. __init__(self: libaster.PartOfMaterialField) -> None
        
        2. __init__(self: libaster.PartOfMaterialField, arg0: List[libaster.Material], arg1: libaster.MeshEntity) -> None
        """
    
    def getMeshEntity(self):
        pass
    
    def getVectorOfMaterial(self):
        pass

# class MaterialField in libaster

class MaterialField(DataStructure):
    pass
    
    # Method resolution order:
    #     MaterialField
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __init__(self, *args, **kwargs):
        """Overloaded function.
        
        1. __init__(self: libaster.MaterialField, arg0: libaster.Mesh) -> None
        
        2. __init__(self: libaster.MaterialField, arg0: Skeleton) -> None
        
        3. __init__(self: libaster.MaterialField, arg0: str, arg1: libaster.Mesh) -> None
        
        4. __init__(self: libaster.MaterialField, arg0: ParallelMesh) -> None
        
        5. __init__(self: libaster.MaterialField, arg0: str, arg1: ParallelMesh) -> None
        """
    
    def addBehaviourOnGroupOfCells(self, behaviour, nameOfGroups):
        """Add behaviour (from DEFI_COMPOR) on group of cells
        
        Arguments:
            behaviour (BehaviourDefinition): Behaviour (from DEFI_COMPOR)
            nameOfGroups (list(str)) : list of names of groups of cells
        """
    
    def addBehaviourOnMesh(self, behaviour):
        """Add behaviour (from DEFI_COMPOR) on mesh
        
        Arguments:
            behaviour (BehaviourDefinition): Behaviour (from DEFI_COMPOR)
        """
    
    def addExternalStateVariable(self, exteVari):
        """Add external state variable in material field
        
        Arguments:
            exteVari (ExternalStateVariablePtr): external state variable
        """
    
    def addMaterialOnGroupOfCells(self, material, nameOfGroups):
        """Add a material properties on list of groups of cells
        
        Arguments:
            material (Material): material properties
            nameOfGroups (list(str)) : list of names of groups of cells
        """
    
    def addMaterialOnMesh(self, material):
        """Add material properties on mesh
        
        Arguments:
            material (Material): material properties
        """
    
    def addMultipleMaterialOnGroupOfCells(self, material, nameOfGroups):
        """Add a vector of multiple material properties on group of cells
        
        Arguments:
            material (list(Material)): list of material properties
            nameOfGroups (list(str)) : list of names of groups of cells
        """
    
    def addMultipleMaterialOnMesh(self, material):
        """Add a vector of multiple material properties on mesh
        
        Arguments:
            material (list(Material)): list of material properties
        """
    
    def build(self):
        """Build material field
        """
    
    def getMesh(self):
        """Get mesh of material field
        
        Returns:
            BaseMesh: mesh
        """
    
    def getVectorOfMaterial(self):
        """Get vector of all the material properties on the material field
        
        Returns:
            list(Material): vector of material properties
        """
    
    def getVectorOfPartOfMaterialField(self):
        """Get vector of all the material properties with mesh entities on the material field
        
        Returns:
            list(PartOfMaterial): vector of material properties with mesh entities
        """
    
    def hasExternalStateVariable(self, *args, **kwargs):
        """Overloaded function.
        
        1. hasExternalStateVariable(self: libaster.MaterialField, exteVariIden: externVarEnumInt) -> bool
        
        
                    Detects the presence of an external state variable
        
                    Arguments:
                        exteVariIden (externVarEnumInt or str): identifier for external state variable
        
                    Returns:
                        bool: True if has external state variables
                    
        
        2. hasExternalStateVariable(self: libaster.MaterialField, exteVariIden: str) -> bool
        
        
                    Detects the presence of an external state variable
        
                    Returns:
                        bool: True if has external state variables
                    
        
        3. hasExternalStateVariable(self: libaster.MaterialField) -> bool
        
        
                    Detects the presence of any external state variable
        
                    Returns:
                        bool: True if has external state variables
        """
    
    def hasExternalStateVariableForLoad(self):
        """Detects the presence of an external state variable for loads
        
        Returns:
            bool: True if has external state variables for loads
        """
    
    def hasExternalStateVariableWithReference(self):
        """Detects the presence of an external state variable with reference value
        
        Returns:
            bool: True if has external state variables with reference value
        """
    
    def setModel(self, model):
        """Set model of the material field
        
        Arguments:
            model (Model): model
        """
    
    def update(self):
        """Update material field
        """

# class Grid in libaster

class Grid(Mesh):
    pass
    
    # Method resolution order:
    #     Grid
    #     Mesh
    #     BaseMesh
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __init__(self, *args, **kwargs):
        """Overloaded function.
        
        1. __init__(self: libaster.Grid) -> None
        
        2. __init__(self: libaster.Grid, arg0: str) -> None
        """

# class MeshesMapping in libaster

class MeshesMapping(DataStructure):
    pass
    
    # Method resolution order:
    #     MeshesMapping
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __init__(self, *args, **kwargs):
        """Overloaded function.
        
        1. __init__(self: libaster.MeshesMapping) -> None
        
        2. __init__(self: libaster.MeshesMapping, arg0: str) -> None
        """
    
    def setFirstMesh(self, arg0):
        pass

# class Skeleton in libaster

class Skeleton(BaseMesh):
    pass
    
    # Method resolution order:
    #     Skeleton
    #     BaseMesh
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __init__(self, *args, **kwargs):
        """Overloaded function.
        
        1. __init__(self: libaster.Skeleton) -> None
        
        2. __init__(self: libaster.Skeleton, arg0: str) -> None
        """

# class DynamicMacroElement in libaster

class DynamicMacroElement(DataStructure):
    pass
    
    # Method resolution order:
    #     DynamicMacroElement
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __init__(self, *args, **kwargs):
        """Overloaded function.
        
        1. __init__(self: libaster.DynamicMacroElement) -> None
        
        2. __init__(self: libaster.DynamicMacroElement, arg0: str) -> None
        """
    
    def getDOFNumbering(self):
        pass
    
    def getDampingMatrix(self):
        pass
    
    def getImpedanceDampingMatrix(self):
        pass
    
    def getImpedanceMassMatrix(self):
        pass
    
    def getImpedanceMatrix(self):
        pass
    
    def getImpedanceStiffnessMatrix(self):
        pass
    
    def getMassMatrix(self):
        pass
    
    def getMechanicalMode(self):
        pass
    
    def getNumberOfNodes(self):
        pass
    
    def getStiffnessMatrixComplex(self):
        pass
    
    def getStiffnessMatrixReal(self):
        pass
    
    def setDampingMatrix(self, arg0):
        pass
    
    def setImpedanceDampingMatrix(self, arg0):
        pass
    
    def setImpedanceMassMatrix(self, arg0):
        pass
    
    def setImpedanceMatrix(self, arg0):
        pass
    
    def setImpedanceStiffnessMatrix(self, arg0):
        pass
    
    def setMassMatrix(self, arg0):
        pass
    
    def setMechanicalMode(self, arg0):
        pass
    
    def setStiffnessMatrix(self, *args, **kwargs):
        """Overloaded function.
        
        1. setStiffnessMatrix(self: libaster.DynamicMacroElement, arg0: libaster.AssemblyMatrixDisplacementComplex) -> bool
        
        2. setStiffnessMatrix(self: libaster.DynamicMacroElement, arg0: libaster.AssemblyMatrixDisplacementReal) -> bool
        """

# class StaticMacroElement in libaster

class StaticMacroElement(DataStructure):
    pass
    
    # Method resolution order:
    #     StaticMacroElement
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __init__(self, *args, **kwargs):
        """Overloaded function.
        
        1. __init__(self: libaster.StaticMacroElement) -> None
        
        2. __init__(self: libaster.StaticMacroElement, arg0: str) -> None
        """

# class CrackShape in libaster

class CrackShape:
    pass
    
    # Method resolution order:
    #     CrackShape
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __init__(self):
        pass
    
    def getCenter(self):
        pass
    
    def getCrackSide(self):
        pass
    
    def getEndPoint(self):
        pass
    
    def getFilletRadius(self):
        pass
    
    def getHalfLength(self):
        pass
    
    def getNormal(self):
        pass
    
    def getSemiMajorAxis(self):
        pass
    
    def getSemiMinorAxis(self):
        pass
    
    def getShape(self):
        pass
    
    def getShapeName(self):
        pass
    
    def getStartingPoint(self):
        pass
    
    def getTangent(self):
        pass
    
    def getVectX(self):
        pass
    
    def getVectY(self):
        pass
    
    def setCylinderCrackShape(self, arg0, arg1, arg2, arg3, arg4):
        pass
    
    def setEllipseCrackShape(self, arg0, arg1, arg2, arg3, arg4, arg5):
        pass
    
    def setHalfLineCrackShape(self, arg0, arg1):
        pass
    
    def setHalfPlaneCrackShape(self, arg0, arg1, arg2):
        pass
    
    def setLineCrackShape(self, arg0, arg1):
        pass
    
    def setNotchCrackShape(self, arg0, arg1, arg2, arg3, arg4):
        pass
    
    def setSegmentCrackShape(self, arg0, arg1):
        pass
    
    def setSquareCrackShape(self, arg0, arg1, arg2, arg3, arg4, arg5, arg6):
        pass

# class CrackTip in libaster

class CrackTip(DataStructure):
    pass
    
    # Method resolution order:
    #     CrackTip
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __init__(self, *args, **kwargs):
        """Overloaded function.
        
        1. __init__(self: libaster.CrackTip) -> None
        
        2. __init__(self: libaster.CrackTip, arg0: str) -> None
        """

# class Crack in libaster

class Crack(DataStructure):
    pass
    
    # Method resolution order:
    #     Crack
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __init__(self, *args, **kwargs):
        """Overloaded function.
        
        1. __init__(self: libaster.Crack) -> None
        
        2. __init__(self: libaster.Crack, arg0: str) -> None
        """
    
    def getCrackTipCellsType(self):
        pass
    
    def getLowerLipGroupName(self):
        pass
    
    def getUpperLipGroupName(self):
        pass

# class GeneralizedModel in libaster

class GeneralizedModel(DataStructure):
    pass
    
    # Method resolution order:
    #     GeneralizedModel
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __init__(self, *args, **kwargs):
        """Overloaded function.
        
        1. __init__(self: libaster.GeneralizedModel) -> None
        
        2. __init__(self: libaster.GeneralizedModel, arg0: str) -> None
        """
    
    def addDynamicMacroElement(self, arg0, arg1):
        pass
    
    def getDynamicMacroElementFromName(self, arg0):
        pass

# class ModelSplitingMethod in libaster

class ModelSplitingMethod:
    """Members:
    
    Centralized
    
    SubDomain
    
    GroupOfCells
    """
    
    # Method resolution order:
    #     ModelSplitingMethod
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __eq__(self, other):
        pass
    
    def __getstate__(self):
        pass
    
    def __hash__(self):
        pass
    
    def __index__(self):
        pass
    
    def __init__(self, value):
        pass
    
    def __int__(self):
        pass
    
    def __ne__(self, other):
        pass
    
    def __repr__(self):
        pass
    
    def __setstate__(self, state):
        pass
    
    def name(self):
        pass
    
    #----------------------------------------------------------------------
    # Data descriptors defined here:
    
    @property
    def __members__(self):
        pass
    
    @property
    def name(self):
        """name(self: handle) -> str
        """
    
    @property
    def value(self):
        pass
    
    #----------------------------------------------------------------------
    # Data and other attributes defined here:
    
    Centralized = 0
    
    GroupOfCells = 2
    
    SubDomain = 1

# class GraphPartitioner in libaster

class GraphPartitioner:
    """Members:
    
    Scotch
    
    Metis
    """
    
    # Method resolution order:
    #     GraphPartitioner
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __eq__(self, other):
        pass
    
    def __getstate__(self):
        pass
    
    def __hash__(self):
        pass
    
    def __index__(self):
        pass
    
    def __init__(self, value):
        pass
    
    def __int__(self):
        pass
    
    def __ne__(self, other):
        pass
    
    def __repr__(self):
        pass
    
    def __setstate__(self, state):
        pass
    
    def name(self):
        pass
    
    #----------------------------------------------------------------------
    # Data descriptors defined here:
    
    @property
    def __members__(self):
        pass
    
    @property
    def name(self):
        """name(self: handle) -> str
        """
    
    @property
    def value(self):
        pass
    
    #----------------------------------------------------------------------
    # Data and other attributes defined here:
    
    Metis = 1
    
    Scotch = 0

# class Model in libaster

class Model(DataStructure):
    pass
    
    # Method resolution order:
    #     Model
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __init__(self, *args, **kwargs):
        """Overloaded function.
        
        1. __init__(self: libaster.Model, arg0: ConnectionMesh) -> None
        
        2. __init__(self: libaster.Model, arg0: str, arg1: ConnectionMesh) -> None
        
        3. __init__(self: libaster.Model, arg0: libaster.BaseMesh) -> None
        
        4. __init__(self: libaster.Model, arg0: libaster.BaseMesh, arg1: bool) -> None
        
        5. __init__(self: libaster.Model, arg0: str, arg1: libaster.FiniteElementDescriptor) -> None
        
        6. __init__(self: libaster.Model, arg0: str, arg1: libaster.FiniteElementDescriptor, arg2: bool) -> None
        """
    
    def addModelingOnGroupOfCells(self, arg0, arg1, arg2):
        pass
    
    def addModelingOnGroupOfNodes(self, arg0, arg1, arg2):
        pass
    
    def addModelingOnMesh(self, arg0, arg1):
        pass
    
    def build(self):
        pass
    
    def existsMultiFiberBeam(self):
        pass
    
    def existsThm(self):
        pass
    
    def getConnectionMesh(self):
        """Return the ConnectionMesh
        
        Returns:
            ConnectionMesh: a pointer to the ConnectionMesh
        """
    
    def getFiniteElementDescriptor(self):
        pass
    
    def getGeometricDimension(self):
        """To know the geometric dimension supported by the model
        
        Returns:
            int: geometric dimension
        """
    
    def getGraphPartitioner(self):
        pass
    
    def getMesh(self):
        """Return the mesh
        
        Returns:
            Mesh: a pointer to the mesh
        """
    
    def getPhysics(self):
        """To know the physics supported by the model
        
        Returns:
            str: Mechanics or Thermal or Acoustic
        """
    
    def getSaneModel(self):
        pass
    
    def getSplittingMethod(self):
        pass
    
    def getTable(self, identifier):
        """Extract a Table from the datastructure.
        
        Arguments:
            identifier (str): Table identifier.
        
        Returns:
            Table: Table stored with the given identifier.
        """
    
    def isAcoustic(self):
        """To know if the model is acoustic or not
        
        Returns:
            bool: True - if the model is acoustic
        """
    
    def isMechanical(self):
        """To know if the model is mechanical or not
        
        Returns:
            bool: True - if the model is mechanical
        """
    
    def isThermal(self):
        """To know if the model is thermal or not
        
        Returns:
            bool: True - if the model is thermal
        """
    
    def isXfem(self):
        pass
    
    def setFrom(self, model):
        """Set a model defined on a ConnectionMesh from an other model
        
        Arguments:
            model (Model): Table identifier.
        """
    
    def setSaneModel(self, arg0):
        pass
    
    def setSplittingMethod(self, *args, **kwargs):
        """Overloaded function.
        
        1. setSplittingMethod(self: libaster.Model, arg0: libaster.ModelSplitingMethod, arg1: libaster.GraphPartitioner) -> None
        
        2. setSplittingMethod(self: libaster.Model, arg0: libaster.ModelSplitingMethod) -> None
        """
    
    def xfemPreconditioningEnable(self):
        pass

# class Physics in libaster

class Physics:
    """Members:
    
    Mechanics
    
    Thermal
    
    Acoustic
    """
    
    # Method resolution order:
    #     Physics
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __eq__(self, other):
        pass
    
    def __getstate__(self):
        pass
    
    def __hash__(self):
        pass
    
    def __index__(self):
        pass
    
    def __init__(self, value):
        pass
    
    def __int__(self):
        pass
    
    def __ne__(self, other):
        pass
    
    def __repr__(self):
        pass
    
    def __setstate__(self, state):
        pass
    
    def name(self):
        pass
    
    #----------------------------------------------------------------------
    # Data descriptors defined here:
    
    @property
    def __members__(self):
        pass
    
    @property
    def name(self):
        """name(self: handle) -> str
        """
    
    @property
    def value(self):
        pass
    
    #----------------------------------------------------------------------
    # Data and other attributes defined here:
    
    Acoustic = 2
    
    Mechanics = 0
    
    Thermal = 1

# class Modelings in libaster

class Modelings:
    """Members:
    
    Axisymmetrical
    
    Tridimensional
    
    TridimensionalAbsorbingBoundary
    
    Planar
    
    PlaneStrain
    
    PlaneStress
    
    DKT
    
    DKTG
    
    PlanarBar
    """
    
    # Method resolution order:
    #     Modelings
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __eq__(self, other):
        pass
    
    def __getstate__(self):
        pass
    
    def __hash__(self):
        pass
    
    def __index__(self):
        pass
    
    def __init__(self, value):
        pass
    
    def __int__(self):
        pass
    
    def __ne__(self, other):
        pass
    
    def __repr__(self):
        pass
    
    def __setstate__(self, state):
        pass
    
    def name(self):
        pass
    
    #----------------------------------------------------------------------
    # Data descriptors defined here:
    
    @property
    def __members__(self):
        pass
    
    @property
    def name(self):
        """name(self: handle) -> str
        """
    
    @property
    def value(self):
        pass
    
    #----------------------------------------------------------------------
    # Data and other attributes defined here:
    
    Axisymmetrical = 0
    
    DKT = 6
    
    DKTG = 7
    
    Planar = 3
    
    PlanarBar = 8
    
    PlaneStrain = 4
    
    PlaneStress = 5
    
    Tridimensional = 1
    
    TridimensionalAbsorbingBoundary = 2

# class PrestressingCable in libaster

class PrestressingCable(DataStructure):
    pass
    
    # Method resolution order:
    #     PrestressingCable
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __init__(self, *args, **kwargs):
        """Overloaded function.
        
        1. __init__(self: libaster.PrestressingCable, arg0: libaster.Model, arg1: libaster.MaterialField, arg2: libaster.ElementaryCharacteristics) -> None
        
        2. __init__(self: libaster.PrestressingCable, arg0: str, arg1: libaster.Model, arg2: libaster.MaterialField, arg3: libaster.ElementaryCharacteristics) -> None
        """
    
    def getElementaryCharacteristics(self):
        pass
    
    def getMaterialField(self):
        pass
    
    def getModel(self):
        """Return the Model.
        
        Returns:
            *Model*: Model object.
        """

# class XfemCrack in libaster

class XfemCrack(DataStructure):
    pass
    
    # Method resolution order:
    #     XfemCrack
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __init__(self, *args, **kwargs):
        """Overloaded function.
        
        1. __init__(self: libaster.XfemCrack, arg0: libaster.Mesh) -> None
        
        2. __init__(self: libaster.XfemCrack, arg0: str, arg1: libaster.Mesh) -> None
        """
    
    def build(self):
        pass
    
    def enrichModelWithXfem(self, arg0):
        pass
    
    def getAuxiliaryGrid(self):
        pass
    
    def getCohesiveCrackTipForPropagation(self):
        pass
    
    def getCrackLipsEntity(self):
        pass
    
    def getCrackShape(self):
        pass
    
    def getCrackTipEntity(self):
        pass
    
    def getDiscontinuityType(self):
        pass
    
    def getDiscontinuousField(self):
        pass
    
    def getEnrichedCells(self):
        pass
    
    def getEnrichedLayersNumber(self):
        pass
    
    def getEnrichmentRadiusZone(self):
        pass
    
    def getEnrichmentType(self):
        pass
    
    def getExistingCrackWithGrid(self):
        pass
    
    def getJunctingCracks(self):
        pass
    
    def getMesh(self):
        pass
    
    def getNormalLevelSetField(self):
        pass
    
    def getNormalLevelSetFunction(self):
        pass
    
    def getTable(self, identifier):
        """Extract a Table from the datastructure.
        
        Arguments:
            identifier (str): Table identifier.
        
        Returns:
            Table: Table stored with the given identifier.
        """
    
    def getTangentialLevelSetField(self):
        pass
    
    def getTangentialLevelSetFunction(self):
        pass
    
    def insertJunctingCracks(self, arg0):
        pass
    
    def setAuxiliaryGrid(self, arg0):
        pass
    
    def setCohesiveCrackTipForPropagation(self, arg0):
        pass
    
    def setCrackLipsEntity(self, arg0):
        pass
    
    def setCrackShape(self, arg0):
        pass
    
    def setCrackTipEntity(self, arg0):
        pass
    
    def setDiscontinuityType(self, arg0):
        pass
    
    def setDiscontinuousField(self, arg0):
        pass
    
    def setEnrichedCells(self, arg0):
        pass
    
    def setEnrichedLayersNumber(self, arg0):
        pass
    
    def setEnrichmentRadiusZone(self, arg0):
        pass
    
    def setEnrichmentType(self, arg0):
        pass
    
    def setExistingCrackWithGrid(self, arg0):
        pass
    
    def setMesh(self, arg0):
        pass
    
    def setNormalLevelSetField(self, arg0):
        pass
    
    def setNormalLevelSetFunction(self, arg0):
        pass
    
    def setPointForJunction(self, arg0):
        pass
    
    def setTangentialLevelSetField(self, arg0):
        pass
    
    def setTangentialLevelSetFunction(self, arg0):
        pass
    
    def update(self):
        pass

# class Result in libaster

class Result(DataStructure):
    pass
    
    # Method resolution order:
    #     Result
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __init__(self, *args, **kwargs):
        """Overloaded function.
        
        1. __init__(self: libaster.Result, arg0: str) -> None
        
        2. __init__(self: libaster.Result, arg0: str, arg1: str) -> None
        """
    
    def addFieldOnNodesDescription(self, arg0):
        pass
    
    def allocate(self, nb_rank):
        """Allocate result
        
        Arguments:
            nb_rank (int):  number of rank to allocate
        """
    
    def build(self, feds= [], fnds= []):
        """Build the result from the name of the result. It stores fields which are setted in c++ or
        created in fortran
        
        Arguments:
            feds (list[FiniteElementDescriptor]) : list of additional finite element descriptor used to
                build FieldOnCells
            fnds (list[FieldOnNodesDescriptionPtr]) : list of additional field description used to
                build FieldOnNodes
        
        Returns:
            bool: *True* if ok.
        """
    
    def getAccessParameters(self):
        """Return the access parameters of the result as Python dict.
        
        Returns:
            dict{str : list[int,float,str]}: Dict of values for each access variable.
        """
    
    def getAllElementaryCharacteristics(self):
        """Return the list of all elementary characteristics used in the result
        
        Returns:
            list[ElementaryCharacteristics]: list of ElementaryCharacteristics.
        """
    
    def getConstantFieldOnCellsChar16(self, name, rank):
        """Get a ConstantFieldOnCellsChar16 from result.
        
        Arguments:
            name (str): symbolic name of the field in the result (ex: 'COMPORTEMENT', ...)
            rank (int): rank to set the field
        
        Returns:
            ConstantFieldOnCellsChar16: field to get
        """
    
    def getConstantFieldOnCellsReal(self, name, rank):
        """Get a ConstantFieldOnCellsReal from result.
        
        Arguments:
            name (str): symbolic name of the field in the result (ex: 'DEPL', 'VITE'...)
            rank (int): rank to set the field
        
        Returns:
            ConstantFieldOnCellsReal: field to get
        """
    
    def getConstantFieldsOnCellsChar16Names(self):
        """Return the names of the contant char16 fields on cells as Python list.
        
        Returns:
            list(str): List of names of the contant fields on cells.
        """
    
    def getConstantFieldsOnCellsRealNames(self):
        """Return the names of the contant real fields on cells as Python list.
        
        Returns:
            list(str): List of names of the contant fields on cells.
        """
    
    def getElementaryCharacteristics(self, *args, **kwargs):
        """Overloaded function.
        
        1. getElementaryCharacteristics(self: libaster.Result, rank: int) -> libaster.ElementaryCharacteristics
        
        
        Get elementary characterictics at the specfied rank
        
        Arguments:
            rank (int): rank
        
        Returns:
            ElementaryCharacteristics: a pointer to elementary characterictics.
                
        
        2. getElementaryCharacteristics(self: libaster.Result) -> libaster.ElementaryCharacteristics
        """
    
    def getFieldOnCellsComplex(self, name, rank):
        """Get a FieldOnCellsComplex from result.
        
        Arguments:
            name (str): symbolic name of the field in the result (ex: 'DEPL', 'VITE'...)
            rank (int): rank to set the field
        
        Returns:
            FieldOnCellsComplex: field to get
        """
    
    def getFieldOnCellsLong(self, name, rank):
        """Get a FieldOnCellsLong from result.
        
        Arguments:
            name (str): symbolic name of the field in the result (ex: 'DEPL', 'VITE'...)
            rank (int): rank to set the field
        
        Returns:
            FieldOnCellsLong: field to get
        """
    
    def getFieldOnCellsReal(self, name, rank):
        """Get a FieldOnCellsReal from result.
        
        Arguments:
            name (str): symbolic name of the field in the result (ex: 'DEPL', 'VITE'...)
            rank (int): rank to set the field
        
        Returns:
            FieldOnCellsReal: field to get
        """
    
    def getFieldOnNodesComplex(self, name, rank):
        """Get a FieldOnNodesComplex from result.
        
        Arguments:
            name (str): symbolic name of the field in the result (ex: 'DEPL', 'VITE'...)
            rank (int): rank to set the field
        
        Returns:
            FieldOnNodesComplex: field to get
        """
    
    def getFieldOnNodesDescriptions(self):
        """Get list of field's description to build internal FieldOnNodes
        
        Returns:
            list[FieldOnNodesDescription]: list of field's description
        """
    
    def getFieldOnNodesReal(self, name, rank):
        """Get a FieldOnNodesReal from result.
        
        Arguments:
            name (str): symbolic name of the field in the result (ex: 'DEPL', 'VITE'...)
            rank (int): rank to set the field
        
        Returns:
            FieldOnNodesReal: field to get
        """
    
    def getFieldsNames(self, *args, **kwargs):
        """Overloaded function.
        
        1. getFieldsNames(self: libaster.Result) -> List[str]
        
        
        Return the list of names of stored fields
        
        Returns:
            list[str]: List of names of stored fields.
                
        
        2. getFieldsNames(self: libaster.Result) -> List[str]
        
        
        Return the list of names of stored fields
        
        Returns:
            list[str]: List of names of stored fields.
        """
    
    def getFieldsOnCellsComplexNames(self):
        """Return the names of the complex fields on cells as Python list.
        
        Returns:
            list(str): List of names of the complex fields on cells.
        """
    
    def getFieldsOnCellsLongNames(self):
        """Return the names of the integer fields on cells as Python list.
        
        Returns:
            list(str): List of names of the integer fields on cells.
        """
    
    def getFieldsOnCellsRealNames(self):
        """Return the names of the real fields on cells as Python list.
        
        Returns:
            list(str): List of names of the real fields on cells.
        """
    
    def getFieldsOnNodesComplexNames(self):
        """Return the names of the complex fields on nodes as Python list.
        
        Returns:
            list(str): List of names of the complex fields on nodes.
        """
    
    def getFieldsOnNodesRealNames(self):
        """Return the names of the real fields on nodes as Python list.
        
        Returns:
            list(str): List of names of the real fields on nodes.
        """
    
    def getFiniteElementDescriptors(self):
        """Get list of finite element descriptor to build internal FieldOnCells
        
        Returns:
            list[FiniteElementDescriptor]: list of finite element descriptor
        """
    
    def getListOfLoads(self, rank):
        """Get list of loads on the specified rank
        
        Arguments:
            rank (int): rank to get
        
        Returns:
            ListOfLoads: a pointer to list of loads.
        """
    
    def getMaterialField(self, *args, **kwargs):
        """Overloaded function.
        
        1. getMaterialField(self: libaster.Result, rank: int) -> libaster.MaterialField
        
        
        Return the material field for the given rank.
        
        Arguments:
            rank (int): rank
        
        Returns:
            MaterialField: Material field.
                      
        
        2. getMaterialField(self: libaster.Result) -> libaster.MaterialField
        """
    
    def getMaterialFields(self):
        """Return the list of all material fields used in the result
        
        Returns:
            list[MaterialField]: list of material field.
        """
    
    def getMesh(self):
        """Return a pointer to mesh
        
        Returns:
            mesh (Mesh): a pointer to the mesh.
        """
    
    def getModel(self, *args, **kwargs):
        """Overloaded function.
        
        1. getModel(self: libaster.Result, rank: int) -> libaster.Model
        
        
        Return the model for the given rank.
        
        Arguments:
            rank (int): rank
        
        Returns:
            Model: Model object.
                      
        
        2. getModel(self: libaster.Result) -> libaster.Model
        """
    
    def getModels(self):
        """Return the list of all models used in the result
        
        Returns:
            list[Model]: list of models.
        """
    
    def getNumberOfRanks(self):
        """Get the number of rank stored in the result
        
        Returns:
            int: number of rank stored.
        """
    
    def getRanks(self):
        """Return the list of ranks used to store fields
        
        Returns:
            list[int]: List of ranks used to store fields.
        """
    
    def getTable(self, identifier):
        """Extract a Table from the datastructure.
        
        Arguments:
            identifier (str): Table identifier.
        
        Returns:
            Table: Table stored with the given identifier.
        """
    
    def getTimeValue(self, rank):
        """Get time at the specified rank
        
        Arguments:
            rank (int):  rank where to save time value
        
        Returns
            float: time value
        """
    
    def hasElementaryCharacteristics(self, *args, **kwargs):
        """Overloaded function.
        
        1. hasElementaryCharacteristics(self: libaster.Result, rank: int) -> bool
        
        
        Test if a elementary characterictics is used at the specfied rank
        
        Arguments:
            rank (int): rank
        
        Returns:
            bool: *True* if at least one elementary characterictics used else *False*.
                
        
        2. hasElementaryCharacteristics(self: libaster.Result) -> bool
        """
    
    def hasListOfLoads(self, *args, **kwargs):
        """Overloaded function.
        
        1. hasListOfLoads(self: libaster.Result, rank: int) -> bool
        
        
        Test if a list of loads is used at the specfied rank
        
        Arguments:
            rank (int): rank
        
        Returns:
            bool: *True* if at least one list of loads is used else *False*.
                
        
        2. hasListOfLoads(self: libaster.Result) -> bool
        """
    
    def hasMaterialField(self, rank):
        """Test if a material field is used at the specfied rank
        
        Arguments:
            rank (int): rank
        
        Returns:
            bool: *True* if at a material field used else *False*.
        """
    
    def hasModel(self, rank):
        """Test if a model is used at the specfied rank
        
        Arguments:
            rank (int): rank
        
        Returns:
            bool: *True* if at a model used else *False*.
        """
    
    def printInfo(self):
        pass
    
    def printListOfFields(self):
        """Print the names of all fields (real, complex, ...) stored in the result.
        """
    
    def printMedFile(self, filename, medname= '', local= True):
        """Print the result in a MED file.
        
        Args:
            filename (str): Path to the output file.
            medname (str): Name of the result in the MED file. (default: "")
            local (bool): Print only the local domain if *True*. (default: True)
        """
    
    def resize(self, nbRanks):
        """Resize the object.
        
        Arguments:
            nbRanks (int): new expected size. Should be greater than the current size,
                otherwise the size is unchanged.
        """
    
    def setElementaryCharacteristics(self, *args, **kwargs):
        """Overloaded function.
        
        1. setElementaryCharacteristics(self: libaster.Result, cara_elem: libaster.ElementaryCharacteristics) -> None
        
        
        Set elementary characterictics on all ranks
        
        Arguments:
            cara_elem (ElementaryCharacteristics): elementary characterictics to set.
                
        
        2. setElementaryCharacteristics(self: libaster.Result, cara_elem: libaster.ElementaryCharacteristics, rank: int) -> None
        
        
        Set elementary characterictics on the specified rank
        
        Arguments:
            cara_elem (ElementaryCharacteristics): elementary characterictics to set.
            rank (int): rank to set
        """
    
    def setField(self, *args, **kwargs):
        """Overloaded function.
        
        1. setField(self: libaster.Result, field: libaster.FieldOnNodesReal, name: str, rank: int) -> None
        
        
        Set a real FieldOnNodes to result.
        
        Arguments:
            field (FieldOnNodesReal): field to set
            name (str): symbolic name of the field in the result (ex: 'DEPL', 'VITE'...)
            rank (int): rank to set the field
                
        
        2. setField(self: libaster.Result, field: libaster.FieldOnNodesComplex, name: str, rank: int) -> None
        
        
        Set a complex FieldOnNodes to result.
        
        Arguments:
            field (FieldOnNodesComplex): field to set
            name (str): symbolic name of the field in the result (ex: 'DEPL', 'VITE'...)
            rank (int): rank to set the field
                
        
        3. setField(self: libaster.Result, field: libaster.FieldOnCellsReal, name: str, rank: int) -> None
        
        
        Set a real FieldOnCells to result
        
        Arguments:
            field (FieldOnCellsReal): field to set
            name (str): symbolic name of the field in the result (ex: 'VARI_ELGA', 'SIEF_ELGA'...)
            rank (int): rank to set the field
                
        
        4. setField(self: libaster.Result, field: libaster.FieldOnCellsComplex, name: str, rank: int) -> None
        
        
        Set a complex FieldOnCells to result
        
        Arguments:
            field (FieldOnCellsComplex): field to set
            name (str): symbolic name of the field in the result (ex: 'VARI_ELGA', 'SIEF_ELGA'...)
            rank (int): rank to set the field
                
        
        5. setField(self: libaster.Result, field: libaster.FieldOnCellsLong, name: str, rank: int) -> None
        
        
        Set a long FieldOnCells to result
        
        Arguments:
            field (FieldOnCellsLong): field to set
            name (str): symbolic name of the field in the result (ex: 'VARI_ELGA', 'SIEF_ELGA'...)
            rank (int): rank to set the field
                
        
        6. setField(self: libaster.Result, field: libaster.ConstantFieldOnCellsChar16, name: str, rank: int) -> None
        
        
        Set a ConstantFieldOnCellsChar16 to result
        
        Arguments:
            field (ConstantFieldOnCellsChar16): field to set
            name (str): symbolic name of the field in the result (ex: 'COMPOR', ...)
            rank (int): rank to set the field
                
        
        7. setField(self: libaster.Result, field: libaster.ConstantFieldOnCellsReal, name: str, rank: int) -> None
        
        
        Set a ConstantFieldOnCellsReal to result
        
        Arguments:
            field (ConstantFieldOnCellsReal): field to set
            name (str): symbolic name of the field in the result (ex: 'COMPOR', ...)
            rank (int): rank to set the field
        """
    
    def setListOfLoads(self, load, rank):
        """Set list of loads on the specified rank
        
        Arguments:
            load (ListOfLoads): list of loads to set.
            rank (int): rank to set
        """
    
    def setMaterialField(self, *args, **kwargs):
        """Overloaded function.
        
        1. setMaterialField(self: libaster.Result, mater: libaster.MaterialField) -> None
        
        
        Set material field on all ranks
        
        Arguments:
            mater (MaterialField): material field to set.
                
        
        2. setMaterialField(self: libaster.Result, mater: libaster.MaterialField, rank: int) -> None
        
        
        Set material field on the specified rank
        
        Arguments:
            mater (MaterialField): material field to set.
            rank (int): rank to set
        """
    
    def setMesh(self, mesh):
        """Set the mesh used by the result.
        
        Arguments:
            mesh (BaseMesh): mesh to set
        """
    
    def setModel(self, *args, **kwargs):
        """Overloaded function.
        
        1. setModel(self: libaster.Result, model: libaster.Model) -> None
        
        
        Set model on all ranks
        
        Arguments:
            model (Model): model to set.
                
        
        2. setModel(self: libaster.Result, model: libaster.Model, rank: int) -> None
        
        
        Set model on the specified rank
        
        Arguments:
            model (Model): model to set
            rank (int): rank to set
        """
    
    def setParameterValue(self, para_name, value, rank):
        """Add theta at the specified rank
        
        Arguments:
            name (float): parameter name to store
            value (float): parameter value to store
            rank (int):  rank where to save time value
        """
    
    def setTimeValue(self, time, rank):
        """Add time at the specified rank
        
        Arguments:
            time (float): time value to save
            rank (int):  rank where to save time value
        """

# class TransientResult in libaster

class TransientResult(Result):
    pass
    
    # Method resolution order:
    #     TransientResult
    #     Result
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __init__(self, *args, **kwargs):
        """Overloaded function.
        
        1. __init__(self: libaster.TransientResult) -> None
        
        2. __init__(self: libaster.TransientResult, arg0: str, arg1: str) -> None
        """

# class LoadResult in libaster

class LoadResult(TransientResult):
    pass
    
    # Method resolution order:
    #     LoadResult
    #     TransientResult
    #     Result
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __init__(self, *args, **kwargs):
        """Overloaded function.
        
        1. __init__(self: libaster.LoadResult) -> None
        
        2. __init__(self: libaster.LoadResult, arg0: str) -> None
        """

# class ThermalResult in libaster

class ThermalResult(TransientResult):
    pass
    
    # Method resolution order:
    #     ThermalResult
    #     TransientResult
    #     Result
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __init__(self, *args, **kwargs):
        """Overloaded function.
        
        1. __init__(self: libaster.ThermalResult) -> None
        
        2. __init__(self: libaster.ThermalResult, arg0: str) -> None
        """

# class CombinedFourierResult in libaster

class CombinedFourierResult(Result):
    pass
    
    # Method resolution order:
    #     CombinedFourierResult
    #     Result
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __init__(self, *args, **kwargs):
        """Overloaded function.
        
        1. __init__(self: libaster.CombinedFourierResult) -> None
        
        2. __init__(self: libaster.CombinedFourierResult, arg0: str) -> None
        """

# class ElasticFourierResult in libaster

class ElasticFourierResult(Result):
    pass
    
    # Method resolution order:
    #     ElasticFourierResult
    #     Result
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __init__(self, *args, **kwargs):
        """Overloaded function.
        
        1. __init__(self: libaster.ElasticFourierResult) -> None
        
        2. __init__(self: libaster.ElasticFourierResult, arg0: str) -> None
        """

# class ThermalFourierResult in libaster

class ThermalFourierResult(Result):
    pass
    
    # Method resolution order:
    #     ThermalFourierResult
    #     Result
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __init__(self, *args, **kwargs):
        """Overloaded function.
        
        1. __init__(self: libaster.ThermalFourierResult) -> None
        
        2. __init__(self: libaster.ThermalFourierResult, arg0: str) -> None
        """

# class MultipleElasticResult in libaster

class MultipleElasticResult(Result):
    pass
    
    # Method resolution order:
    #     MultipleElasticResult
    #     Result
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __init__(self, *args, **kwargs):
        """Overloaded function.
        
        1. __init__(self: libaster.MultipleElasticResult) -> None
        
        2. __init__(self: libaster.MultipleElasticResult, arg0: str) -> None
        """

# class NonLinearResult in libaster

class NonLinearResult(TransientResult):
    pass
    
    # Method resolution order:
    #     NonLinearResult
    #     TransientResult
    #     Result
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __init__(self, *args, **kwargs):
        """Overloaded function.
        
        1. __init__(self: libaster.NonLinearResult) -> None
        
        2. __init__(self: libaster.NonLinearResult, arg0: str) -> None
        """
    
    def setContact(self, *args, **kwargs):
        """Overloaded function.
        
        1. setContact(self: libaster.NonLinearResult, arg0: libaster.Contact) -> None
        
        2. setContact(self: libaster.NonLinearResult, arg0: libaster.Contact, arg1: int) -> None
        """

# class PhysicalProblem in libaster

class PhysicalProblem:
    pass
    
    # Method resolution order:
    #     PhysicalProblem
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __init__(self, *args, **kwargs):
        """Overloaded function.
        
        1. __init__(self: libaster.PhysicalProblem, arg0: libaster.Model, arg1: libaster.MaterialField) -> None
        
        2. __init__(self: libaster.PhysicalProblem, arg0: libaster.Model, arg1: libaster.MaterialField, arg2: libaster.ElementaryCharacteristics) -> None
        """
    
    def addDirichletBC(self, *args, **kwargs):
        """Overloaded function.
        
        1. addDirichletBC(self: libaster.PhysicalProblem, arg0: libaster.DirichletBC) -> None
        
        2. addDirichletBC(self: libaster.PhysicalProblem, arg0: libaster.DirichletBC, arg1: libaster.Function) -> None
        
        3. addDirichletBC(self: libaster.PhysicalProblem, arg0: libaster.DirichletBC, arg1: libaster.Formula) -> None
        
        4. addDirichletBC(self: libaster.PhysicalProblem, arg0: libaster.DirichletBC, arg1: libaster.Function2D) -> None
        """
    
    def addLoad(self, *args, **kwargs):
        """Overloaded function.
        
        1. addLoad(self: libaster.PhysicalProblem, arg0: libaster.MechanicalLoadReal) -> None
        
        2. addLoad(self: libaster.PhysicalProblem, arg0: libaster.MechanicalLoadReal, arg1: libaster.Function) -> None
        
        3. addLoad(self: libaster.PhysicalProblem, arg0: libaster.MechanicalLoadReal, arg1: libaster.Formula) -> None
        
        4. addLoad(self: libaster.PhysicalProblem, arg0: libaster.MechanicalLoadReal, arg1: libaster.Function2D) -> None
        
        5. addLoad(self: libaster.PhysicalProblem, arg0: libaster.MechanicalLoadFunction) -> None
        
        6. addLoad(self: libaster.PhysicalProblem, arg0: libaster.MechanicalLoadFunction, arg1: libaster.Function) -> None
        
        7. addLoad(self: libaster.PhysicalProblem, arg0: libaster.MechanicalLoadFunction, arg1: libaster.Formula) -> None
        
        8. addLoad(self: libaster.PhysicalProblem, arg0: libaster.MechanicalLoadFunction, arg1: libaster.Function2D) -> None
        
        9. addLoad(self: libaster.PhysicalProblem, arg0: libaster.MechanicalLoadComplex) -> None
        
        10. addLoad(self: libaster.PhysicalProblem, arg0: libaster.MechanicalLoadComplex, arg1: libaster.Function) -> None
        
        11. addLoad(self: libaster.PhysicalProblem, arg0: libaster.MechanicalLoadComplex, arg1: libaster.Formula) -> None
        
        12. addLoad(self: libaster.PhysicalProblem, arg0: libaster.MechanicalLoadComplex, arg1: libaster.Function2D) -> None
        
        13. addLoad(self: libaster.PhysicalProblem, arg0: ParallelMechanicalLoad<ConstantFieldOnCells<double> >) -> None
        
        14. addLoad(self: libaster.PhysicalProblem, arg0: ParallelMechanicalLoad<ConstantFieldOnCells<double> >, arg1: libaster.Function) -> None
        
        15. addLoad(self: libaster.PhysicalProblem, arg0: ParallelMechanicalLoad<ConstantFieldOnCells<double> >, arg1: libaster.Formula) -> None
        
        16. addLoad(self: libaster.PhysicalProblem, arg0: ParallelMechanicalLoad<ConstantFieldOnCells<double> >, arg1: libaster.Function2D) -> None
        
        17. addLoad(self: libaster.PhysicalProblem, arg0: ParallelMechanicalLoad<ConstantFieldOnCells<JeveuxString<24> > >) -> None
        
        18. addLoad(self: libaster.PhysicalProblem, arg0: ParallelMechanicalLoad<ConstantFieldOnCells<JeveuxString<24> > >, arg1: libaster.Function) -> None
        
        19. addLoad(self: libaster.PhysicalProblem, arg0: ParallelMechanicalLoad<ConstantFieldOnCells<JeveuxString<24> > >, arg1: libaster.Formula) -> None
        
        20. addLoad(self: libaster.PhysicalProblem, arg0: ParallelMechanicalLoad<ConstantFieldOnCells<JeveuxString<24> > >, arg1: libaster.Function2D) -> None
        
        21. addLoad(self: libaster.PhysicalProblem, arg0: ParallelThermalLoad<ConstantFieldOnCells<double> >) -> None
        
        22. addLoad(self: libaster.PhysicalProblem, arg0: ParallelThermalLoad<ConstantFieldOnCells<double> >, arg1: libaster.Function) -> None
        
        23. addLoad(self: libaster.PhysicalProblem, arg0: ParallelThermalLoad<ConstantFieldOnCells<double> >, arg1: libaster.Formula) -> None
        
        24. addLoad(self: libaster.PhysicalProblem, arg0: ParallelThermalLoad<ConstantFieldOnCells<double> >, arg1: libaster.Function2D) -> None
        
        25. addLoad(self: libaster.PhysicalProblem, arg0: ParallelThermalLoad<ConstantFieldOnCells<JeveuxString<24> > >) -> None
        
        26. addLoad(self: libaster.PhysicalProblem, arg0: ParallelThermalLoad<ConstantFieldOnCells<JeveuxString<24> > >, arg1: libaster.Function) -> None
        
        27. addLoad(self: libaster.PhysicalProblem, arg0: ParallelThermalLoad<ConstantFieldOnCells<JeveuxString<24> > >, arg1: libaster.Formula) -> None
        
        28. addLoad(self: libaster.PhysicalProblem, arg0: ParallelThermalLoad<ConstantFieldOnCells<JeveuxString<24> > >, arg1: libaster.Function2D) -> None
        
        29. addLoad(self: libaster.PhysicalProblem, arg0: libaster.ThermalLoadReal) -> None
        
        30. addLoad(self: libaster.PhysicalProblem, arg0: libaster.ThermalLoadReal, arg1: libaster.Function) -> None
        
        31. addLoad(self: libaster.PhysicalProblem, arg0: libaster.ThermalLoadReal, arg1: libaster.Formula) -> None
        
        32. addLoad(self: libaster.PhysicalProblem, arg0: libaster.ThermalLoadReal, arg1: libaster.Function2D) -> None
        
        33. addLoad(self: libaster.PhysicalProblem, arg0: libaster.ThermalLoadFunction) -> None
        
        34. addLoad(self: libaster.PhysicalProblem, arg0: libaster.ThermalLoadFunction, arg1: libaster.Function) -> None
        
        35. addLoad(self: libaster.PhysicalProblem, arg0: libaster.ThermalLoadFunction, arg1: libaster.Formula) -> None
        
        36. addLoad(self: libaster.PhysicalProblem, arg0: libaster.ThermalLoadFunction, arg1: libaster.Function2D) -> None
        
        37. addLoad(self: libaster.PhysicalProblem, arg0: libaster.AcousticLoadComplex) -> None
        
        38. addLoad(self: libaster.PhysicalProblem, arg0: libaster.AcousticLoadComplex, arg1: libaster.Function) -> None
        
        39. addLoad(self: libaster.PhysicalProblem, arg0: libaster.AcousticLoadComplex, arg1: libaster.Formula) -> None
        
        40. addLoad(self: libaster.PhysicalProblem, arg0: libaster.AcousticLoadComplex, arg1: libaster.Function2D) -> None
        """
    
    def computeBehaviourProperty(self, COMPORTEMENT, SIGM_INIT= 'NON', INFO= 1):
        """Create constant fields on cells for behaviour (COMPOR, CARCRI and MULCOM)
        
        Arguments:
            COMPORTEMENT (list[dict]): keywords as provided to STAT_NON_LINE/COMPORTEMENT
            SIGM_INIT (str): "OUI" if there is an initial stress field
            INFO (int): level of verbosity, 1 to have description of behaviour or 0 to be quiet
        """
    
    def computeDOFNumbering(self):
        """Build DOF numbering from the model and loads
        
        Returns:
            Bool: True if success
        """
    
    def computeListOfLoads(self):
        """Build the list of loads from the added loads
        
        Returns:
            Bool: True if success
        """
    
    def computeReferenceExternalStateVariables(self):
        """Compute field for external state variables reference value
        
        Returns:
            FieldOnCells: field for external state variables reference values
        """
    
    def getBehaviourProperty(self):
        """Return the behaviour properties
        
        Returns:
            BehaviourPropertyPtr: a pointer to the behaviour properties
        """
    
    def getCodedMaterial(self):
        """Return the coded material
        
        Returns:
            CodedMaterialPtr: a pointer to the coded material
        """
    
    def getDOFNumbering(self):
        """Return the DOF numbering
        
        Returns:
            BaseDOFNumberingPtr: a pointer to the DOF numbering
        """
    
    def getElementaryCharacteristics(self):
        """Return the elementary charateristics
        
        Returns:
            ElementaryCharacteristicsPtr: a pointer to the elementary charateristics
        """
    
    def getExternalStateVariables(self, time):
        """Get the field for external state variables
        
        Arguments:
            time [float] : time value to evaluate values
        
        Returns:
            FieldOnCellsRealPtr : external values
        """
    
    def getListOfLoads(self):
        """Return list of loads.
        
        Returns:
            ListOfLoadsPtr: a pointer to list of loads
        """
    
    def getMaterialField(self):
        """Return the material field
        
        Returns:
            MaterialFieldPtr: a pointer to the material field
        """
    
    def getMesh(self):
        """Return the mesh
        
        Returns:
            MeshPtr: a pointer to the mesh
        """
    
    def getModel(self):
        """Return the model
        
        Returns:
            ModelPtr: a pointer to the model
        """
    
    def getReferenceExternalStateVariables(self):
        """Get the field of reference values for external state variables
        
        Returns:
            FieldOnCellsRealPtr : field of reference values
        """
    
    def setDOFNumbering(self, dofNume):
        """Set the DOF numbering
        
        Arguments:
            BaseDOFNumberingPtr: a pointer to the DOF numbering
        """

# class Glossary in libaster

class Glossary:
    pass
    
    # Method resolution order:
    #     Glossary
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __init__(self,  *args, **kwargs):
        """Initialize self.  See help(type(self)) for accurate signature.
        """
    
    def getComponent(self, arg0):
        pass
    
    def getModeling(self, arg0):
        pass
    
    def getPhysics(self, arg0):
        pass

# built-in function getGlossary in libaster

def getGlossary():
    pass

# class CyclicSymmetryMode in libaster

class CyclicSymmetryMode(DataStructure):
    pass
    
    # Method resolution order:
    #     CyclicSymmetryMode
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __init__(self, *args, **kwargs):
        """Overloaded function.
        
        1. __init__(self: libaster.CyclicSymmetryMode) -> None
        
        2. __init__(self: libaster.CyclicSymmetryMode, arg0: str) -> None
        """

# class FullResult in libaster

class FullResult(Result):
    pass
    
    # Method resolution order:
    #     FullResult
    #     Result
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __init__(self, *args, **kwargs):
        """Overloaded function.
        
        1. __init__(self: libaster.FullResult, arg0: str, arg1: str) -> None
        
        2. __init__(self: libaster.FullResult, arg0: str) -> None
        """
    
    def getDOFNumbering(self):
        pass
    
    def setDOFNumbering(self, arg0):
        pass

# class ModeResult in libaster

class ModeResult(FullResult):
    pass
    
    # Method resolution order:
    #     ModeResult
    #     FullResult
    #     Result
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __init__(self, *args, **kwargs):
        """Overloaded function.
        
        1. __init__(self: libaster.ModeResult) -> None
        
        2. __init__(self: libaster.ModeResult, arg0: str) -> None
        """
    
    def getDOFNumbering(self):
        pass
    
    def getMassMatrix(self):
        pass
    
    def getStiffnessMatrix(self):
        pass
    
    def setMassMatrix(self, *args, **kwargs):
        """Overloaded function.
        
        1. setMassMatrix(self: libaster.ModeResult, arg0: libaster.AssemblyMatrixDisplacementReal) -> None
        
        2. setMassMatrix(self: libaster.ModeResult, arg0: libaster.AssemblyMatrixTemperatureReal) -> None
        
        3. setMassMatrix(self: libaster.ModeResult, arg0: libaster.AssemblyMatrixDisplacementComplex) -> None
        
        4. setMassMatrix(self: libaster.ModeResult, arg0: libaster.GeneralizedAssemblyMatrixComplex) -> None
        """
    
    def setStiffnessMatrix(self, *args, **kwargs):
        """Overloaded function.
        
        1. setStiffnessMatrix(self: libaster.ModeResult, arg0: libaster.AssemblyMatrixDisplacementReal) -> None
        
        2. setStiffnessMatrix(self: libaster.ModeResult, arg0: libaster.AssemblyMatrixTemperatureReal) -> None
        
        3. setStiffnessMatrix(self: libaster.ModeResult, arg0: libaster.AssemblyMatrixDisplacementComplex) -> None
        
        4. setStiffnessMatrix(self: libaster.ModeResult, arg0: libaster.AssemblyMatrixPressureReal) -> None
        
        5. setStiffnessMatrix(self: libaster.ModeResult, arg0: libaster.AssemblyMatrixPressureReal) -> None
        
        6. setStiffnessMatrix(self: libaster.ModeResult, arg0: libaster.GeneralizedAssemblyMatrixReal) -> None
        """
    
    def setStructureInterface(self, arg0):
        pass

# class ModeResultComplex in libaster

class ModeResultComplex(ModeResult):
    pass
    
    # Method resolution order:
    #     ModeResultComplex
    #     ModeResult
    #     FullResult
    #     Result
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __init__(self, *args, **kwargs):
        """Overloaded function.
        
        1. __init__(self: libaster.ModeResultComplex) -> None
        
        2. __init__(self: libaster.ModeResultComplex, arg0: str) -> None
        """
    
    def setDampingMatrix(self, arg0):
        pass
    
    def setStiffnessMatrix(self, *args, **kwargs):
        """Overloaded function.
        
        1. setStiffnessMatrix(self: libaster.ModeResultComplex, arg0: libaster.AssemblyMatrixDisplacementReal) -> bool
        
        2. setStiffnessMatrix(self: libaster.ModeResultComplex, arg0: libaster.AssemblyMatrixDisplacementComplex) -> bool
        
        3. setStiffnessMatrix(self: libaster.ModeResultComplex, arg0: libaster.AssemblyMatrixTemperatureReal) -> bool
        
        4. setStiffnessMatrix(self: libaster.ModeResultComplex, arg0: libaster.AssemblyMatrixPressureReal) -> bool
        
        5. setStiffnessMatrix(self: libaster.ModeResultComplex, arg0: libaster.GeneralizedAssemblyMatrixReal) -> bool
        
        6. setStiffnessMatrix(self: libaster.ModeResultComplex, arg0: libaster.GeneralizedAssemblyMatrixComplex) -> bool
        """
    
    def setStructureInterface(self, arg0):
        pass

# class AcousticModeResult in libaster

class AcousticModeResult(FullResult):
    pass
    
    # Method resolution order:
    #     AcousticModeResult
    #     FullResult
    #     Result
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __init__(self, *args, **kwargs):
        """Overloaded function.
        
        1. __init__(self: libaster.AcousticModeResult) -> None
        
        2. __init__(self: libaster.AcousticModeResult, arg0: str) -> None
        """
    
    def setStiffnessMatrix(self, arg0):
        pass

# class BucklingModeResult in libaster

class BucklingModeResult(FullResult):
    pass
    
    # Method resolution order:
    #     BucklingModeResult
    #     FullResult
    #     Result
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __init__(self, *args, **kwargs):
        """Overloaded function.
        
        1. __init__(self: libaster.BucklingModeResult) -> None
        
        2. __init__(self: libaster.BucklingModeResult, arg0: str) -> None
        """
    
    def getStiffnessMatrix(self):
        pass
    
    def setStiffnessMatrix(self, *args, **kwargs):
        """Overloaded function.
        
        1. setStiffnessMatrix(self: libaster.BucklingModeResult, arg0: libaster.AssemblyMatrixDisplacementReal) -> bool
        
        2. setStiffnessMatrix(self: libaster.BucklingModeResult, arg0: libaster.AssemblyMatrixDisplacementComplex) -> bool
        
        3. setStiffnessMatrix(self: libaster.BucklingModeResult, arg0: libaster.AssemblyMatrixTemperatureReal) -> bool
        
        4. setStiffnessMatrix(self: libaster.BucklingModeResult, arg0: libaster.AssemblyMatrixPressureReal) -> bool
        
        5. setStiffnessMatrix(self: libaster.BucklingModeResult, arg0: libaster.GeneralizedAssemblyMatrixReal) -> bool
        
        6. setStiffnessMatrix(self: libaster.BucklingModeResult, arg0: libaster.GeneralizedAssemblyMatrixComplex) -> bool
        """

# class GeneralizedResultReal in libaster

class GeneralizedResultReal(DataStructure):
    pass
    
    # Method resolution order:
    #     GeneralizedResultReal
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __init__(self,  *args, **kwargs):
        """Initialize self.  See help(type(self)) for accurate signature.
        """

# class GeneralizedResultComplex in libaster

class GeneralizedResultComplex(DataStructure):
    pass
    
    # Method resolution order:
    #     GeneralizedResultComplex
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __init__(self,  *args, **kwargs):
        """Initialize self.  See help(type(self)) for accurate signature.
        """

# class TransientGeneralizedResult in libaster

class TransientGeneralizedResult(GeneralizedResultReal):
    pass
    
    # Method resolution order:
    #     TransientGeneralizedResult
    #     GeneralizedResultReal
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __init__(self, *args, **kwargs):
        """Overloaded function.
        
        1. __init__(self: libaster.TransientGeneralizedResult) -> None
        
        2. __init__(self: libaster.TransientGeneralizedResult, arg0: str) -> None
        """
    
    def getDOFNumbering(self):
        pass
    
    def getGeneralizedDOFNumbering(self):
        pass
    
    def setDOFNumbering(self, arg0):
        pass
    
    def setGeneralizedDOFNumbering(self, arg0):
        pass

# class HarmoGeneralizedResult in libaster

class HarmoGeneralizedResult(GeneralizedResultComplex):
    pass
    
    # Method resolution order:
    #     HarmoGeneralizedResult
    #     GeneralizedResultComplex
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __init__(self, *args, **kwargs):
        """Overloaded function.
        
        1. __init__(self: libaster.HarmoGeneralizedResult) -> None
        
        2. __init__(self: libaster.HarmoGeneralizedResult, arg0: str) -> None
        """
    
    def getDOFNumbering(self):
        pass
    
    def getGeneralizedDOFNumbering(self):
        pass
    
    def setDOFNumbering(self, arg0):
        pass
    
    def setGeneralizedDOFNumbering(self, arg0):
        pass

# class ElasticResult in libaster

class ElasticResult(Result):
    pass
    
    # Method resolution order:
    #     ElasticResult
    #     Result
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __init__(self, *args, **kwargs):
        """Overloaded function.
        
        1. __init__(self: libaster.ElasticResult) -> None
        
        2. __init__(self: libaster.ElasticResult, arg0: str) -> None
        """

# class MeshCoordinatesField in libaster

class MeshCoordinatesField(DataStructure):
    pass
    
    # Method resolution order:
    #     MeshCoordinatesField
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __add__(self, *args, **kwargs):
        """Overloaded function.
        
        1. __add__(self: libaster.MeshCoordinatesField, arg0: libaster.MeshCoordinatesField) -> libaster.MeshCoordinatesField
        
        2. __add__(self: libaster.MeshCoordinatesField, arg0: libaster.FieldOnNodesReal) -> libaster.MeshCoordinatesField
        
        3. __add__(self: libaster.FieldOnNodesReal, arg0: libaster.MeshCoordinatesField) -> libaster.MeshCoordinatesField
        """
    
    def __getitem__(self, idx):
        """Return the coordinate at index *idx* in the vector.
        
        The value is the same as *getValues()[idx]* without creating the entire vector.
        
        Returns:
            float: Values of the *idx*-th coordinate.
        """
    
    def __iadd__(self, arg0):
        pass
    
    def __imul__(self, arg0):
        pass
    
    def __init__(self, arg0):
        pass
    
    def __isub__(self, arg0):
        pass
    
    def __mul__(self, arg0):
        pass
    
    def __neg__(self):
        pass
    
    def __rmul__(self, arg0):
        pass
    
    def __sub__(self, arg0):
        pass
    
    def duplicate(self):
        """Return a copy of MeshCoordinatesField object
        
        Returns:
            MeshCoordinatesField : MeshCoordinatesField object
        """
    
    def getValues(self):
        """Return a list of values of the coordinates as (x1, y1, z1, x2, y2, z2...)
        
        Returns:
            list[float]: List of coordinates (size = 3 * number of nodes).
        """
    
    def size(self):
        """Return the size of the field
        
        Returns:
            int : number of values of MeshCoordinatesField object
        """
    
    def updateValuePointers(self):
        """Update values of internal pointer.
        """

# class FullTransientResult in libaster

class FullTransientResult(FullResult):
    pass
    
    # Method resolution order:
    #     FullTransientResult
    #     FullResult
    #     Result
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __init__(self, *args, **kwargs):
        """Overloaded function.
        
        1. __init__(self: libaster.FullTransientResult, arg0: str) -> None
        
        2. __init__(self: libaster.FullTransientResult) -> None
        """

# class FullHarmonicResult in libaster

class FullHarmonicResult(FullResult):
    pass
    
    # Method resolution order:
    #     FullHarmonicResult
    #     FullResult
    #     Result
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __init__(self, *args, **kwargs):
        """Overloaded function.
        
        1. __init__(self: libaster.FullHarmonicResult, arg0: str) -> None
        
        2. __init__(self: libaster.FullHarmonicResult) -> None
        """

# class FullHarmonicAcousticResult in libaster

class FullHarmonicAcousticResult(FullResult):
    pass
    
    # Method resolution order:
    #     FullHarmonicAcousticResult
    #     FullResult
    #     Result
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __init__(self, *args, **kwargs):
        """Overloaded function.
        
        1. __init__(self: libaster.FullHarmonicAcousticResult, arg0: str) -> None
        
        2. __init__(self: libaster.FullHarmonicAcousticResult) -> None
        """

# class FluidStructureModalBasis in libaster

class FluidStructureModalBasis(DataStructure):
    pass
    
    # Method resolution order:
    #     FluidStructureModalBasis
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __init__(self, *args, **kwargs):
        """Overloaded function.
        
        1. __init__(self: libaster.FluidStructureModalBasis) -> None
        
        2. __init__(self: libaster.FluidStructureModalBasis, arg0: str) -> None
        """
    
    def getTable(self, identifier):
        """Extract a Table from the datastructure.
        
        Arguments:
            identifier (str): Table identifier.
        
        Returns:
            Table: Table stored with the given identifier.
        """

# class GeneralizedModeResult in libaster

class GeneralizedModeResult(FullResult):
    pass
    
    # Method resolution order:
    #     GeneralizedModeResult
    #     FullResult
    #     Result
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __init__(self, *args, **kwargs):
        """Overloaded function.
        
        1. __init__(self: libaster.GeneralizedModeResult, arg0: str) -> None
        
        2. __init__(self: libaster.GeneralizedModeResult) -> None
        """
    
    def getDampingMatrix(self):
        pass
    
    def getGeneralizedDOFNumbering(self):
        pass
    
    def getStiffnessMatrix(self):
        pass
    
    def setDampingMatrix(self, arg0):
        pass
    
    def setGeneralizedDOFNumbering(self, arg0):
        pass
    
    def setStiffnessMatrix(self, *args, **kwargs):
        """Overloaded function.
        
        1. setStiffnessMatrix(self: libaster.GeneralizedModeResult, arg0: libaster.GeneralizedAssemblyMatrixReal) -> bool
        
        2. setStiffnessMatrix(self: libaster.GeneralizedModeResult, arg0: libaster.GeneralizedAssemblyMatrixComplex) -> bool
        """

# class ParallelMesh in libaster

class ParallelMesh(BaseMesh):
    pass
    
    # Method resolution order:
    #     ParallelMesh
    #     BaseMesh
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __init__(self, *args, **kwargs):
        """Overloaded function.
        
        1. __init__(self: libaster.ParallelMesh) -> None
        
        2. __init__(self: libaster.ParallelMesh, arg0: str) -> None
        """
    
    def getCells(self, group_name= ''):
        """Return the list of the indexes of the cells that belong to a group of cells.
        
        Arguments:
            group_name (str): Name of the local group.
        
        Returns:
            list[int]: Indexes of the cells of the local group.
        """
    
    def getCellsRank(self):
        """Return the rank of the processor which owns the cells
        
        Returns:
            list[int]: MPI-Rank of the owners of the cells
        """
    
    def getGroupsOfCells(self, local= False):
        """Return the list of the existing (local or global) groups of cells.
        
        Arguments:
            local=false (bool): search in local or global groups
        
        Returns:
            list[str]: List of (local or global) groups names (stripped).
        """
    
    def getGroupsOfNodes(self, local= False):
        """Return the list of the existing (local or global) groups of nodes.
        
        Arguments:
            local=false (bool): search in local or global groups
        
        Returns:
            list[str]: List of (local or global) groups names (stripped).
        """
    
    def getInnerCells(self):
        """Return the list of the indexes of the inner cells in the mesh
        
        Returns:
            list[int]: Indexes of the cells.
        """
    
    def getInnerNodes(self):
        """Return the list of the indexes of the inner nodes in the mesh
        
        Returns:
            list[int]: Indexes of the nodes.
        """
    
    def getNodesRank(self):
        """Return the rank of the processor which owns the nodes
        
        Returns:
            list[int]: MPI-Rank of the owners of the nodes
        """
    
    def getOuterCells(self):
        """Return the list of the indexes of the outer cells in the mesh
        
        Returns:
            list[int]: Indexes of the cells.
        """
    
    def getOuterNodes(self):
        """Return the list of the indexes of the outer nodes in the mesh
        
        Returns:
            list[int]: Indexes of the nodes.
        """
    
    def hasGroupOfCells(self, group_name, local= False):
        """The global group exists in the mesh
        
        Arguments:
            group_name (str): Name of the global group.
            local=false (bool): search in local or global groups
        
        Returns:
            bool: *True* if exists, *False* otherwise.
        """
    
    def hasGroupOfNodes(self, group_name, local= False):
        """The (local or global) group exists in the mesh
        
        Arguments:
            group_name (str): Name of the (local or global) group.
            local=false (bool): search local or global groups
        
        Returns:
            bool: *True* if exists, *False* otherwise.
        """

# class ParallelDOFNumbering in libaster

class ParallelDOFNumbering(BaseDOFNumbering):
    pass
    
    # Method resolution order:
    #     ParallelDOFNumbering
    #     BaseDOFNumbering
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __init__(self, *args, **kwargs):
        """Overloaded function.
        
        1. __init__(self: libaster.ParallelDOFNumbering) -> None
        
        2. __init__(self: libaster.ParallelDOFNumbering, arg0: str) -> None
        """
    
    def getComponentAssociatedToRow(self, row, local= False):
        """Returns the component name associated to a dof index.
        
        - If the row is associated to a physical DOF, the name of the component is returned.
        
        - If the row is associated to a Lagrange multiplier DOF for a Dirichlet boundary
          condition, the name of the component which is constrained by the multiplier is
          returned, precedeed by 'LAGR:', e.g. 'LAGR:DX'.
        
        - If the row is associated to a Lagrange multiplier DOF for a multipoint-constraint
          (MPC) implying several DOF, 'LAGR:MPC' is returned (since no component can be
          identified).
        
        Arguments:
            node (int): Index of the node.
            local (bool): row in local or global numbering
        
        Returns:
            str: component names.
        """
    
    def getComponents(self):
        """Returns all the component names assigned in the numbering.
        
        Returns:
            str: component names.
        """
    
    def getComponentsAssociatedToNode(self, node, local= False):
        """Returns the components name associated to a node index.
        
        Arguments:
            node (int): Index of the node.
            local (bool): local or parallel request
        
        Returns:
            str: component names.
        """
    
    def getGhostRows(self, local= True):
        """Returns the indexes of the ghost DOFs.
        
        Arguments:
            local (bool): local or global numbering
        
        Returns:
            int: indexes of the ghost DOFs.
        """
    
    def getNodeAssociatedToRow(self, row, local= False):
        """Returns the node index associated to a dof index.
        
        Arguments:
            row (int): Index of the dof.
            local (bool, optional): not used (default: false).
        
        Returns:
            int: index of the dof.
        """
    
    def getNumberOfDofs(self, local= False):
        """Returns the number of DOFs.
        
        Arguments:
            local (bool): local or parallel request
        
        Returns:
            int: number of DOFs.
        """
    
    def getRowAssociatedToNodeComponent(self, node, component, local= False):
        """Returns the index of the dof associated to a node.
        
        Arguments:
            node (int): Index of the node.
            component (str): name of the component
            local (bool, optional): not used (default: false).
        
        Returns:
            int: index of the dof.
        """
    
    def getRowsAssociatedToLagrangeMultipliers(self, local= False):
        """Returns the indexes of the Lagrange multipliers dof.
        
        Arguments:
            local (bool, optional): not used (default: false).
        
        Returns:
            int: indexes of the Lagrange multipliers dof.
        """
    
    def getRowsAssociatedToPhysicalDofs(self, local= False):
        """Returns the indexes of the physical dof.
        
        Arguments:
            local (bool, optional): not used (default: false).
        
        Returns:
            int: indexes of the physical dof.
        """
    
    def globalToLocalRow(self, glob):
        """Returns the local number of a global DOF.
        
        Arguments:
            glob (int): global DOF number
        
        Returns:
            int: local number of the DOF.
        """
    
    def isRowAssociatedToPhysical(self, row, local= False):
        """If the row is associated to a physical DOF, return True
        
        If the row is associated to a Lagrange multiplier DOF for a Dirichlet boundary
          condition, return False
        
        Arguments:
            row (int): Index of the dof.
            local (bool, optional): not used (default: false).
        
        Returns:
            int: index of the dof.
        """
    
    def localToGlobalRow(self, loc):
        """Returns the global number of a local DOF.
        
        Arguments:
            loc (int): local DOF number
        
        Returns:
            int: global number of the DOF.
        """
    
    def useLagrangeMultipliers(self):
        """Lagrange multipliers are used for BC or MPC.
        
        Returns:
            bool: *True* if used, *False* otherwise.
        """
    
    def useSingleLagrangeMultipliers(self):
        """Single Lagrange multipliers are used for BC or MPC.
        
        Returns:
            bool: *True* if used, *False* otherwise.
        """

# class ParallelMechanicalLoadReal in libaster

class ParallelMechanicalLoadReal(DataStructure):
    pass
    
    # Method resolution order:
    #     ParallelMechanicalLoadReal
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __init__(self, *args, **kwargs):
        """Overloaded function.
        
        1. __init__(self: libaster.ParallelMechanicalLoadReal, arg0: libaster.MechanicalLoadReal, arg1: libaster.Model) -> None
        
        2. __init__(self: libaster.ParallelMechanicalLoadReal, arg0: str, arg1: libaster.MechanicalLoadReal, arg2: libaster.Model) -> None
        """
    
    def getFiniteElementDescriptor(self):
        pass
    
    def getModel(self):
        pass

# class ParallelMechanicalLoadFunction in libaster

class ParallelMechanicalLoadFunction(DataStructure):
    pass
    
    # Method resolution order:
    #     ParallelMechanicalLoadFunction
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __init__(self, *args, **kwargs):
        """Overloaded function.
        
        1. __init__(self: libaster.ParallelMechanicalLoadFunction, arg0: libaster.MechanicalLoadFunction, arg1: libaster.Model) -> None
        
        2. __init__(self: libaster.ParallelMechanicalLoadFunction, arg0: str, arg1: libaster.MechanicalLoadFunction, arg2: libaster.Model) -> None
        """
    
    def getFiniteElementDescriptor(self):
        pass
    
    def getModel(self):
        pass

# class ParallelThermalLoadReal in libaster

class ParallelThermalLoadReal(DataStructure):
    pass
    
    # Method resolution order:
    #     ParallelThermalLoadReal
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __init__(self, *args, **kwargs):
        """Overloaded function.
        
        1. __init__(self: libaster.ParallelThermalLoadReal, arg0: libaster.ThermalLoadReal, arg1: libaster.Model) -> None
        
        2. __init__(self: libaster.ParallelThermalLoadReal, arg0: str, arg1: libaster.ThermalLoadReal, arg2: libaster.Model) -> None
        """
    
    def getFiniteElementDescriptor(self):
        pass
    
    def getModel(self):
        pass

# class ParallelThermalLoadFunction in libaster

class ParallelThermalLoadFunction(DataStructure):
    pass
    
    # Method resolution order:
    #     ParallelThermalLoadFunction
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __init__(self, *args, **kwargs):
        """Overloaded function.
        
        1. __init__(self: libaster.ParallelThermalLoadFunction, arg0: libaster.ThermalLoadFunction, arg1: libaster.Model) -> None
        
        2. __init__(self: libaster.ParallelThermalLoadFunction, arg0: str, arg1: libaster.ThermalLoadFunction, arg2: libaster.Model) -> None
        """
    
    def getFiniteElementDescriptor(self):
        pass
    
    def getModel(self):
        pass

# class ParallelFiniteElementDescriptor in libaster

class ParallelFiniteElementDescriptor(FiniteElementDescriptor):
    pass
    
    # Method resolution order:
    #     ParallelFiniteElementDescriptor
    #     FiniteElementDescriptor
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __init__(self,  *args, **kwargs):
        """Initialize self.  See help(type(self)) for accurate signature.
        """
    
    def getJoins(self):
        """Return the vector of joins between the curent domain and the others subdomains.
        
        Returns:
            list: joins between subdomains.
        """

# class ConnectionMesh in libaster

class ConnectionMesh(BaseMesh):
    pass
    
    # Method resolution order:
    #     ConnectionMesh
    #     BaseMesh
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __init__(self, *args, **kwargs):
        """Overloaded function.
        
        1. __init__(self: libaster.ConnectionMesh, arg0: libaster.ParallelMesh, arg1: List[str], arg2: List[str]) -> None
        
        2. __init__(self: libaster.ConnectionMesh, arg0: str, arg1: libaster.ParallelMesh, arg2: List[str], arg3: List[str]) -> None
        """
    
    def getCells(self, group_name= ''):
        """Return the list of the indexes of the cells that belong to a group of cells.
        
        Arguments:
            group_name (str): Name of the local group.
        
        Returns:
            list[int]: Indexes of the cells of the local group.
        """
    
    def getGroupsOfCells(self, local= False):
        """Return the list of the existing groups of cells.
        
        Returns:
            list[str]: List of groups names (stripped).
        """
    
    def getGroupsOfNodes(self, local= False):
        """Return the list of the existing groups of nodes.
        
        Returns:
            list[str]: List of groups names (stripped).
        """
    
    def getNodesGlobalNumbering(self):
        """Return a tuple of the nodes of the mesh with a global numbering
        
        Returns:
            tuple[int]: list of nodes with global numbering
        """
    
    def getNodesLocalNumbering(self):
        """Return a tuple of the nodes of the mesh with a local numbering.
        The local numbering is the one coming from the owner of the node,
        hence some nodes can have the same local numbering
        
        Returns:
            tuple[int]: list of nodes with local numbering
        """
    
    def getParallelMesh(self):
        """Return a pointer to the ParallelMesh used to built it.
        
        Returns:
            ParallelMeshPtr: pointer to the ParallelMesh
        """
    
    def hasGroupOfCells(self, name, local= False):
        """Allows to know if the given group of cells is present in the mesh
        
        Arguments:
            name (str): name of the group of cell
        
        Returns:
            bool: True if the group is present
        """
    
    def hasGroupOfNodes(self, name, local= False):
        """Allows to know if the given group of nodes is present in the mesh
        
        Arguments:
            name (str): name of the group of nodes
        
        Returns:
            bool: True if the group is present
        """

# class ResultNaming in libaster

class ResultNaming:
    pass
    
    # Method resolution order:
    #     ResultNaming
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __init__(self,  *args, **kwargs):
        """Initialize self.  See help(type(self)) for accurate signature.
        """
    
    #----------------------------------------------------------------------
    # Static methods defined here:

# class ListOfFloats in libaster

class ListOfFloats(DataStructure):
    pass
    
    # Method resolution order:
    #     ListOfFloats
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __init__(self, *args, **kwargs):
        """Overloaded function.
        
        1. __init__(self: libaster.ListOfFloats) -> None
        
        2. __init__(self: libaster.ListOfFloats, arg0: str) -> None
        """
    
    def getValues(self):
        pass
    
    def setVectorValues(self, arg0):
        pass
    
    #----------------------------------------------------------------------
    # Data descriptors defined here:
    
    @property
    def size(self):
        pass

# class ListOfIntegers in libaster

class ListOfIntegers(DataStructure):
    pass
    
    # Method resolution order:
    #     ListOfIntegers
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __init__(self, *args, **kwargs):
        """Overloaded function.
        
        1. __init__(self: libaster.ListOfIntegers) -> None
        
        2. __init__(self: libaster.ListOfIntegers, arg0: str) -> None
        """
    
    def getValues(self):
        pass
    
    def setVectorValues(self, arg0):
        pass
    
    #----------------------------------------------------------------------
    # Data descriptors defined here:
    
    @property
    def size(self):
        pass

# class EmpiricalModeResult in libaster

class EmpiricalModeResult(Result):
    pass
    
    # Method resolution order:
    #     EmpiricalModeResult
    #     Result
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __init__(self, *args, **kwargs):
        """Overloaded function.
        
        1. __init__(self: libaster.EmpiricalModeResult) -> None
        
        2. __init__(self: libaster.EmpiricalModeResult, arg0: str) -> None
        """

# class EvolutionParameter in libaster

class EvolutionParameter:
    pass
    
    # Method resolution order:
    #     EvolutionParameter
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __init__(self, result, fieldName):
        """Constructor of object
        
        Arguments:
            result (TransientResult): transient result to define external state variable
            fieldName (str): field in transient result to define external state variable
        """
    
    def setLeftExtension(self, typeExtension):
        """Set type of the extension to the left of the function to shift the results
        
        Arguments:
            typeExtension (str): type of extension ('CONSTANT', 'EXCLU', 'LINEAIRE')
        """
    
    def setRightExtension(self, typeExtension):
        """Set type of the extension to the right of the function to shift the results
        
        Arguments:
            typeExtension (str): type of extension ('CONSTANT', 'EXCLU', 'LINEAIRE')
        """
    
    def setTimeFunction(self, *args, **kwargs):
        """Overloaded function.
        
        1. setTimeFunction(self: libaster.EvolutionParameter, formula: libaster.Formula) -> None
        
        
                    Set function to shift results
        
                    Arguments:
                        formula (Formula): formula
                    
        
        2. setTimeFunction(self: libaster.EvolutionParameter, function: libaster.Function) -> None
        
        
                    Set function to shift results
        
                    Arguments:
                        function (Function): function
        """

# class ExternalStateVariable in libaster

class ExternalStateVariable:
    pass
    
    # Method resolution order:
    #     ExternalStateVariable
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __init__(self, *args, **kwargs):
        """Overloaded function.
        
        1. __init__(self: libaster.ExternalStateVariable, arg0: str, arg1: libaster.BaseMesh) -> None
        
        2. __init__(self: libaster.ExternalStateVariable, arg0: str, arg1: libaster.BaseMesh, arg2: str) -> None
        
        3. __init__(self: libaster.ExternalStateVariable, arg0: externVarEnumInt, arg1: libaster.BaseMesh) -> None
        
        4. __init__(self: libaster.ExternalStateVariable, arg0: externVarEnumInt, arg1: libaster.BaseMesh, arg2: str) -> None
        """
    
    def getField(self):
        """Get the field of values
        """
    
    def getTransientResult(self):
        """Get the transient result
        """
    
    def setEvolutionParameter(self, evolutionParameter):
        """Define evolution parameters for values of external state variable
        
        Arguments:
            evolutionParameter (EvolutionParameter): object EvolutionParameter to define
        """
    
    def setField(self, field):
        """Define constant value in time for external state variable
        
        Arguments:
            field (field): field to define value
        """
    
    def setReferenceValue(self, value):
        """Set reference value for external state variable
        
        Arguments:
            value (float): reference value
        """

# class ExternalStateVariablesResult in libaster

class ExternalStateVariablesResult(TransientResult):
    pass
    
    # Method resolution order:
    #     ExternalStateVariablesResult
    #     TransientResult
    #     Result
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __init__(self, *args, **kwargs):
        """Overloaded function.
        
        1. __init__(self: libaster.ExternalStateVariablesResult) -> None
        
        2. __init__(self: libaster.ExternalStateVariablesResult, arg0: str) -> None
        """

# built-in function createEnthalpy in libaster

def createEnthalpy(rho_cp_func, beta_func):
    """Integrate the rho_cp function by adding a point at T=0 K to be sure \\
    to always manipulate a positive enthalpy.
    
    Arguments:
        rhoc_cp_func[Function]: Function of RHO_CP
        beta_func[Function]: Function of BETA to modify (add value at T=0K)
    """

# built-in function petscFinalize in libaster

def petscFinalize():
    """Stops the PETSc interface.
    """

# built-in function _petscInitializeWithOptions in libaster

def _petscInitializeWithOptions(options):
    """Starts the PETSc interface with options.
    
    Arguments:
        options[str]: PETSc options
    """

# built-in function assemblyMatrixToPetsc in libaster

def assemblyMatrixToPetsc(*args, **kwargs):
    """Overloaded function.
    
    1. assemblyMatrixToPetsc(matr: libaster.AssemblyMatrixDisplacementReal) -> object
    
    
    Convert a *AssemblyMatrix* object to a PETSc *Mat* object.
    
    Arguments:
        matr (*AssemblyMatrix*): code_aster matrix.
    
    Returns:
        *Mat*: PETSc matrix.
            
    
    2. assemblyMatrixToPetsc(matr: libaster.AssemblyMatrixTemperatureReal) -> object
    
    
    Convert a *AssemblyMatrix* object to a PETSc *Mat* object.
    
    Arguments:
        matr (*AssemblyMatrix*): code_aster matrix.
    
    Returns:
        *Mat*: PETSc matrix.
    """

# class BehaviourProperty in libaster

class BehaviourProperty(DataStructure):
    pass
    
    # Method resolution order:
    #     BehaviourProperty
    #     DataStructure
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __init__(self, *args, **kwargs):
        """Overloaded function.
        
        1. __init__(self: libaster.BehaviourProperty) -> None
        
        2. __init__(self: libaster.BehaviourProperty, arg0: libaster.Model, arg1: libaster.MaterialField) -> None
        """
    
    def getBehaviourField(self):
        """Return a pointer to the field for behaviour.
        
        Returns:
            ConstantFieldOnCellsChar16Ptr: behaviour.
        """
    
    def getConvergenceCriteria(self):
        """Return a pointer to the field for convergence criteria.
        
        Returns:
            ConstantFieldOnCellsRealPtr: convergence criteria.
        """
    
    def getMaterialField(self):
        """Return a pointer to the material field.
        
        Returns:
            MaterialFieldPtr: material field setted.
        """
    
    def getModel(self):
        """Return a pointer to the model.
        
        Returns:
            ModelPtr: model setted.
        """
    
    def getMultipleBehaviourField(self):
        """Return a pointer to the field for multiple behaviour like cristals.
        
        Returns:
            ConstantFieldOnCellsChar16Ptr: multiple behaviour.
        """

# class CodedMaterial in libaster

class CodedMaterial:
    pass
    
    # Method resolution order:
    #     CodedMaterial
    #     pybind11_builtins.pybind11_object
    #     builtins.object
    
    # Methods defined here:
    
    def __init__(self, *args, **kwargs):
        """Overloaded function.
        
        1. __init__(self: libaster.CodedMaterial, arg0: libaster.MaterialField, arg1: libaster.Model) -> None
        
        2. __init__(self: libaster.CodedMaterial, arg0: str, arg1: libaster.MaterialField, arg2: libaster.Model) -> None
        """
    
    def allocate(self, force= False):
        pass
    
    def constant(self):
        pass
    
    def getCodedMaterialField(self):
        pass

# built-in function setFortranLoggingLevel in libaster

def setFortranLoggingLevel(level):
    """Set level of logging for fortran code.
    
    Arguments:
        level[int]: Level of logging
    """

# built-in function resetFortranLoggingLevel in libaster

def resetFortranLoggingLevel():
    """Reset level of logging for fortran code (level = 0).
    """
