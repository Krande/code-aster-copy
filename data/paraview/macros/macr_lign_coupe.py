from __future__ import annotations
from packaging.version import Version
import enum

# state file generated using paraview version 5.11.2
import paraview

PV_VERSION = Version(paraview.__version__)
paraview.compatibility.major = 5
paraview.compatibility.minor = 11

#### import the simple module from the paraview
from paraview.simple import (
    Glyph,
    GetActiveSource,
    GetScalarBar,
    Calculator,
    GetActiveView,
    GetColorTransferFunction,
    ColorBy,
    HideScalarBarIfNotNeeded,
    Render,
    Show,
    Hide,
    UpdatePipeline,
)

#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()


PREFIX_CUT = "CUT_"


class AxisFields(enum.Enum):
    VECT_X = f"{PREFIX_CUT}AXE_X"
    VECT_Y = f"{PREFIX_CUT}AXE_Y"
    VECT_Z = f"{PREFIX_CUT}AXE_Z"


AXIS_LOCATION_MASK = f"{PREFIX_CUT}AXE_MASK"

AXIS_TO_COLOR = {
    AxisFields.VECT_X: [1.0, 0.0, 0.0],
    AxisFields.VECT_Y: [1.0, 1.0, 0.0],
    AxisFields.VECT_Z: [0.0, 1.0, 0.0],
}
WHITE = [1.0, 1.0, 1.0]
BLACK = [.0, .0, .0]


def get_available_point_arrays(source) -> dict[str, "PointData"]:
    return {point_data.Name: point_data for point_data in source.PointData if point_data.Name.startswith(PREFIX_CUT)}


def get_available_components(point_fields_data, field_name: str) -> list[str]:
    ai = point_fields_data.GetArrayInformation(field_name)
    return [ai.GetComponentName(compo_num) for compo_num in range(ai.GetNumberOfComponents())]


def main():
    base_source = GetActiveSource()
    view = GetActiveView()
    available_points_arrays = get_available_point_arrays(source=base_source)
    available_axis_arrays = tuple((field.value for field in AxisFields if field.value in available_points_arrays))
    available_non_axis_arrays = list((field for field in available_points_arrays if field not in available_axis_arrays))

    if AXIS_LOCATION_MASK in available_non_axis_arrays:
        available_non_axis_arrays.remove(AXIS_LOCATION_MASK)

    xmin, xmax, ymin, ymax, zmin, zmax = base_source.GetDataInformation().DataInformation.GetBounds()
    max_bound = max((xmax-xmin, ymax - ymin, zmax - zmin))

    point_fields_data = base_source.PointData.GetFieldData()

    parent_filter = base_source

    line_source = Calculator(registrationName="LINES", Input=parent_filter)
    line_display = Show(line_source, view, "UnstructuredGridRepresentation")
    line_display.SetRepresentationType("Wireframe")
    line_display.RenderLinesAsTubes = 1
    line_display.LineWidth = 10.0
    # change solid color
    line_display.AmbientColor = WHITE
    line_display.DiffuseColor = WHITE
    point_source = Calculator(registrationName="POINTS", Input=parent_filter)
    point_display = Show(point_source, view, "UnstructuredGridRepresentation")
    point_display.SetRepresentationType("Points")
    point_display.RenderPointsAsSpheres = 1
    point_display.PointSize = 15.0
    point_display.AmbientColor = BLACK
    point_display.DiffuseColor = BLACK

    # TODO add Point and Line source to enhance visualisation
    if available_axis_arrays:
        glyph_scale = Calculator(registrationName="SET GLYPH SIZE", Input=parent_filter)

        glyph_scale.Function = f"{max_bound / 10.} * {AXIS_LOCATION_MASK}"
        glyph_scale.ResultArrayName = "GLYPH_SCALE"

        parent_filter = glyph_scale

    for field_name in available_axis_arrays:
        # create a new 'Glyph'
        glyph = Glyph(
            registrationName=field_name,
            Input=parent_filter,
            GlyphType="Arrow",
        )
        glyph.OrientationArray = ["POINTS", field_name]
        glyph.GlyphTransform = "Transform2"
        glyph.GlyphMode = "All Points"
        glyph.ScaleArray = ['POINTS', 'GLYPH_SCALE']
        glyph.ScaleFactor = 1
        # get display properties
        glyph_display = Show(glyph, view, "UnstructuredGridRepresentation")
        # disable scalar coloring
        ColorBy(glyph_display, None)
        color = AXIS_TO_COLOR[AxisFields(field_name)]
        # change solid color
        glyph_display.AmbientColor = color
        glyph_display.DiffuseColor = color
        Hide(glyph, view)

    for field_name in available_non_axis_arrays:
        components = get_available_components(point_fields_data=point_fields_data, field_name=field_name)
        # create a new 'Extract Component'
        for compo_name in components:
            output_name = f"{field_name}_{compo_name}"
            extract_compo = Calculator(registrationName=output_name, Input=parent_filter)
            extract_compo.Function = output_name
            extract_compo.ResultArrayName = output_name
            extract_compo.UpdatePipeline()
            # show data in view
            extract_compo_display = Show(extract_compo, view, "UnstructuredGridRepresentation")
            if PV_VERSION >= Version("5.13.0"):
                extract_compo_display.DisableLighting = 1
            extract_compo_display.SetRepresentationType("Wireframe")
            extract_compo_display.RenderLinesAsTubes = 1
            extract_compo_display.RenderPointsAsSpheres = 1
            extract_compo_display.LineWidth = 10.0
            extract_compo_display.PointSize = 20.0

            extract_compo_display.GaussianRadius = 0.2
            extract_compo_display.ShaderPreset = "Sphere"

            # set scalar coloring
            ColorBy(extract_compo_display, ("POINTS", output_name))

            # get separate color transfer function/color map
            extract_compo_lut = GetColorTransferFunction(output_name, extract_compo_display)

            # Apply a preset using its name.
            extract_compo_lut.ApplyPreset("Turbo", True)
            extract_compo_lut.NumberOfTableValues = 12
            extract_compo_lut.RescaleOnVisibilityChange = 1
            extract_compo_lut.AutomaticRescaleRangeMode = 'Clamp and update every timestep'

            # Scalar bar settings
            colorbar = GetScalarBar(extract_compo_lut, view)
            colorbar.ScalarBarLength = 0.7
            colorbar.DrawTickMarks = 0
            colorbar.DrawTickLabels = 0
            colorbar.AddRangeLabels = 0
            colorbar.AddRangeAnnotations = 1
            colorbar.AutomaticAnnotations = 1
            colorbar.WindowLocation = 'Lower Right Corner'

            # Hide the scalar bar for this color map if no visible data is colored by it.
            HideScalarBarIfNotNeeded(extract_compo_lut, view)

            # rescale color and/or opacity maps used to include current data range
            extract_compo_display.RescaleTransferFunctionToDataRange(True)

            # show color bar/color legend
            extract_compo_display.SetScalarBarVisibility(view, True)
            Hide(extract_compo, view)

    Render(view)
    UpdatePipeline()


main()
