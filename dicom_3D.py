import vtk

def ShowDicomVtk3D(dicompath, output_stl_path):
    render = vtk.vtkRenderer()
    renWin = vtk.vtkRenderWindow()
    ir = vtk.vtkRenderWindowInteractor()
    ir.SetRenderWindow(renWin)
    renWin.AddRenderer(render)

    style = vtk.vtkInteractorStyleTrackballCamera()
    ir.SetInteractorStyle(style)

    reader = vtk.vtkDICOMImageReader()
    reader.SetDirectoryName(dicompath)

    contourfilter = vtk.vtkContourFilter()
    contourfilter.SetInputConnection(reader.GetOutputPort())
    contourfilter.SetValue(0, 300)

    normal = vtk.vtkPolyDataNormals()
    normal.SetInputConnection(contourfilter.GetOutputPort())
    normal.SetFeatureAngle(60)

    conMapper = vtk.vtkPolyDataMapper()
    conMapper.SetInputConnection(normal.GetOutputPort())
    conMapper.ScalarVisibilityOff()

    conActor = vtk.vtkActor()
    conActor.SetMapper(conMapper)
    render.AddActor(conActor)

    boxFilter = vtk.vtkOutlineFilter()
    boxFilter.SetInputConnection(reader.GetOutputPort())

    boxMapper = vtk.vtkPolyDataMapper()
    boxMapper.SetInputConnection(boxFilter.GetOutputPort())

    boxActor = vtk.vtkActor()
    boxActor.SetMapper(boxMapper)
    boxActor.GetProperty().SetColor(255, 255, 255)
    render.AddActor(boxActor)

    camera = vtk.vtkCamera()
    camera.SetViewUp(0, 0, -1)
    camera.SetPosition(0, 1, 0)
    camera.SetFocalPoint(0, 0, 0)
    camera.ComputeViewPlaneNormal()
    camera.Dolly(1.5)

    render.SetActiveCamera(camera)
    render.ResetCamera()

    # 设置 conActor 的颜色
    conActor.GetProperty().SetColor(1.0, 0.0, 0.0)  # 使用 RGB 值表示颜色，这里是红色

    render.AddActor(conActor)

    # 设置 boxActor 的颜色
    boxActor.GetProperty().SetColor(0.0, 1.0, 0.0)  # 使用 RGB 值表示颜色，这里是绿色

    render.AddActor(boxActor)

    ir.Initialize()
    # 保存 STL 文件
    stl_writer = vtk.vtkSTLWriter()
    stl_writer.SetFileName(output_stl_path)
    stl_writer.SetInputConnection(normal.GetOutputPort())
    stl_writer.Write()

    ir.Start()


# dcm文件路径
path = "/Users/serena-mo/Documents/医学影像文件/202310_anonymous"
# STL 文件保存路径
output_stl_path = "./3D_dicom.stl"
ShowDicomVtk3D(path, output_stl_path)
