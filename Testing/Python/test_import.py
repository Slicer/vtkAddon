def test_basic_import():
    from vtk import vtkAddon

    # Make sure it has one of the attributes we are expecting
    assert hasattr(vtkAddon, 'vtkCurveGenerator')
