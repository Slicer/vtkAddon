# vtkAddon

## Overview

This module contains general-purpose features that may be integrated into VTK library in the future.

## Contact

Questions regarding this project should be posted on 3D Slicer forum: https://discourse.slicer.org

## History

The original vtkAddon library was created and integrated into Slicer by Andras Lasso (PerkLab, Queen's University)
in March 2015, and maintained by the Slicer community since then.

In July 2015, Bill Lorensen copied the vtkAddon sources into a [dedicated repository][lorensen-vtkAddon] and added
support to build the library as a VTK remote module. However, since this repository duplicated content already
maintained in the Slicer source tree, it was eventually deprecated and left unmaintained.

Finally, in March 2020, this GitHub project was created by extracting the history from the main Slicer repository
and the project was [integrated][slicer-vtkAddon-pr] into Slicer as a remote project, enabling consistent updates
while decoupling development from the core Slicer source tree.

### VTK module Python wrapping

The module `vtkMacroKitPythonWrap` was first developed by Jean-Christophe Fillion-Robin (Kitware) within the Slicer
repository starting in 2011 (See [CMake/vtkMacroKitPythonWrap.cmake][slicer-vtkMacroKitPythonWrap] history). Later
in 2015, this module was moved into the dedicated `vtkAddon` repository  to support use in other projects.

In September 2020, support for VTK versions earlier than 8.2 was removed ([a64d854](https://github.com/Slicer/vtkAddon/commit/a64d854)),
and compatibility with VTK 8.90 and newer was added ([4a3a79b3](https://github.com/Slicer/vtkAddon/pull/12/commits/4a3a79b3d4ffbc39a33a6803f8da80a94bb270c6))
to align with modern VTK development practices.

## License

This project is distributed under the BSD-style Slicer license allowing academic and commercial use without any restrictions. Please see the [LICENSE](LICENSE) file for details.

[lorensen-vtkAddon]: https://github.com/lorensen/vtkAddon
[slicer-vtkAddon-pr]: https://github.com/Slicer/Slicer/pull/4765
[slicer-vtkMacroKitPythonWrap]: https://github.com/Slicer/Slicer/commits/5b3e3c8b3eed806bcfdf8df57aefdd6d8ce52fac/CMake/vtkMacroKitPythonWrap.cmake
