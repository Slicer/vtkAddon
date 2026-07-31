#!/usr/bin/env python3
"""Prune a wrapped module's Python Init.data DEPENDS to the modules it truly needs.

Background
----------
Slicer wraps its C++ classes for Python through vtkMacroKitPythonWrap, which
seeds each module's wrapping dependencies from the whole VTK link line
(``${VTK_LIBRARIES}``). The generated ``<Module>Init.data`` therefore lists a
``DEPENDS`` on every VTK python module -- 125 of them -- even though a module
such as vtkAddon only subclasses a handful of VTK base classes.

That over-listing was harmless while ``import vtk`` eagerly loaded all of VTK.
With lazy VTK loading it is not: a wrapped subclass whose direct base class has
not been imported is built without that base and silently loses the inherited
methods (the vtkOrientedGridTransform.SetDisplacementGridData regression). The
fix restores the base-class imports -- but importing all 125 VTK modules at
startup defeats the point of lazy loading.

A wrapped subclass only needs its **direct base class** registered first; each
VTK module in turn imports its own bases transitively. So the minimal, correct
set of VTK DEPENDS is exactly the modules that own the direct base classes of
the classes this module wraps.

What this does
--------------
Reads the class list and DEPENDS from an Init.data file and the merged wrapping
hierarchy (``Class : Super ; header ; owning_module``). Keeps every
``vtkmodules.X`` line whose X owns a direct base class of a wrapped class, drops
the rest, and passes all other DEPENDS lines (cross-module Slicer dependencies
such as MRMLCorePython) through unchanged. Conservative by construction: only
spurious VTK modules are removed.
"""
import argparse
import sys


def normalize(name):
    """Reduce a hierarchy class token to a plain lookup key.

    Drops template parameters (``Foo<T>`` -> ``Foo``) and rejects typedef-style
    supers (``= something``) that are not real base classes.
    """
    name = name.strip()
    lt = name.find('<')
    if lt != -1:
        name = name[:lt]
    if name.startswith('='):
        return ''
    return name.strip()


def parse_hierarchy(path):
    """Return (super_of, module_of) mapping over every class in the file."""
    super_of = {}
    module_of = {}
    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            fields = [part.strip() for part in line.split(';')]
            if len(fields) < 3:
                continue
            decl, _header, module = fields[0], fields[1], fields[2]
            if ' : ' in decl:
                cls, sup = decl.split(' : ', 1)
            else:
                cls, sup = decl, None
            cls = normalize(cls)
            if not cls:
                continue
            module_of[cls] = module
            if sup:
                super_of[cls] = normalize(sup)
    return super_of, module_of


def read_init_data(path):
    """Split an Init.data file into (target, class_list, depends_list)."""
    with open(path) as handle:
        lines = [line.rstrip('\n') for line in handle]
    target = lines[0] if lines else ''
    classes = []
    depends = []
    in_depends = False
    for line in lines[1:]:
        if line.strip() == 'DEPENDS':
            in_depends = True
            continue
        if in_depends:
            if line.strip():
                depends.append(line.strip())
        else:
            if line.strip():
                classes.append(line.strip())
    return target, classes, depends


def base_owner_modules(classes, super_of, module_of, own_target):
    """Owning modules of the direct base classes of the wrapped classes."""
    owners = set()
    for cls in classes:
        sup = super_of.get(cls)
        if not sup:
            continue
        owner = module_of.get(sup)
        if owner and owner != own_target:
            owners.add(owner)
    return owners


def prune(init_data_path, hierarchy_path, output_path):
    target, classes, depends = read_init_data(init_data_path)
    super_of, module_of = parse_hierarchy(hierarchy_path)
    owners = base_owner_modules(classes, super_of, module_of, target)

    kept = []
    dropped = []
    for dep in depends:
        if dep.startswith('vtkmodules.'):
            module = dep[len('vtkmodules.'):]
            if module in owners:
                kept.append(dep)
            else:
                dropped.append(dep)
        else:
            # Cross-module Slicer dependency (e.g. MRMLCorePython): keep as-is.
            kept.append(dep)

    with open(output_path, 'w') as handle:
        handle.write(target + '\n')
        for cls in classes:
            handle.write(cls + '\n')
        if kept:
            handle.write('DEPENDS\n')
            for dep in kept:
                handle.write(dep + '\n')

    sys.stderr.write(
        '[prune-depends] %s: kept %d, dropped %d VTK module(s)\n'
        % (target, len(kept), len(dropped)))


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('init_data', help='input <Module>Init.data (full DEPENDS)')
    parser.add_argument('hierarchy', help='merged <Module>Hierarchy.txt')
    parser.add_argument('output', help='pruned Init.data to write')
    args = parser.parse_args(argv)
    prune(args.init_data, args.hierarchy, args.output)


if __name__ == '__main__':
    main()
