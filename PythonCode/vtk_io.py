# -*- coding: utf-8 -*-
#
'''
I/O for VTK <https://www.vtk.org/wp-content/uploads/2015/04/file-formats.pdf>.

.. moduleauthor:: Nico Schlömer <nico.schloemer@gmail.com>

Adapted from https://github.com/nschloe/meshio
'''
import logging
import numpy

# from .__about__ import __version__


# https://www.vtk.org/doc/nightly/html/vtkCellType_8h_source.html
vtk_to_meshio_type = {
    0: 'empty',
    1: 'vertex',
    # 2: 'poly_vertex',
    3: 'line',
    # 4: 'poly_line',
    5: 'triangle',
    # 6: 'triangle_strip',
    # 7: 'polygon',
    # 8: 'pixel',
    9: 'quad',
    10: 'tetra',
    # 11: 'voxel',
    12: 'hexahedron',
    13: 'wedge',
    14: 'pyramid',
    15: 'penta_prism',
    16: 'hexa_prism',
    21: 'line3',
    22: 'triangle6',
    23: 'quad8',
    24: 'tetra10',
    25: 'hexahedron20',
    26: 'wedge15',
    27: 'pyramid13',
    28: 'quad9',
    29: 'hexahedron27',
    30: 'quad6',
    31: 'wedge12',
    32: 'wedge18',
    33: 'hexahedron24',
    34: 'triangle7',
    35: 'line4',
    #
    # 60: VTK_HIGHER_ORDER_EDGE,
    # 61: VTK_HIGHER_ORDER_TRIANGLE,
    # 62: VTK_HIGHER_ORDER_QUAD,
    # 63: VTK_HIGHER_ORDER_POLYGON,
    # 64: VTK_HIGHER_ORDER_TETRAHEDRON,
    # 65: VTK_HIGHER_ORDER_WEDGE,
    # 66: VTK_HIGHER_ORDER_PYRAMID,
    # 67: VTK_HIGHER_ORDER_HEXAHEDRON,
    }
meshio_to_vtk_type = {v: k for k, v in vtk_to_meshio_type.items()}


# These are all VTK data types. One sometimes finds 'vtktypeint64', but
# this is ill-formed.
vtk_to_numpy_dtype = {
    'bit': numpy.dtype('bool'),
    'unsigned_char': numpy.dtype('uint8'),
    'char': numpy.dtype('int8'),
    'unsigned_short': numpy.dtype('uint16'),
    'short': numpy.dtype('int16'),
    'unsigned_int': numpy.dtype('uint32'),
    'int': numpy.dtype('int32'),
    'unsigned_long': numpy.dtype('int64'),
    'long': numpy.dtype('int64'),
    'float': numpy.dtype('float32'),
    'double': numpy.dtype('float64'),
    }

numpy_to_vtk_dtype = {v: k for k, v in vtk_to_numpy_dtype.items()}


def write(filename,
          points,
          cells,
          point_data=None,
          cell_data=None,
          field_data=None,
          write_binary=True
          ):
    if not write_binary:
        logging.warning('VTK ASCII files are only meant for debugging.')

    point_data = {} if point_data is None else point_data
    cell_data = {} if cell_data is None else cell_data
    field_data = {} if field_data is None else field_data

    with open(filename, 'wb') as f:
        f.write('# vtk DataFile Version 4.0\n'.encode('utf-8'))
        f.write('vtk output\n'.encode('utf-8'))
        # f.write('written Sby meshio v{}\n'.format(__version__).encode('utf-8'))
        f.write(('BINARY\n' if write_binary else 'ASCII\n').encode('utf-8'))
        f.write('DATASET UNSTRUCTURED_GRID\n'.encode('utf-8'))

        # write points and cells
        _write_points(f, points, write_binary)
        _write_cells(f, cells, write_binary)

        # write point data
        if point_data:
            num_points = len(points)
            f.write('POINT_DATA {}\n'.format(num_points).encode('utf-8'))
            _write_field_data(f, point_data, write_binary)

        # write cell data
        if cell_data:
            total_num_cells = sum([len(c) for c in cells.values()])
            cell_data_raw = raw_from_cell_data(cell_data)
            f.write('CELL_DATA {}\n'.format(total_num_cells).encode('utf-8'))
            _write_field_data(f, cell_data_raw, write_binary)

    return


def _write_points(f, points, write_binary):
    f.write(
        'POINTS {} {}\n'.format(
            len(points), numpy_to_vtk_dtype[points.dtype]
            ).encode('utf-8'))

    if write_binary:
        # Binary data must be big endian, see
        # <https://www.vtk.org/Wiki/VTK/Writing_VTK_files_using_python#.22legacy.22>.
        points.astype(points.dtype.newbyteorder('>')).tofile(f, sep='')
    else:
        # ascii
        points.tofile(f, sep=' ')
    f.write('\n'.encode('utf-8'))
    return


def _write_cells(f, cells, write_binary):
    total_num_cells = sum([len(c) for c in cells.values()])
    total_num_idx = sum([numpy.prod(c.shape) for c in cells.values()])
    # For each cell, the number of nodes is stored
    total_num_idx += total_num_cells
    f.write(
        'CELLS {} {}\n'.format(total_num_cells, total_num_idx)
        .encode('utf-8'))
    if write_binary:
        for key in cells:
            n = cells[key].shape[1]
            d = numpy.column_stack([
                numpy.full(len(cells[key]), n), cells[key]
                ]).astype(numpy.dtype('>i4'))
            f.write(d.tostring())
        if write_binary:
            f.write('\n'.encode('utf-8'))
    else:
        # ascii
        for key in cells:
            n = cells[key].shape[1]
            for cell in cells[key]:
                f.write((' '.join([
                    '{}'.format(idx)
                    for idx in numpy.concatenate([[n], cell])
                    ]) + '\n').encode('utf-8'))

    # write cell types
    f.write('CELL_TYPES {}\n'.format(total_num_cells).encode('utf-8'))
    if write_binary:
        for key in cells:
            d = numpy.full(
                len(cells[key]), meshio_to_vtk_type[key]
                ).astype(numpy.dtype('>i4'))
            f.write(d.tostring())
        f.write('\n'.encode('utf-8'))
    else:
        # ascii
        for key in cells:
            for _ in range(len(cells[key])):
                f.write(
                    '{}\n'.format(meshio_to_vtk_type[key]).encode('utf-8')
                    )
    return


def _write_field_data(f, data, write_binary):
    f.write((
        'FIELD FieldData {}\n'.format(len(data))
        ).encode('utf-8'))
    for name, values in data.items():
        if len(values.shape) == 1:
            num_tuples = values.shape[0]
            num_components = 1
        else:
            assert len(values.shape) == 2, \
                'Only one and two-dimensional field data supported.'
            num_tuples = values.shape[0]
            num_components = values.shape[1]
        f.write(('{} {} {} {}\n'.format(
            name, num_components, num_tuples,
            numpy_to_vtk_dtype[values.dtype]
            )).encode('utf-8'))
        if write_binary:
            values.astype(values.dtype.newbyteorder('>')).tofile(f, sep='')
        else:
            # ascii
            values.tofile(f, sep=' ')
            # numpy.savetxt(f, points)
        f.write('\n'.encode('utf-8'))
    return