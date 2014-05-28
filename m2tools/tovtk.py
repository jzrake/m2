import command
import checkpoint
import autolog


class VtkExporterImageData(object):

    @autolog.logmethod
    def run(self, args):
        import vtk
        import numpy as np

        for filename in args.filenames:
            self.checkpoint = checkpoint.Checkpoint(filename)
            self.run_file('B', filename.replace('.h5', '.B.vtk'), True)
            self.run_file('v', filename.replace('.h5', '.v.vtk'), True)
            self.run_file('p', filename.replace('.h5', '.p.vtk'), False)
            self.run_file('d', filename.replace('.h5', '.d.vtk'), False)
            self.checkpoint.close()

    @autolog.logmethod
    def run_file(self, field, output_name, vectors=False):
        import vtk
        import numpy as np

        Nx, Ny, Nz = self.checkpoint.domain_resolution[1:]

        if vectors:
            data = np.zeros([3, Nx, Ny, Nz])
            for i in range(3):
                data[i] = self.checkpoint.cell_primitive[field + str(i+1)][:]
                data[i] -= data[i].min()
                data[i] /= data[i].max()
        else:
            data = self.checkpoint.cell_primitive[field][:]

        data_matrix = np.zeros(data.shape, dtype=np.uint8, order='F')
        data_matrix[...] = data * 256
 
        dataImporter = vtk.vtkImageImport()
        ata_string = data_matrix.tostring()

        dataImporter.CopyImportVoidPointer(data_matrix, data_matrix.nbytes)
        dataImporter.SetDataScalarTypeToUnsignedChar()
        dataImporter.SetNumberOfScalarComponents(3 if vectors else 1)
        dataImporter.SetDataExtent (0, Nx-1, 0, Ny-1, 0, Nz-1)
        dataImporter.SetWholeExtent(0, Nx-1, 0, Ny-1, 0, Nz-1)
        dataImporter.SetScalarArrayName(field)

        writer = vtk.vtkDataSetWriter()
        writer.SetInputConnection(dataImporter.GetOutputPort())
        writer.SetFileName(output_name)
        writer.SetFileTypeToBinary()
        writer.Update()
        writer.Write()



class VtkExporterStructuredGrid(object):

    @autolog.logmethod
    def run(self, args):
        for filename in args.filenames:
            self.checkpoint = checkpoint.Checkpoint(filename)
            self.run_file(filename.replace('.h5', '.vtk'))
            self.checkpoint.close()

    @autolog.logmethod
    def run_file(self, output_name):
        import itertools

        outfile = open(output_name, 'w')
        X1, X2, X3 = self.checkpoint.cell_edge_coordinates
        n1, n2, n3 = self.checkpoint.domain_resolution[1:]
        num_points = (n1+1) * (n2+1) * (n3+1)
        num_cells = n1 * n2 * n3

        header = \
"""# vtk DataFile Version 3.1 
m2 astrophysical MHD code
ASCII
DATASET STRUCTURED_GRID
DIMENSIONS %(n1)d %(n2)d %(n3)d
POINTS %(num_points)d FLOAT
"""

        outfile.write(header % dict(n1=n1+1, n2=n2+1, n3=n3+1, num_points=num_points))
        for x1, x2, x3 in itertools.product(X1, X2, X3):
            outfile.write("%+10.8e %+10.8e %+10.8e\n" % (x1, x2, x3))

        outfile.write("\nCELL_DATA %d\n" % num_cells)
        for f in ['p', 'd']:
            outfile.write("SCALARS %s FLOAT\n" % f)
            outfile.write("LOOKUP_TABLE default\n")
            for z in self.checkpoint.cell_primitive[f][:].flat:
                outfile.write("%+10.8e\n" % z)

        v1 = self.checkpoint.cell_primitive['v1'][:].flat
        v2 = self.checkpoint.cell_primitive['v2'][:].flat
        v3 = self.checkpoint.cell_primitive['v3'][:].flat
        outfile.write("\nVECTORS velocity-field FLOAT\n")
        for b1, b2, b3 in zip(v1, v2, v3):
            outfile.write("%+10.8e %+10.8e %+10.8e\n" % (b1, b2, b3))

        B1 = self.checkpoint.cell_primitive['B1'][:].flat
        B2 = self.checkpoint.cell_primitive['B2'][:].flat
        B3 = self.checkpoint.cell_primitive['B3'][:].flat
        outfile.write("\nVECTORS magnetic-field FLOAT\n")
        for b1, b2, b3 in zip(B1, B2, B3):
            outfile.write("%+10.8e %+10.8e %+10.8e\n" % (b1, b2, b3))


class VtkExporterUnstructuredGrid(object):

    @autolog.logmethod
    def run(self, args):
        for filename in args.filenames:
            self.run_file(args, filename)

    @autolog.logmethod
    def run_file(self, args, filename):
        self.checkpoint = checkpoint.Checkpoint(filename)
        self.points = self.get_points()
        self.cells = self.get_cells()
        outfile = open(filename.replace('.h5', '.vtk'), 'w')

        outfile.write("# vtk DataFile Version 3.1\n")
        outfile.write("m2: astrophysical MHD code\n")
        outfile.write("ASCII\n")
        outfile.write("DATASET UNSTRUCTURED_GRID\n\n")
        outfile.write("POINTS %d FLOAT\n" % len(self.points[1]))
        
        for p in self.points[1]:
            outfile.write("%+10.8e %+10.8e %+10.8e\n" % p)

        outfile.write("\nCELLS %d %d\n" % (len(self.cells), 9*len(self.cells)))
        for c in self.cells:
            outfile.write("8" + (" %d"*8) % c + "\n")

        outfile.write("\nCELL_TYPES %d\n" % len(self.cells))
        for c in self.cells:
            outfile.write("12\n")

        outfile.write("\nCELL_DATA %d\n" % len(self.cells))
        for field in 'pd':
            outfile.write("\nSCALARS %s FLOAT\n" % field)
            outfile.write("LOOKUP_TABLE default\n")
            for z in self.checkpoint.cell_primitive[field][:].flat:
                outfile.write("%+10.8e\n" % z)

        for field in 'Bv':
            outfile.write("\nVECTORS %s FLOAT\n" % field)
            B1 = self.checkpoint.cell_primitive[field+'1'][:].flat
            B2 = self.checkpoint.cell_primitive[field+'2'][:].flat
            B3 = self.checkpoint.cell_primitive[field+'3'][:].flat
            for b1, b2, b3 in zip(B1, B2, B3):
                outfile.write("%+10.8e %+10.8e %+10.8e\n" % (b1, b2, b3))

        self.checkpoint.close()

    @autolog.logmethod
    def get_points(self):
        import itertools
        import numpy as np
        point_num = { }
        point_val = [ ]
        X1, X2, X3 = self.checkpoint.cell_edge_coordinates
        coords = itertools.product(enumerate(X1),
                                   enumerate(X2),
                                   enumerate(X3))
        for n, ((i,x1), (j,x2), (k,x3)) in enumerate(coords):
            point_num[(i, j, k)] = n
            point_val.append((x1, x2, x3))
        return (point_num, point_val)

    @autolog.logmethod
    def get_points_spherical(self):
        """ not really functioning yet """
        import itertools
        import numpy as np
        point_num = { }
        point_val = [ ]
        X1, X2, X3 = self.checkpoint.cell_edge_meshgrid
        Y1 = X1 * np.sin(X2) * np.cos(X3)
        Y2 = X1 * np.sin(X2) * np.sin(X3)
        Y3 = X1 * np.cos(X2) * (1 + 0*X3)
        n1, n2, n3 = self.checkpoint.domain_resolution[1:]
        indices = itertools.product(range(n1+1), range(n2+1), range(n3+1))
        coords = zip(Y1.flat, Y2.flat, Y3.flat)
        for n, ((i, j, k), (x1, x2, x3)) in enumerate(zip(indices, coords)):
            point_num[(i, j, k)] = n
            point_val.append((x1, x2, x3))
        return (point_num, point_val)

    @autolog.logmethod
    def get_cells(self):
        import itertools
        cells = [ ]
        P = self.points[0]
        n1, n2, n3 = self.checkpoint.domain_resolution[1:]
        indices = itertools.product(range(n1), range(n2), range(n3))
        for i,j,k in indices:
            cells.append((P[(i+0,j+0,k+0)],
                          P[(i+1,j+0,k+0)],
                          P[(i+1,j+1,k+0)],
                          P[(i+0,j+1,k+0)],
                          P[(i+0,j+0,k+1)],
                          P[(i+1,j+0,k+1)],
                          P[(i+1,j+1,k+1)],
                          P[(i+0,j+1,k+1)]))
        return cells



class VtkExporterStructuredGridPoints(object):

    @autolog.logmethod
    def run(self, args):
        for filename in args.filenames:
            self.run_file(args, filename)

    @autolog.logmethod
    def run_file(self, args, filename):
        self.checkpoint = checkpoint.Checkpoint(filename)
        n1, n2, n3 = self.checkpoint.domain_resolution[1:]
        points = self.checkpoint.get_cell_center_coordinates()
        if self.checkpoint.geometry == 'spherical':
            slicer = checkpoint.Slicer()
            self.checkpoint.set_selection(slicer[1:,1:,0:])

        outfile = open(filename.replace('.h5', '.vtk'), 'w')
        outfile.write("# vtk DataFile Version 3.1\n")
        outfile.write("m2: astrophysical MHD code\n")
        outfile.write("ASCII\n")
        outfile.write("DATASET STRUCTURED_GRID\n")
        outfile.write("DIMENSIONS %d %d %d\n" % (n1, n2, n3))

        outfile.write("\nPOINTS %d FLOAT\n" % (n1 * n2 * n3))
        for p in points:
            outfile.write("%+10.8e %+10.8e %+10.8e\n" % p)

        outfile.write("\nPOINT_DATA %d\n" % (n1 * n2 * n3))

        for field in 'pd':
            outfile.write("\nSCALARS %s FLOAT\n" % field)
            outfile.write("LOOKUP_TABLE default\n")
            for z in self.checkpoint.get_field(field).flat:
                outfile.write("%+10.8e\n" % z)

        for field in 'Bv':
            outfile.write("\nVECTORS %s FLOAT\n" % field)
            f1, f2, f3 = self.checkpoint.get_vector_field(field)
            for b1, b2, b3 in zip(f1.flat, f2.flat, f3.flat):
                outfile.write("%+10.8e %+10.8e %+10.8e\n" % (b1, b2, b3))
        self.checkpoint.close()



class VtkExporterStructuredGridCells(object):

    @autolog.logmethod
    def run(self, args):
        for filename in args.filenames:
            self.run_file(args, filename)

    @autolog.logmethod
    def run_file(self, args, filename):
        self.checkpoint = checkpoint.Checkpoint(filename)
        n1, n2, n3 = self.checkpoint.domain_resolution[1:]
        points = self.checkpoint.get_cell_edge_coordinates()
        if self.checkpoint.geometry == 'spherical':
            slicer = checkpoint.Slicer()
            self.checkpoint.set_selection(slicer[1:,1:,0:])

        outfile = open(filename.replace('.h5', '.vtk'), 'w')
        outfile.write("# vtk DataFile Version 3.1\n")
        outfile.write("m2: astrophysical MHD code\n")
        outfile.write("ASCII\n")
        outfile.write("DATASET STRUCTURED_GRID\n")
        outfile.write("DIMENSIONS %d %d %d\n" % (n1+1, n2+1, n3+1))

        outfile.write("\nPOINTS %d FLOAT\n" % ((n1+1)*(n2+1)*(n3+1)))
        for p in points:
            outfile.write("%+10.8e %+10.8e %+10.8e\n" % p)
        outfile.write("\nCELL_DATA %d\n" % (n1 * n2 * n3))

        for field in 'pd':
            outfile.write("\nSCALARS %s FLOAT\n" % field)
            outfile.write("LOOKUP_TABLE default\n")
            for z in self.checkpoint.get_field(field).flat:
                outfile.write("%+10.8e\n" % z)

        for field in 'Bv':
            outfile.write("\nVECTORS %s FLOAT\n" % field)
            f1, f2, f3 = self.checkpoint.get_vector_field(field)
            for b1, b2, b3 in zip(f1.flat, f2.flat, f3.flat):
                outfile.write("%+10.8e %+10.8e %+10.8e\n" % (b1, b2, b3))
        self.checkpoint.close()



class VtkDataExporterCommand(command.Command):
    _alias = 'tovtk'
    _exporters = {'structured-cells': VtkExporterStructuredGridCells,
                  'structured-points': VtkExporterStructuredGridPoints,
                  'structured': VtkExporterStructuredGrid,
                  'unstructured': VtkExporterUnstructuredGrid,
                  'image': VtkExporterImageData}

    def configure_parser(self, parser):
        parser.add_argument("filenames", nargs='+')
        parser.add_argument("--exporter",
                            default='structured-cells',
                            choices=self._exporters.keys())

    def run(self, args):
        exporter = self._exporters[args.exporter]()
        exporter.run(args) 


