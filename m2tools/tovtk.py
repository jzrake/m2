import command
import checkpoint
import autolog


class VtkDataExporterVtk(object):

    @autolog.logmethod
    def run(self, args):
        import vtk
        import numpy as np

        for filename in args.filenames:
            self.checkpoint = checkpoint.Checkpoint(filename)
            self.output_file('B', filename.replace('.h5', '.B.vtk'), True)
            self.output_file('v', filename.replace('.h5', '.v.vtk'), True)
            self.output_file('p', filename.replace('.h5', '.p.vtk'), False)
            self.output_file('d', filename.replace('.h5', '.d.vtk'), False)
            self.checkpoint.close()

    @autolog.logmethod
    def output_file(self, field, output_name, vectors=False):
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



class VtkDataExporterByhand(object):

    @autolog.logmethod
    def run(self, args):
        for filename in args.filenames:
            self.checkpoint = checkpoint.Checkpoint(filename)
            self.output_file(filename.replace('.h5', '.vtk'))
            self.checkpoint.close()

    @autolog.logmethod
    def output_file(self, output_name):
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



class VtkDataExporterCommand(command.Command):
    _alias = 'tovtk'

    def configure_parser(self, parser):
        parser.add_argument("filenames", nargs='+')
        parser.add_argument("--byhand", action='store_true')

    def run(self, args):
        if args.byhand:
            exporter = VtkDataExporterByhand()
        else:
            exporter = VtkDataExporterVtk()
        exporter.run(args) 


