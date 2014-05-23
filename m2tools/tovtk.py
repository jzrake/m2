import command
import checkpoint
import autolog

class VtkDataExporterCommand(command.Command):
    _alias = 'tovtk'

    def configure_parser(self, parser):
        parser.add_argument("filenames", nargs='+')

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
        data_string = data_matrix.tostring()

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

