import command
import checkpoint
import autolog



class VtkStructuredGridExporter(object):

    @autolog.logmethod
    def run(self, args):
        for filename in args.filenames:
            self.run_file(args, filename)

    @autolog.logmethod
    def run_file(self, args, filename):
        import numpy as np

        chkpt = checkpoint.Checkpoint(filename)
        n1, n2, n3 = chkpt.domain_resolution[1:]
        p = np.zeros([n1+1, n2+1, n3+1, 3], dtype=np.float32)
        p[...,0], p[...,1], p[...,2] = chkpt.get_cell_edge_coordinates() 
        p = p.swapaxes(0, 2) # VTK reads arrays in Fortan order

        # if spherical then lower r and theta boundaries have an extra index
        if chkpt.geometry == 'spherical':
            slicer = checkpoint.Slicer()
            chkpt.set_selection(slicer[1:,1:,0:])

        outfile = open(filename.replace('.h5', '.vtk'), 'w')
        outfile.write("# vtk DataFile Version 3.1\n")
        outfile.write("m2: astrophysical MHD code\n")
        outfile.write("ASCII\n" if args.ascii else "BINARY\n")
        outfile.write("DATASET STRUCTURED_GRID\n")
        outfile.write("DIMENSIONS %d %d %d\n" % (n1+1, n2+1, n3+1))
        outfile.write("\nPOINTS %d FLOAT\n" % ((n1+1)*(n2+1)*(n3+1)))

        # assumes machine is little-endian (VTK uses big-endian)
        if args.ascii:
            for P in zip(p[...,0].flat,
                         p[...,1].flat,
                         p[...,2].flat):
                outfile.write("%+10.8e %+10.8e %+10.8e\n" % P)
        else:
            p.byteswap().tofile(outfile)

        outfile.write("\nCELL_DATA %d\n" % (n1 * n2 * n3))

        for field in 'pd':
            outfile.write("\nSCALARS %s FLOAT\n" % field)
            outfile.write("LOOKUP_TABLE default\n")
            f = chkpt.get_field(field).astype(np.float32)
            f = f.swapaxes(0, 2)
            if args.ascii:
                for z in f.flat:
                    outfile.write("%+10.8e\n" % z)
            else:
                f.byteswap().tofile(outfile)

        for field in 'Bv':
            outfile.write("\nVECTORS %s FLOAT\n" % field)
            f = np.zeros([n1, n2, n3, 3], dtype=np.float32)
            f[...,0], f[...,1], f[...,2] = chkpt.get_vector_field(field)
            f = f.swapaxes(0, 2)
            if args.ascii:
                for F in zip(f[...,0].flat,
                             f[...,1].flat,
                             f[...,2].flat):
                    outfile.write("%+10.8e %+10.8e %+10.8e\n" % F)
            else:
                f.byteswap().tofile(outfile)

        outfile.close()
        chkpt.close()



class VtkDataExporterCommand(command.Command):
    _alias = 'tovtk'

    def configure_parser(self, parser):
        parser.add_argument("filenames", nargs='+')
        parser.add_argument("--ascii", action='store_true')

    def run(self, args):
        exporter = VtkStructuredGridExporter()
        exporter.run(args) 


