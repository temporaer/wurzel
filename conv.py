import vtk
import numpy as np
import numpy_support as ns

def dat2mhd(fn):
		with open("L2_17aug.dat") as fd:
				D = np.fromfile(file=fd, dtype=np.uint8).reshape((256,256,120) ).astype("float32")/255.0

		D = np.log(D+1)

		from scipy.ndimage.interpolation import zoom
		D = zoom(D, [1,1,256.0/120.0])

		flat_d      = D.transpose(2,1,0).flatten()
		vtk_d_array = ns.numpy_to_vtk(flat_d)

		image = vtk.vtkImageData() 

		points= image.GetPointData()
		points.SetScalars(vtk_d_array)

		image.SetDimensions(D.shape)

		image.Update()

		w = vtk.vtkMetaImageWriter()
		w.SetFileName("bla.hdr")
		w.SetInput(image)
		w.Write()


def mhd2npy(fn):
		# back the other way:
		r     = vtk.vtkMetaImageReader()
		r.SetFileName(fn)
		r.Update()
		image = r.GetOutput()
		vtk_data = image.GetPointData().GetScalars()
		numpy_data = ns.vtk_to_numpy(vtk_data)
		dims = image.GetDimensions()
		numpy_data = numpy_data.reshape(dims[2], dims[1], dims[0])
		numpy_data = numpy_data.transpose(2,1,0)
		return numpy_data
		# numpy_data[x,y,z]...

if __name__ == "__main__":
		#dat2mhd("L2_17aug.dat")
		mhd2npy("/home/local/l_schulz/tmp/build-itk/InsightApplications-3.20.0/build/Overlay.mhd")
