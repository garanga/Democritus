/*
 * vtk.h
 *
 *  Created on: Aug 24, 2017
 *      Author: pavel
 */

#ifndef VTK_H_
#define VTK_H_

#include "vtkAutoInit.h"
VTK_MODULE_INIT(vtkRenderingOpenGL2);
VTK_MODULE_INIT(vtkInteractionStyle);

#include "vtkPoints.h"

#include "vtkPolyData.h"
#include "vtkPolyDataMapper.h"

#include "vtkQuad.h"
#include "vtkCellArray.h"
#include "vtkUnstructuredGrid.h"

#include "vtkSphereSource.h"
#include "vtkGlyph3D.h"

#include "vtkDataSetMapper.h"
#include "vtkActor.h"
#include "vtkProperty.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"

#include "vtkXMLUnstructuredGridWriter.h"
#include "vtkXMLUnstructuredGridReader.h"

#include "vtkXMLPolyDataWriter.h"

#include "vtkAlgorithm.h"

#include "vtkDoubleArray.h"
#include "vtkPointData.h"

#endif /* VTK_H_ */
