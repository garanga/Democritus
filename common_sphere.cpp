/*
 * common_sphere.cpp
 *
 *  Created on: Aug 23, 2017
 *      Author: pavel
 */

#include "common.h"

#include "Sphere.h"               /* <===== replace this line if necessary */
#include "vtk.h"


#include <iostream>
#include <fstream>
#include <string>
#include <vector>


double lx, ly, x_0, y_0;
unsigned int no_of_particles;
double Time = 0.0;
std::ofstream fphase("phase.dat"), fenergy("energy.dat");

int nstep;
double timestep;
int nprint, nenergy;
std::ofstream flast("lastframe.dat");
Vector G;
std::vector<Sphere> particle;     /* <===== replace this line if necessary */


// These are used only here
void   init_system(char* fname);
double total_kinetic_energy();

int main(int argc, char* argv[])
{


    lx = 100.0; ly = 100.0;

//    std::cout << lx << "\t" << ly << std::endl;

    fenergy.precision(10);

//    init_system("../init_hopper/closed_hopper.random");
    init_system("../SandGlass/SandGlass.random");

    std::cout << G << std::endl;

    init_algorithm();
//
    phase_plot1(0);
//


//    std::cout << "\n"
//              << "\n"
//              << "Initial position" << "\n"
//              << "Particle 111:" << "\n"
//              << "predicted rtd0_ " << "\t" << particle[111].pos () << "\t"
//              << "predicted rtd2_ " << "\t" << particle[111].rtd2() << "\t"
//              << "predicted rtd3_ " << "\t" << particle[111].rtd3() << "\t"
//              << "predicted rtd4_ " << "\t" << particle[111].rtd4() << "\t"
//              << "\n"
//              << "\n";

    for(int i=0; i<nstep; i++)
    {

        std::cout << "Step No. " << i+1 << std::endl;

        std::cout << "energy: " << total_kinetic_energy() << std::endl;

        step();
        if((i+1)%nprint==0)
        {
            std::cout << "phase_plot: " << i+1 << "  " << particle.size() << " particles\n";
//            phase_plot (fphase);
//            std::cout << "./vtk/test" + std::to_string(i) + ".vtp" << std::endl;
            phase_plot1(i+1);
        }
//        if((i+1)%nenergy==0){
//            fenergy << Time << "\t" << total_kinetic_energy() << std::endl;
//        }
    }
//    phase_plot(flast);
}

void integrate()
{
    for(unsigned int i=0; i<particle.size(); i++)
    {

        if(particle[i].ptype()==0)
        {
//            std::cout << "\t" << "Integrate mobile particle No. " << i+1 << std::endl;
            particle[i].set_force_to_zero();

//            std::cout << "Particle No: " << i << "\t";

            particle[i].predict(timestep);


//            if (i == 320)
//            {
//                std::cout << "\n"
//                          << "\n"
//                          << "Predict info" << "\n"
//                          << "Particle 111:" << "\n"
//                          << "predicted rtd0_ " << "\t" << particle[i].pos () << "\t"
//                          << "predicted rtd2_ " << "\t" << particle[i].rtd2() << "\t"
//                          << "predicted rtd3_ " << "\t" << particle[i].rtd3() << "\t"
//                          << "predicted rtd4_ " << "\t" << particle[i].rtd4() << "\t"
//                          << "\n"
//                          << "\n";
//
//            }

        }
        else
        {
//            std::cout << "\t" << "Integrate not-mobile particle No. " << i+1 << std::endl;
            particle[i].boundary_conditions(i, timestep, Time);
        }
    }

//    std::cout << "Calculate forces" << std::endl;
    make_forces();
//    std::cout << "Correct" << std::endl;
    for(unsigned int i=0; i<particle.size(); i++)
    {
        if(particle[i].ptype()==0)
        {

//            std::cout << i << "\t" << "Correct!!!!!!!!!!!!!" << std::endl;
//            if (i == 111)
//            {
//                std::cout << "\n"
//                          << "\n"
//                          << "info before correct" << "\n"
//                          << "Particle 111:" << "\n"
//                          << "predicted rtd0_ " << "\t" << particle[i].pos () << "\t"
//                          << "predicted rtd2_ " << "\t" << particle[i].rtd2() << "\t"
//                          << "\n"
//                          << "\n";
//
//            }
//            if (i==320)
//            {
//                particle[i].correct(timestep,true);
//            }
//            else
//            {
                particle[i].correct(timestep);
//            }

//            if (i == 320)
//            {
//                std::cout << "\n"
//                          << "\n"
//                          << "Corrected info" << "\n"
//                          << "Particle 111:" << "\n"
//                          << "predicted rtd0_ " << "\t" << particle[i].pos () << "\t"
//                          << "predicted rtd2_ " << "\t" << particle[i].rtd2() << "\t"
//                          << "predicted rtd3_ " << "\t" << particle[i].rtd2() << "\t"
//                          << "predicted rtd4_ " << "\t" << particle[i].rtd2() << "\t"
//                          << "\n"
//                          << "\n";

//            }
        }
    }

//    std::cout << "Set periodic boundary conditions" << std::endl;

//    for(unsigned int i=0; i<particle.size(); i++)
//    {
//        particle[i].periodic_bc(x_0, y_0, lx, ly);
//    }

//    std::cout << "Increase time step" << std::endl;

    Time+=timestep;
}

void
init_system(char* fname)
{
    std::ifstream fin(fname);

    while(fin.peek() == '#')
    {
        std::string type;
        fin >> type;

        if     (type=="#gravity:")
        {
            fin >> G.x() >> G.y() >> G.phi();
            fin.ignore(100,'\n');
            std::cout << "gravity: " << G << std::endl;
        }
        else if(type=="#Time:")
        {
            fin >> Time;
            fin.ignore(100,'\n');
            std::cout << "Time: " << Time << std::endl;
        }
        else if(type=="#nstep:")
        {
            fin >> nstep;
            fin.ignore(100,'\n');
            std::cout << "nstep: " << nstep << std::endl;
        }
        else if(type=="#timestep:")
        {
            fin >> timestep;
            fin.ignore(100,'\n');
            std::cout << "timestep: " << timestep << std::endl;
        }
        else if(type=="#nprint:")
        {
            fin >> nprint;
            fin.ignore(100,'\n');
            std::cout << "nprint: " << nprint << std::endl;
        }
        else if(type=="#nenergy:")
        {
            fin >> nenergy;
            fin.ignore(100,'\n');
            std::cout << "nenergy: " << nenergy << std::endl;
        }
        else if(type=="#lx:")
        {
            fin >> lx;
            fin.ignore(100,'\n');
            std::cout << "lx: " << lx << std::endl;
        }
        else if(type=="#ly:")
        {
            fin >> ly;
            fin.ignore(100,'\n');
            std::cout << "ly: " << ly << std::endl;
        }
        else if(type=="#x_0:")
        {
            fin >> x_0;
            fin.ignore(100,'\n');
            std::cout << "x_0: " << x_0 << std::endl;
        }
        else if(type=="#y_0:")
        {
            fin >> y_0;
            fin.ignore(100,'\n');
            std::cout << "y_0: " << y_0 << std::endl;
        }
        else
        {
            std::cerr << "init: unknown global property: " << type << std::endl;
            abort();
        }
    }

    while(fin)
    {
        Sphere s;
        fin >> s;
        if(fin)
        {
            particle.push_back(s);
        }
    }
    no_of_particles = particle.size();
    std::cout << no_of_particles << " particles read\n" << std::flush;
}

double
total_kinetic_energy()
{
    double sum = 0.0;
    for(unsigned int i=0; i<particle.size(); i++)
    {
        if(particle[i].ptype()==0)
        {
            sum+=particle[i].kinetic_energy();
        }
    }
    return sum;
}

void phase_plot(std::ostream& os)
{
    os << "#NewFrame\n";
    os << "#no_of_particles: " << no_of_particles << std::endl;
    os << "#compressed: no\n";
    os << "#type: SphereXYPhiVxVyOmegaRMFixed25\n";
    os << "#gravity: " << G.x() << " " << G.y() << " " << G.phi() << std::endl;
    os << "#Time: " << Time << std::endl;
    os << "#timestep: " << timestep << std::endl;
    os << "#EndOfHeader\n";
    for(unsigned int i=0; i<particle.size(); i++)
    {
        os << particle[i];
    }
    os << std::flush;
}

void phase_plot1(int i)
{
    vtkSmartPointer<vtkPoints> points =
        vtkSmartPointer<vtkPoints>::New();

    vtkSmartPointer<vtkDoubleArray> radii =
        vtkSmartPointer<vtkDoubleArray>::New();
    radii->SetNumberOfValues(particle.size());

    int cnt=0;
    for (auto& p : particle)
    {
        points -> InsertNextPoint(p.x(),p.y(),0.0);
        radii->SetValue(cnt, p.r());
        ++cnt;
    }

    vtkSmartPointer<vtkPolyData> polyData =
        vtkSmartPointer<vtkPolyData>::New();

    polyData->SetPoints(points);

    polyData->GetPointData()->SetScalars(radii);

    // Write file
    vtkSmartPointer<vtkXMLPolyDataWriter> writer =
            vtkSmartPointer<vtkXMLPolyDataWriter>::New();
//    writer->SetFileName("./vtk/test.vtp");

    std::string name = "./vtk/test-" + std::to_string(i) + ".vtp";

    writer->SetFileName(name.c_str());

//    writer->SetFileName("./vtk/test-" + std::to_string(i) + ".vtp");
    writer->SetInputData(polyData);
    writer->SetDataModeToAscii();
    writer->Write();
}

//    // Set color
//
//    points->InsertNextPoint( 0,0,0);
//    points->InsertNextPoint( 5,0,0);
//    points->InsertNextPoint(10,0,0);
//    vtkSmartPointer<vtkUnsignedCharArray> colors =
//        vtkSmartPointer<vtkUnsignedCharArray>::New();
//    colors->SetName("colors");
//    colors->SetNumberOfComponents(3);
//    unsigned char r[3] = {255,0,0};
//    unsigned char g[3] = {0,255,0};
//    unsigned char b[3] = {0,0,255};
//    colors->InsertNextTypedTuple(r);
//    colors->InsertNextTypedTuple(g);
//    colors->InsertNextTypedTuple(b);
//    polyData->GetPointData()->SetScalars(colors);
//    glyph3D->SetColorModeToColorByScalar();

//    // Draw spheres
//
//    vtkSmartPointer<vtkSphereSource> sphereSource =
//        vtkSmartPointer<vtkSphereSource> :: New();//
//    sphereSource->SetRadius(1.0);
//
//    vtkSmartPointer<vtkGlyph3D> glyph3D =
//        vtkSmartPointer<vtkGlyph3D> :: New();//
//    glyph3D->SetScaleModeToScaleByScalar();
//    glyph3D->SetScaleFactor(1.0)//
//    glyph3D->SetSourceConnection(sphereSource->GetOutputPort());//
//    glyph3D->SetInputData(polyData);//
//    glyph3D->ScalingOff();
//    glyph3D->Update();
//
//    // Create a mapper and actor
//
//    vtkSmartPointer<vtkPolyDataMapper> mapper =
//        vtkSmartPointer<vtkPolyDataMapper>::New();
//    mapper->SetInputConnection(glyph3D->GetOutputPort());
//
//    vtkSmartPointer<vtkActor> actor =
//        vtkSmartPointer<vtkActor>::New();
//    actor->SetMapper(mapper);
//
//    // Visualize
//
//    vtkSmartPointer<vtkRenderer> renderer =
//       vtkSmartPointer<vtkRenderer>::New();
//    renderer->AddActor(actor);
//    renderer->SetBackground(1,1,1); // Background color white
//
//    vtkSmartPointer<vtkRenderWindow> renderWindow =
//       vtkSmartPointer<vtkRenderWindow>::New();
//    renderWindow->AddRenderer(renderer);
//
//    vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor =
//       vtkSmartPointer<vtkRenderWindowInteractor>::New();
//    renderWindowInteractor->SetRenderWindow(renderWindow);
//
//    renderWindow->Render();
//    renderWindowInteractor->Start();
