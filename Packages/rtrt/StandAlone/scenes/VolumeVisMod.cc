#include <Packages/rtrt/Core/TimeObj.h>
#include <Packages/rtrt/Core/Array1.h>
#include <Packages/rtrt/Core/BrickArray3.h>
#include <Packages/rtrt/Core/Camera.h>
#include <Packages/rtrt/Core/CatmullRomSpline.h>
#include <Packages/rtrt/Core/Group.h>
#include <Packages/rtrt/Core/Light.h>
#include <Packages/rtrt/Core/Phong.h>
#include <Packages/rtrt/Core/LambertianMaterial.h>
#include <Packages/rtrt/Core/Scene.h>
#include <Packages/rtrt/Core/VolumeVis.h>
#include <Packages/rtrt/Core/VolumeVisDpy.h>
#include <Packages/rtrt/Core/CutPlane.h>
#include <Packages/rtrt/Core/PlaneDpy.h>
#include <Core/Thread/Thread.h>
#include <nrrd.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <stdlib.h>
#include <vector>
#include <string>

using namespace std;
using namespace rtrt;
using SCIRun::Thread;

void get_material(Array1<Color> &matls, Array1<AlphaPos> &alphas) {

  matls.add(Color(0,0,1));
  matls.add(Color(0,0.4,1));
  matls.add(Color(0,0.8,1));
  matls.add(Color(0,1,0.8));
  matls.add(Color(0,1,0.4));
  matls.add(Color(0,1,0));
  matls.add(Color(0.4,1,0));
  matls.add(Color(0.8,1,0));
  matls.add(Color(1,0.9176,0));
  matls.add(Color(1,0.8,0));
  matls.add(Color(1,0.4,0));
  matls.add(Color(1,0,0));

  alphas.add(AlphaPos(0       , 0));  // pos 0
  alphas.add(AlphaPos(28.0/255, 0));  // pos 28
  alphas.add(AlphaPos(64.0/255, 0.1));// pos 64
  alphas.add(AlphaPos(100.0/255,0));  // pos 100
  alphas.add(AlphaPos(156.0/255,0));  // pos 156
  alphas.add(AlphaPos(192.0/255,1));  // pos 192
  alphas.add(AlphaPos(228.0/255,0));  // pos 192 
  alphas.add(AlphaPos(1,        0));  // pos 192
}

void get_material2(Array1<Color> &matls, Array1<AlphaPos> &alphas) {

  float div = 1.0/255;
  matls.add(Color(255, 255, 255) * div);
  matls.add(Color(255, 255, 180) * div);
  matls.add(Color(255, 247, 120) * div);   
  matls.add(Color(255, 228, 80) * div);
  matls.add(Color(255, 204, 55) * div);   
  matls.add(Color(255, 163, 20) * div);
  matls.add(Color(255, 120, 0) * div);   
  matls.add(Color(230, 71, 0) * div);
  matls.add(Color(200, 41, 0) * div);   
  matls.add(Color(153, 18, 0) * div);
  matls.add(Color(102, 2, 0) * div);   
  matls.add(Color(52, 0, 0) * div);
  matls.add(Color(0, 0, 0) * div);


  alphas.add(AlphaPos(0       , 0));  // pos 0
  alphas.add(AlphaPos(28.0/255, 0));  // pos 28
  alphas.add(AlphaPos(64.0/255, 0.1));// pos 64
  alphas.add(AlphaPos(100.0/255,0));  // pos 100
  alphas.add(AlphaPos(156.0/255,0));  // pos 156
  alphas.add(AlphaPos(192.0/255,1));  // pos 192
  alphas.add(AlphaPos(228.0/255,0));  // pos 192 
  alphas.add(AlphaPos(1,        0));  // pos 192
}

VolumeVis *create_volume_from_nrrd(char *filename,
				   bool override_data_min, double data_min_in,
				   bool override_data_max, double data_max_in,
				   double spec_coeff, double ambient,
				   double diffuse, double specular,
				   VolumeVisDpy *dpy)
{
  BrickArray3<float> data;
  float data_min = FLT_MAX;
  float data_max = -FLT_MAX;
  Point minP, maxP;
  // Do the nrrd stuff
  Nrrd *n = nrrdNew();
  // load the nrrd in
  if (nrrdLoad(n,filename)) {
    char *err = biffGet(NRRD);
    cerr << "Error reading nrrd "<< filename <<": "<<err<<"\n";
    free(err);
    biffDone(NRRD);
    return 0;
  }
  // check to make sure the dimensions are good
  if (n->dim != 3) {
    cerr << "VolumeVisMod error: nrrd->dim="<<n->dim<<"\n";
    cerr << "  Can only deal with 3-dimensional scalar fields... sorry.\n";
    return 0;
  }
  // convert the type to floats if you need to
  nrrdBigInt num_elements = nrrdElementNumber(n);
  cerr << "Number of data members = " << num_elements << endl;
  if (n->type != nrrdTypeFloat) {
    cerr << "Converting type from ";
    switch(n->type) {
    case nrrdTypeUnknown: cerr << "nrrdTypeUnknown"; break;
    case nrrdTypeChar: cerr << "nrrdTypeChar"; break;
    case nrrdTypeUChar: cerr << "nrrdTypeUChar"; break;
    case nrrdTypeShort: cerr << "nrrdTypeShort"; break;
    case nrrdTypeUShort: cerr << "nrrdTypeUShort"; break;
    case nrrdTypeInt: cerr << "nrrdTypeInt"; break;
    case nrrdTypeUInt: cerr << "nrrdTypeUInt"; break;
    case nrrdTypeLLong: cerr << "nrrdTypeLLong"; break;
    case nrrdTypeULLong: cerr << "nrrdTypeULLong"; break;
    case nrrdTypeDouble: cerr << "nrrdTypeDouble"; break;
    default: cerr << "Unknown!!";
    }
    cerr << " to nrrdTypeFloat\n";
    Nrrd *new_n = nrrdNew();
    nrrdConvert(new_n, n, nrrdTypeFloat);
    // since the data was copied blow away the memory for the old nrrd
    nrrdNuke(n);
    n = new_n;
    cerr << "Number of data members = " << num_elements << endl;
  }
  // get the dimensions
  int nx, ny, nz;
  nx = n->axis[0].size;
  ny = n->axis[1].size;
  nz = n->axis[2].size;
  cout << "dim = (" << nx << ", " << ny << ", " << nz << ")\n";
  cout << "total = " << nz * ny * nz << endl;
  cout << "spacing = " << n->axis[0].spacing << " x "<<n->axis[1].spacing<< " x "<<n->axis[2].spacing<< endl;
  data.resize(nx,ny,nz); // resize the bricked data
  // get the physical bounds
  minP = Point(0,0,0);
  maxP = Point((nx - 1) * n->axis[0].spacing,
	       (ny - 1) * n->axis[1].spacing,
	       (nz - 1) * n->axis[2].spacing);
  // lets normalize the dimensions to 1
  Vector size = maxP - minP;
  // find the biggest dimension
  double max_dim = Max(Max(size.x(),size.y()),size.z());
  maxP = ((maxP-minP)/max_dim).asPoint();
  minP = Point(0,0,0);
  // copy the data into the brickArray
  cerr << "Number of data members = " << num_elements << endl;
  float *p = (float*)n->data; // get the pointer to the raw data
  for (int z = 0; z < nz; z++)
    for (int y = 0; y < ny; y++)
      for (int x = 0; x < nx; x++) {
	float val = *p++;
	data(x,y,z) = val;
	// also find the min and max
	if (val < data_min)
	  data_min = val;
	else if (val > data_max)
	  data_max = val;
      }
#if 0
  // compute the min and max of the data
  double dmin,dmax;
  nrrdMinMaxFind(&dmin,&dmax,n);
  data_min = (float)dmin;
  data_max = (float)dmax;
#endif
  // delete the memory that is no longer in use
  nrrdNuke(n);


  // override the min and max if it was passed in
  if (override_data_min)
    data_min = data_min_in;
  if (override_data_max)
    data_max = data_max_in;

  cout << "minP = "<<minP<<", maxP = "<<maxP<<endl;

  return new VolumeVis(data, data_min, data_max,
		       nx, ny, nz,
		       minP, maxP,
		       spec_coeff, ambient, diffuse,
		       specular, dpy);
}  

VolumeVis *create_volume_default(int scene_type, double val,
				 int nx, int ny, int nz,
				 bool override_data_min, double data_min_in,
				 bool override_data_max, double data_max_in,
				 double spec_coeff, double ambient,
				 double diffuse, double specular,
				 VolumeVisDpy *dpy)
{
  BrickArray3<float> data;

  // make sure dimensions are good
  if (!nx) nx = 1;
  if (!ny) ny = 1;
  if (!nz) nz = 1;
  // resize the bricked data
  data.resize(nx,ny,nz);
  // fill the data with values
  for (int x = 0; x < nx; x++) {
    for (int y = 0; y < ny; y++) {
      for (int z = 0; z < nz; z++) {
	switch (scene_type) {
	case 0:
	  data(x,y,z) = z;
	  break;
	case 1:
	  data(x,y,z) = y*z;
	  break;
	case 2:
	  data(x,y,z) = x*y*z;
	  break;
	case 3:
	  data(x,y,z) = y+z;
	  break;
	case 4:
	  data(x,y,z) = x+y+z;
	  break;
	case 5:
	  data(x,y,z) = val;
	  break;
	}
      }
    }
  }
  // compute the min and max of the data
  float data_min = 0;
  float data_max;
  switch (scene_type) {
  case 0:
    data_max = (nz-1);
    break;
  case 1:
    data_max = (ny-1) * (nz-1);
    break;
  case 2:
    data_max = (nx-1) * (ny-1) * (nz-1);
    break;
  case 3:
    data_max = (ny-1) + (nz-1);
    break;
  case 4:
    data_max = (nx-1) + (ny-1) + (nz-1);
    break;
  case 5:
    data_min = val;
    data_max = val;
    break;
  }
  // set the physical dimensions of the data
  Point minP(0,0,0), maxP(1,1,1);
  cout << "dim = (" << nx << ", " << ny << ", " << nz << ")\n";
  cout << "total = " << nz * ny * nz << endl;

  cout << "minP = "<<minP<<", maxP = "<<maxP<<endl;

  // override the min and max if it was passed in
  if (override_data_min)
    data_min = data_min_in;
  if (override_data_max)
    data_max = data_max_in;

  return new VolumeVis(data, data_min, data_max,
		       nx, ny, nz,
		       minP, maxP,
		       spec_coeff, ambient, diffuse,
		       specular, dpy);
}

#define RAINBOW_COLOR_MAP 0
#define INVERSE_BLACK_BODY_COLOR_MAP 1

extern "C" 
Scene* make_scene(int argc, char* argv[], int /*nworkers*/)
{
  int nx = 20;
  int ny = 30;
  int nz = 40;
  int scene_type = 0;
  vector<string> data_files;
  bool cut=false;
  bool do_phong = true;
  int ncolors=256;
  float t_inc = 0.01;
  double spec_coeff = 64;
  double ambient = 0.5;
  double diffuse = 1.0;
  double specular = 1.0;
  float val=1;
  bool override_data_min = false;
  bool override_data_max = false;
  float data_min_in = -1;
  float data_max_in = -1;
  float frame_rate = 3;
  int color_map_type = RAINBOW_COLOR_MAP;
  Color bgcolor(1.,1.,1.);
  
  for(int i=1;i<argc;i++){
    if(strcmp(argv[i], "-dim")==0){
      i++;
      if(sscanf(argv[i], "%dx%dx%d", &nx, &ny, &nz) != 3){
	cerr << "Error parsing dimensions: " << argv[i] << '\n';
	cerr << "dim = (" << nx << ", " << ny << ", " << nz << ")\n";
	exit(1);
      } 
    }
    else if (strcmp(argv[i], "-type")==0) {
      i++;
      if (strcmp(argv[i], "nrrd") == 0) {
	scene_type = 6;
	i++;
	data_files.push_back(argv[i]);
      } else if (strcmp(argv[i], "nrrdlist") == 0) {
	scene_type = 6;
	// open up the file and then suck in all the files
	i++;
	cout << "Reading nrrd file list from " << argv[i] << endl;
	ifstream in(argv[i]);
	while (in) {
	  string file;
	  in >> file;
	  data_files.push_back(file);
	  cout << "Nrrd file: "<<file<<"\n";
	}
	cout << "Read "<<data_files.size()<<" nrrd files.\n";
      } else {
	scene_type = atoi(argv[i]);
	if (scene_type == 5) {
	  i++;
	  val = atof(argv[i]);
	}
      }
    } else if(strcmp(argv[i], "-cut")==0){
      cut=true;
    } else if(strcmp(argv[i], "-lam")==0){
      do_phong = false;
    } else if(strcmp(argv[i], "-ncolors")==0){
      i++;
      ncolors = atoi(argv[i]);
    } else if(strcmp(argv[i], "-tinc")==0){
      i++;
      t_inc = atof(argv[i]);
    } else if(strcmp(argv[i], "-min")==0){
      i++;
      data_min_in = atof(argv[i]);
      override_data_min = true;
    } else if(strcmp(argv[i], "-max")==0){
      i++;
      data_max_in = atof(argv[i]);
      override_data_max = true;
    } else if(strcmp(argv[i], "-spow")==0){
      i++;
      spec_coeff = atof(argv[i]);
    } else if(strcmp(argv[i], "-ambient")==0){
      i++;
      ambient = atof(argv[i]);
    } else if(strcmp(argv[i], "-diffuse")==0){
      i++;
      diffuse = atof(argv[i]);
    } else if(strcmp(argv[i], "-specular")==0){
      i++;
      specular = atof(argv[i]);
    } else if(strcmp(argv[i], "-rate")==0){
      i++;
      frame_rate = atof(argv[i]);
    } else if(strcmp(argv[i], "-colormap")==0){
      i++;
      if (strcmp(argv[i], "rainbow")==0){
	color_map_type = RAINBOW_COLOR_MAP;
      } else if (strcmp(argv[i], "inversebb")==0){
	color_map_type = RAINBOW_COLOR_MAP;
      } else {
	cerr << "Unknown color map type.  Using rainbow.\n";
	color_map_type = RAINBOW_COLOR_MAP;
      }
    } else if(strcmp(argv[i], "-bgcolor")==0){
      float r,g,b;
      r = atof(argv[++i]);
      g = atof(argv[++i]);
      b = atof(argv[++i]);
      bgcolor = Color(r,g,b);
    } else {
      cerr << "Unknown option: " << argv[i] << '\n';
      cerr << "Valid options for scene: " << argv[0] << '\n';
      cerr << " -dim [int]x[int]x[int]";
      cerr << " -type [int or \"nrrd\"]\n";
      cerr << "\t\t0 - data = z\n";
      cerr << "\t\t1 - data = y*z\n";
      cerr << "\t\t2 - data = x*y*z\n";
      cerr << "\t\t3 - data = y+z\n";
      cerr << "\t\t4 - data = x+y+z\n";
      cerr << "\t\t5 [val=float] - data = val\n";
      cerr << "\t\tnrrd [path to nrrd file] (-dim parameter ignored)\n";
      cerr << " -cut - turn on the cutting plane\n";
      cerr << " -lam - use a lambertian surface\n";
      cerr << " -ncolors [int] - the size of the transfer function\n";
      cerr << " -tinc [float] - the number of samples per unit\n";
      cerr << " -spow [float] - the spectral exponent\n";
      cerr << " -ambient [float] - the ambient factor\n";
      cerr << " -diffuse [float] - the diffuse factor\n";
      cerr << " -specular [float] - the specular factor\n";
      cerr << " -rate [float] - frame rate of the time steps\n";
      cerr << " -colormap [string] - \"rainbow\" or \"inverse\"\n";
      cerr << " -bgcolor [float] [float] [float] - the three floats are r, g, b\n";
      return 0;
    }
  }

  // check the parameters
  if (scene_type == 6)
    if (data_files.size() == 0)
      // the file was not set
      scene_type = 0;
  
  Camera cam(Point(0.5,0.5,3), Point(0.5,0.5,0.5),
	     Vector(0,1,0), 40.0);
  
  double ambient_scale=1.0;
  
  //  double bgscale=0.5;
  //  Color bgcolor(bgscale*108/255., bgscale*166/255., bgscale*205/255.);
  
  Array1<Color> matls;
  Array1<AlphaPos> alphas;
  switch (color_map_type) {
  case INVERSE_BLACK_BODY_COLOR_MAP:
    get_material2(matls,alphas);
    break;
  case RAINBOW_COLOR_MAP:
  default:
    get_material(matls,alphas);
    break;
  }
  VolumeVisDpy *dpy = new VolumeVisDpy(matls, alphas, ncolors, t_inc);

  // Generate the data
  Object *obj;

  if (scene_type == 6) {
    if (data_files.size() > 1) {
      // create a timeobj
      TimeObj *timeobj = new TimeObj(frame_rate);
      for(unsigned int i = 0; i < data_files.size(); i++) {
	char *myfile = strdup(data_files[i].c_str());
	Object *volume = (Object*)create_volume_from_nrrd
	  (myfile, override_data_min, data_min_in,
	   override_data_max, data_max_in, spec_coeff,
	   ambient, diffuse, specular, dpy);
	// add the volume to the timeobj
	if (myfile)
	  free(myfile);
	if (volume == 0)
	  // there was a problem
	  continue;
	timeobj->add(volume);
      }
      obj = (Object*)timeobj;
    } else {
      char *myfile = strdup(data_files[0].c_str());
      obj = (Object*)create_volume_from_nrrd
	(myfile, override_data_min, data_min_in,
	 override_data_max, data_max_in, spec_coeff,
	 ambient, diffuse, specular, dpy);
      if (myfile)
	free(myfile);
    }
  } else {
    obj = (Object*) create_volume_default
	(scene_type, val, nx, ny, nz, override_data_min, data_min_in,
	 override_data_max, data_max_in, spec_coeff, ambient,
	 diffuse, specular, dpy);
  }

  (new Thread(dpy, "VolumeVis display thread"))->detach();

  if(cut){
    PlaneDpy* pd=new PlaneDpy(Vector(0,0,1), Point(0,0,0));
    obj=(Object*)new CutPlane(obj, pd);
    new Thread(pd, "Cutting plane display thread");
  }
  Group* all = new Group();
  all->add(obj);

  Plane groundplane ( Point(-500, 300, 0), Vector(7, -3, 2) );
  Color cup(0.9, 0.7, 0.3);
  Color cdown(0.0, 0.0, 0.2);

  Scene* scene=new Scene(all, cam,
			 bgcolor, cdown, cup, groundplane, 
			 ambient_scale);
  
  scene->add_light(new Light(Point(500,-300,300), Color(.8,.8,.8), 0));
  scene->shadow_mode=1;
  scene->attach_display(dpy);
  
  return scene;
}



