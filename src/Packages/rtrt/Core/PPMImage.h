
#ifndef PPMIMAGE_H
#define PPMIMAGE_H 1

#include <Packages/rtrt/Core/Color.h>
#include <Packages/rtrt/Core/Array2.h>

#include <vector>
#include <string>
#include <iostream>
#include <fstream>

#include <stdio.h>
#include <stdlib.h>

namespace rtrt {
class PPMImage;
}

namespace SCIRun {
void Pio(Piostream&, rtrt::PPMImage&);
}

namespace rtrt {

using std::ifstream;
using std::ofstream;
using std::cout;
using std::cerr;
using std::endl;
using std::vector;
using std::string;

class PPMImage
{

 protected:
  unsigned            u_,v_;
  unsigned            max_;
  bool                valid_;
  vector<rtrt::Color> image_;
  bool                flipped_;

  void eat_comments_and_whitespace(ifstream &str);

 public:
  PPMImage() {} // for Pio.
  PPMImage(const string& s, bool flip=false) 
    : flipped_(flip) 
  { 
    valid_ = read_image(s.c_str());
  }
  PPMImage(int nu, int nv, bool flip=false) 
    : u_(nu), v_(nv), valid_(false), flipped_(flip) 
  {
    image_.resize(u_*v_);
  }

  virtual ~PPMImage() {}

  friend void SCIRun::Pio(SCIRun::Piostream&, PPMImage&);

  unsigned get_width() { return u_; }
  unsigned get_height() { return v_; }
  unsigned get_size() { return max_; }

  void get_dimensions_and_data(Array2<rtrt::Color> &c, int &nu, int &nv) {
    if (valid_) {
      c.resize(u_+2,v_+2);  // image size + slop for interpolation
      nu=u_;
      nv=v_;
      for (unsigned v=0; v<v_; ++v)
        for (unsigned u=0; u<u_; ++u)
          c(u,v)=image_[v*u_+u];
    } else {
      c.resize(0,0);
      nu=0;
      nv=0;
    }
  }

  bool valid() { return valid_; }

  rtrt::Color &operator()(unsigned u, unsigned v)
  {
    if (v>=v_) v=v_-1;
    if (u>=u_) u=u_-1;
    return image_[v*u_+u];
  }

  const rtrt::Color &operator()(unsigned u, unsigned v) const
  {
    if (v>=v_) v=v_-1;
    if (u>=u_) u=u_-1;
    return image_[v*u_+u];
  }

  bool write_image(const char* filename, int bin=1)
  {
    ofstream outdata(filename);
    if (!outdata.is_open()) {
      cerr << "PPMImage: ERROR: I/O fault: couldn't write image file: "
	   << filename << "\n";
      return false;
    }
    if (bin)
      outdata << "P6\n# PPM binary image created with rtrt\n";
    else
      outdata << "P3\n# PPM ASCII image created with rtrt\n";

    outdata << u_ << " " << v_ << "\n";
    outdata << "255\n";

    unsigned char c[3];
    if (bin) {
      for(unsigned v=0;v<v_;++v){
	for(unsigned u=0;u<u_;++u){
	  c[0]=(unsigned char)(image_[v*u_+u].red()*255);
	  c[1]=(unsigned char)(image_[v*u_+u].green()*255);
	  c[2]=(unsigned char)(image_[v*u_+u].blue()*255);
	  outdata.write((char *)c, 3);
	}
      }
    } else {
      int count=0;
      for(unsigned v=0;v<v_;++v){
	for(unsigned u=0;u<u_;++u, ++count){
	  if (count == 5) { outdata << "\n"; count=0; }
	  outdata << (int)(image_[v*u_+u].red()*255) << " ";
	  outdata << (int)(image_[v*u_+u].green()*255) << " ";
	  outdata << (int)(image_[v*u_+u].blue()*255) << " ";
	}
      }
    }
    return true;
  }

  bool read_image(const char* filename)
  {
    ifstream indata(filename);
    unsigned char color[3];
    string token;
    
    if (!indata.is_open()) {
      cerr << "PPMImage: ERROR: I/O fault: no such file: " 
           << filename << "\n";
      valid_ = false;
      return false;
    }
    
    indata >> token; // P6
    if (token != "P6" && token != "P3") {
      cerr << "PPMImage: WARNING: format error: file not a PPM: "
           << filename << "\n";
    }

    cerr << "PPMImage: reading image: " << filename;
    if (flipped_)
      cerr << " (flipped!)";
    cerr << endl;

    eat_comments_and_whitespace(indata);
    indata >> u_ >> v_;
    eat_comments_and_whitespace(indata);
    indata >> max_;
    eat_comments_and_whitespace(indata);
    image_.resize(u_*v_);
    if (token == "P6") {
      for(unsigned v=0;v<v_;++v){
	for(unsigned u=0;u<u_;++u){
	  indata.read((char*)color, 3);
          if (flipped_) {
            image_[(v_-v-1)*u_+u]=rtrt::Color(color[0]/(double)max_,
                                              color[1]/(double)max_,
                                              color[2]/(double)max_);
          } else {
            image_[v*u_+u]=rtrt::Color(color[0]/(double)max_,
                                       color[1]/(double)max_,
                                       color[2]/(double)max_);
          }
	}
      }    
    } else { // P3
      int r, g, b;
      for(unsigned v=0;v<v_;++v){
	for(unsigned u=0;u<u_;++u){
	  indata >> r >> g >> b;
          if (flipped_) {
            image_[(v_-v-1)*u_+u]=rtrt::Color(r/(double)max_,
                                              g/(double)max_,
                                              b/(double)max_);
          } else {
            image_[v*u_+u]=rtrt::Color(r/(double)max_,
                                       g/(double)max_,
                                       b/(double)max_);
          }
	}
      }    
    }
    valid_ = true;
    return true;
  }
};

} // end namespace

#endif
