
#ifndef UINTAH_HOMEBREW_VarLabel_H
#define UINTAH_HOMEBREW_VarLabel_H

#include <string>

namespace Uintah {
  namespace Grid {
     class TypeDescription;
    
    /**************************************
      
      CLASS
        VarLabel
      
        Short Description...
      
      GENERAL INFORMATION
      
        VarLabel.h
      
        Steven G. Parker
        Department of Computer Science
        University of Utah
      
        Center for the Simulation of Accidental Fires and Explosions (C-SAFE)
      
        Copyright (C) 2000 SCI Group
      
      KEYWORDS
        VarLabel
      
      DESCRIPTION
        Long description...
      
      WARNING
      
      ****************************************/
    
    class VarLabel {
    public:
       VarLabel(const std::string&, const TypeDescription*);
    private:
       std::string d_name;
       const TypeDescription* d_td;

       VarLabel(const VarLabel&);
       VarLabel& operator=(const VarLabel&);
    };

    
  } // end namespace Grid
} // end namespace Uintah

//
// $Log$
// Revision 1.1  2000/04/19 05:26:15  sparker
// Implemented new problemSetup/initialization phases
// Simplified DataWarehouse interface (not finished yet)
// Made MPM get through problemSetup, but still not finished
//
//

#endif

