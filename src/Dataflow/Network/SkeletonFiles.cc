/* SkeletonFiles.cc */

namespace PSECore {
namespace Dataflow {

char component_skeleton[] = \
"/*\n"
" *  %s.cc:\n" /* component name */
" *\n"
" *  Written by:\n"
" *   %s\n" /* author name */
" *   %s\n"   /* today's date */
" *\n"
" */\n"
"\n"
"#include <PSECore/Dataflow/Module.h>\n"
"#include <SCICore/Malloc/Allocator.h>\n"
"\n"
"#include <%s/share/share.h>\n" /* package name */
"\n"
"namespace %s {\n" /* package name */
"namespace Modules {\n"
"\n"
"using namespace PSECore::Dataflow;\n"
"\n"
"class %sSHARE %s : public Module {\n" /* package name, component name */
"" 
"public:\n"
"  %s(const clString& id);\n" /* component name */
"\n"
"  virtual ~%s();\n" /* component name */
"\n"
"  virtual void execute();\n"
"\n"
"  virtual void tcl_command(TCLArgs&, void*);\n"
"};\n"
"\n"
"extern \"C\" %sSHARE Module* make_%s(const clString& id) {\n" /* package name, component name */
"  return new %s(id);\n" /* component name */
"}\n"
"\n"
"%s::%s(const clString& id)\n" /* component name, component name */
"  : Module(\"%s\", id, Source)\n" /* component name */
"{\n"
"}\n"
"\n"
"%s::~%s()" /* component name, component name */
"{\n"
"}\n"
"\n"
"void %s::execute()" /* component name */
"{\n"
"}\n"
"\n"
"void %s::tcl_command(TCLArgs& args, void* userdata)\n" /* component name */
"{\n"
"  Module::tcl_command(args, userdata);\n"
"}\n"
"\n"
"} // End namespace Modules\n"
"} // End namespace %s\n" /* package name */
"\n"
"\n";

char gui_skeleton[] = \
"itcl_class %s_%s_%s {\n" /* package name, category name, component name */
"    inherit Module\n"
"    constructor {config} {\n"
"        set name %s\n" /* component name */
"        set_defaults\n"
"    }\n"
"\n"
"    method set_defaults {} {\n"
"    }\n"
"\n"
"    method ui {} {\n"
"        set w .ui[modname]\n"
"        if {[winfo exists $w]} {\n"
"            raise $w\n"
"            return\n"
"        }\n"
"        toplevel $w\n"
"    }\n"
"}\n"
"\n"
"\n";

char dllentry_skeleton[] = \
"/* DllEntry.cc */\n"
"\n"
"#ifdef _WIN32\n"
"\n"
"#include <afxwin.h>\n"
"#include <stdio.h>\n"
"\n"
"BOOL APIENTRY DllMain(HANDLE hModule, \n"
"                      DWORD  ul_reason_for_call, \n"
"                      LPVOID lpReserved)\n"
"{\n"
"#ifdef DEBUG\n"
"  char reason[100]=\"\\0\";\n"
"  printf(\"\\n*** %%sd.dll is initializing {%%s,%%d} ***\\n\",__FILE__,__LINE__);\n" /* package name */
"  printf(\"*** hModule = %%d ***\\n\",hModule);\n"
"  switch (ul_reason_for_call){\n"
"    case DLL_PROCESS_ATTACH:sprintf(reason,\"DLL_PROCESS_ATTACH\"); break;\n"
"    case DLL_THREAD_ATTACH:sprintf(reason,\"DLL_THREAD_ATTACH\"); break;\n"
"    case DLL_THREAD_DETACH:sprintf(reason,\"DLL_THREAD_DETACH\"); break;\n"
"    case DLL_PROCESS_DETACH:sprintf(reason,\"DLL_PROCESS_DETACH\"); break;\n"
"  }\n"
"  printf(\"*** ul_reason_for_call = %%s ***\\n\",reason);\n"
"#endif\n"
"  return TRUE;\n"
"}\n"
"\n"
"#endif\n"
"\n"
"\n";

char share_skeleton[] = \
"/* share.h */\n"
"\n"
"#undef %sSHARE\n" /* package name */
"\n"
"#ifdef _WIN32\n"
"  #if defined(BUILD_%s)\n" /* package name */
"    #define %sSHARE __declspec(dllexport)\n" /* package name */
"  #else\n"
"    #define %sSHARE __declspec(dllimport)\n" /* package name */
"  #endif\n" 
"#else\n" 
"  #define %sSHARE\n" /* package name */
"#endif\n"
"\n"
"\n";

char package_submk_skeleton[] = \
"include $(OBJTOP_ABS)/scripts/largeso_prologue.mk\n"
"\n"
"SRCDIR := %s\n" /* package name */
"\n"
"SUBDIRS := $(SRCDIR)/GUI $(SRCDIR)/Datatypes \\\n"
"        $(SRCDIR)/Modules\n"
"\n"
"include $(OBJTOP_ABS)/scripts/recurse.mk\n"
"\n"
"PSELIBS := PSECore SCICore\n"
"LIBS := $(TK_LIBRARY) $(GL_LIBS) -lm\n"
"\n"
"include $(OBJTOP_ABS)/scripts/largeso_epilogue.mk\n"
"\n"
"\n";

char modules_submk_skeleton[] = \
"# *** NOTE ***\n"
"#\n"
"# Do not remove or modify the comment line:\n"
"#\n"
"# #[INSERT NEW ?????? HERE]\n"
"#\n"
"# It is required by the Component Wizard to properly edit this file.\n"
"# if you want to edit this file by hand, see the \"Create A New Component\"\n"
"# documentation on how to do it correctly.\n"
"\n"
"SRCDIR := %s/Modules\n" /* package name */
"\n"
"SUBDIRS := \\\n"
"#[INSERT NEW CATEGORY DIR HERE]\n"
"\n"
"include $(OBJTOP_ABS)/scripts/recurse.mk\n"
"\n"
"\n";

char category_submk_skeleton[] = \
"# *** NOTE ***\n"
"#\n"
"# Do not remove or modify the comment line:\n"
"#\n"
"# #[INSERT NEW ?????? HERE]\n"
"#\n"
"# It is required by the Component Wizard to properly edit this file.\n"
"# if you want to edit this file by hand, see the \"Create A New Component\"\n"
"# documentation on how to do it correctly.\n"
"\n"
"include $(OBJTOP_ABS)/scripts/smallso_prologue.mk\n"
"\n"
"SRCDIR   := %s/Modules/%s\n" /* package name, category name */
"\n"
"SRCS     += \\\n"
"#[INSERT NEW CODE FILE HERE]\n"
"\n"
"PSELIBS := PSECore/Datatypes PSECore/Dataflow \\\n"
"        SCICore/Persistent SCICore/Containers SCICore/Util \\\n"
"        SCICore/Exceptions SCICore/Thread SCICore/TclInterface \\\n"
"        SCICore/Geom SCICore/Datatypes SCICore/Geometry \\\n"
"        SCICore/TkExtensions\n"
"LIBS := $(TK_LIBRARY) $(GL_LIBS) -lm\n"
"\n"
"include $(OBJTOP_ABS)/scripts/smallso_epilogue.mk\n"
"\n"
"\n";

char datatypes_submk_skeleton[] = \
"# *** NOTE ***\n"
"#\n"
"# Do not remove or modify the comment line:\n"
"#\n"
"# #[INSERT NEW ?????? HERE]\n"
"#\n"
"# It is required by the Component Wizard to properly edit this file.\n"
"# if you want to edit this file by hand, see the \"Create A New Component\"\n"
"# documentation on how to do it correctly.\n"
"\n"
"include $(OBJTOP_ABS)/scripts/smallso_prologue.mk\n"
"\n"
"SRCDIR   := %s/Datatypes\n" /* package name */
"\n"
"SRCS     += \\\n"
"#[INSERT NEW CODE FILE HERE]\n"
"\n"
"PSELIBS :=\n"
"LIBS :=\n"
"\n"
"include $(OBJTOP_ABS)/scripts/smallso_epilogue.mk"
"\n"
"\n";

char gui_submk_skeleton[] = \
"# *** NOTE ***\n"
"#\n"
"# Do not remove or modify the comment line:\n"
"#\n"
"# #[INSERT NEW ?????? HERE]\n"
"#\n"
"# It is required by the Component Wizard to properly edit this file.\n"
"# if you want to edit this file by hand, see the \"Create A New Component\"\n"
"# documentation on how to do it correctly.\n"
"\n"
"SRCDIR := %s/GUI\n" /* package name */
"\n"
"ALLTARGETS := $(ALLTARGETS) $(SRCDIR)/tclIndex\n"
"\n"
"$(SRCDIR)/tclIndex: \\\n"
"#[INSERT NEW TCL FILE HERE]\n"
"\n"
"CLEANPROGS := $(CLEANPROGS) $(SRCDIR)/tclIndex\n"
"\n"
"\n";

} // Dataflow
} // PSECore








