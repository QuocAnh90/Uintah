#
#  The contents of this file are subject to the University of Utah Public
#  License (the "License"); you may not use this file except in compliance
#  with the License.
#  
#  Software distributed under the License is distributed on an "AS IS"
#  basis, WITHOUT WARRANTY OF ANY KIND, either express or implied. See the
#  License for the specific language governing rights and limitations under
#  the License.
#  
#  The Original Source Code is SCIRun, released March 12, 2001.
#  
#  The Original Source Code was developed by the University of Utah.
#  Portions created by UNIVERSITY are Copyright (C) 2001, 1994 
#  University of Utah. All Rights Reserved.
#

# example of Writer
# by Samsonov Alexei
# October 2000 
# NOTE: if file to be splitted, uncomment corresponding lines in the file

catch {rename SCIRun_DataIO_MatrixWriter ""}

itcl_class SCIRun_DataIO_MatrixWriter {
    inherit Module
    constructor {config} {
	set name MatrixWriter
	set_defaults
    }
    method set_defaults {} {
	global $this-filetype $this-confirm $this-split env
	set $this-filetype Binary
	set $this-split 0
      	set $this-confirm 1
	puts "conf is: $env(SCI_CONFIRM_OVERWRITE)"
	if { [info exists env(SCI_CONFIRM_OVERWRITE)] && 
	     ([string equal 0 $env(SCI_CONFIRM_OVERWRITE)] ||
	      [string equal -nocase no $env(SCI_CONFIRM_OVERWRITE)]) } {
	    set $this-confirm 0
	}

    }
    method overwrite {} {
	global $this-confirm $this-filetype
	if {[info exists $this-confirm] && [info exists $this-filename] && \
		[set $this-confirm] && [file exists [set $this-filename]] } {
	    set value [tk_messageBox -type yesno -parent . \
			   -icon warning -message \
			   "File [set $this-filename] already exists.\n Would you like to overwrite it?"]
	    if [string equal "no" $value] { return 0 }
	}
	return 1
    }

    method ui {} {
	global env
	set w .ui[modname]
	if {[winfo exists $w]} {
	    set child [lindex [winfo children $w] 0]
	    # $w withdrawn by $child's procedures
	    raise $child
	    return;
	}
	
	toplevel $w
	set initdir ""

	# place to put preferred data directory
	# it's used if $this-filename is empty
	if {[info exists env(PSE_DATA)]} {
	    set initdir $env(PSE_DATA)
	}
	
	if { $initdir==""} {
	    if {[info exists env(SCI_DATA)]} {
		set initdir $env(SCI_DATA)
	    } elseif {[info exists env(SCIRUN_DATA)]} {
		set initdir $env(SCIRUN_DATA)
	    }
	}

	#######################################################
	# to be modified for particular reader

	# extansion to append if no extension supplied by user
	set defext ".mat"
	
	# name to appear initially
	set defname "MyMatrix"
	set title "Save matrix file"

	# file types to appers in filter box
	set types {
	    {{Matrix}        {.mat} }
	    {{All Files}       {.*} }
	}
	
	######################################################
	
	makeSaveFilebox \
		-parent $w \
		-filevar $this-filename \
		-command "$this-c needexecute; destroy $w" \
		-cancel "destroy $w" \
		-title $title \
		-filetypes $types \
	        -initialfile $defname \
		-initialdir $initdir \
		-defaultextension $defext \
		-formatvar $this-filetype \
		-splitvar $this-split \
	        -confirmvar $this-confirm
    }
}
