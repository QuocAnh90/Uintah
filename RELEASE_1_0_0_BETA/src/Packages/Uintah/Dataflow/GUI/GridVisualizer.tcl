#  GridVisualizer.tcl
#  Written by:
#   James Bigler
#   Department of Computer Science
#   University of Utah
#   Aug. 2000
#  Copyright (C) 2000 SCI Group

itcl_class Uintah_Visualization_GridVisualizer {
    inherit Module

    protected var_list ""
    protected mat_list ""
    protected type_list ""
    protected var_val_list {}
    protected time_list {}
    protected num_materials 0
    protected num_colors 240

    protected matrix_types {"Determinant" "Trace" "Norm"}
    protected vector_types {"length" "length2" "x" "y" "z"}
    protected num_m_type
    protected num_v_type
    # this represents the number of graphs made
    # when a new graph is created then this number is incremented
    # Thus repeatidly punching graph will continue to make new ones
    # without replacing old ones.
    protected graph_id 0

    constructor {config} {
	set name GridVisualizer
	set_defaults
    }
    
    method set_defaults {} {
	#the grid colors
	global $this-level1_grid_color
	set $this-level1_grid_color red
	global $this-level2_grid_color
	set $this-level2_grid_color green
	global $this-level3_grid_color
	set $this-level3_grid_color yellow
	global $this-level4_grid_color
	set $this-level4_grid_color magenta
	global $this-level5_grid_color
	set $this-level5_grid_color cyan
	global $this-level6_grid_color
	set $this-level6_grid_color blue
	
	#the node colors
	global $this-level1_node_color
	set $this-level1_node_color red
	global $this-level2_node_color
	set $this-level2_node_color green
	global $this-level3_node_color
	set $this-level3_node_color yellow
	global $this-level4_node_color
	set $this-level4_node_color magenta
	global $this-level5_node_color
	set $this-level5_node_color cyan
	global $this-level6_node_color
	set $this-level6_node_color blue
	global $this-nl
	set $this-nl 0

	#plane selection
	global $this-plane_on
	set $this-plane_on 0
	global $this-node_select_on
	set $this-node_select_on 0
	$this-c needexecute
	global $this-var_orientation
	set $this-var_orientation 0
	global $this-index_l
	set $this-index_l 0
	global $this-index_x
	set $this-index_x 0
	global $this-index_y
	set $this-index_y 0
	global $this-index_z
	set $this-index_z 0

	#sphere stuff
	global $this-radius
	set $this-radius 0.01
	global $this-polygons
	set $this-polygons 8

	#var graphing stuff
	global $this-curr_var
	set var_list ""

	# selection stuff
	set num_m_type [llength $matrix_types]
	set num_v_type [llength $vector_types]
    }
    method make_color_menu {w v n} {
	$w add radiobutton \
		-variable $v \
		-command $n \
		-label red \
		-value red
	$w add radiobutton \
		-variable $v \
		-command $n \
		-label green \
		-value green
	$w add radiobutton \
		-variable $v \
		-command $n \
		-label blue \
		-value blue
	$w add radiobutton \
		-variable $v \
		-command $n \
		-label yellow \
		-value yellow
	$w add radiobutton \
		-variable $v \
		-command $n \
		-label cyan \
		-value cyan
	$w add radiobutton \
		-variable $v \
		-command $n \
		-label magenta \
		-value magenta

    }
    # this is used for seting up the colors used for
    # the levels and grid points
    method setup_color {w text nl n} {
	for {set i 1} { $i <= $nl} {incr i} {
	    set st "$w.l$i"
	    append st "$text"
#	    puts $st
	    menubutton $st -text "Level [expr $i - 1] $text color" \
		    -menu $st.list -relief groove
	    pack $st -side top -anchor w -padx 2 -pady 2
	    
	    menu $st.list
	    set var "$this-level$i"
	    append var "_$text"
	    append var "_color"
	    make_color_menu $st.list $var $n
	}
    }
    method isVisible {} {
	if {[winfo exists .ui[modname]]} {
	    return 1
	} else {
	    return 0
	}
    }
    method Rebuild {} {
	set w .ui[modname]

	$this destroyFrames
	$this makeFrames $w
    }
    method build {} {
	set w .ui[modname]

	$this makeFrames $w
    }
    method destroyFrames {} {
	set w .ui[modname]
	destroy $w.colormenus
	destroy $w.vars
	set $var_list ""
    }
    method makeFrames {w} {
	$this buildColorMenus $w
	$this buildVarFrame $w
    }
    method graph {name} {
    }
    # this sets all the variables, var_root_j, with val
    method mat_sel_sub { var_root number val} {
	for {set j 0} { $j < $number} {incr j} {
	    set tail "_$j"
	    set $var_root$tail $val
	}
    }
    # called when SelAll or SelNone is evoked from the top level
    method mat_sel { var_root number val type} {
	for {set j 0} { $j < $number} {incr j} {
	    set tail "_$j"
	    switch $type {
		"matrix3" {
		    mat_sel_sub $var_root$tail $num_m_type $val
		}
		"vector" {
		    mat_sel_sub $var_root$tail $num_v_type $val
		}
		"scaler" {
		    set $var_root$tail $val
		}
	    }
	}
    }	
    # creates the material selection menu
    # it generates sub menus for matrix3 and vector types
    method make_mat_menu {w mat i type} {
	set fname "$w.mat$i"
	menubutton $fname -text "Material" \
		-menu $fname.list -relief groove
	pack $fname -side right -padx 2 -pady 2
	
	menu $fname.list
	set var_o $i
	append var_o "_o[set $this-var_orientation]"
	$fname.list add command -label "Sel All" \
		-command "$this mat_sel $this-matvar_$var_o $mat 1 $type"
	$fname.list add command -label "Sel None" \
		-command "$this mat_sel $this-matvar_$var_o $mat 0 $type"
	for {set j 0} { $j < $mat} {incr j} {
	    set var $var_o
	    append var "_$j"
#	    puts "***var = $var"
	    if {$type == "matrix3"} {
		$fname.list add cascade -label "Mat $j" \
			-menu $fname.list.types$var
		menu $fname.list.types$var
		$fname.list.types$var add command -label "Sel All" \
			-command "$this mat_sel_sub $this-matvar_$var \
			$num_m_type 1"
		$fname.list.types$var add command -label "Sel None" \
			-command "$this mat_sel_sub $this-matvar_$var \
			$num_m_type 0"
		for {set k 0} { $k < $num_m_type} {incr k} {
		    set var2 $var
		    append var2 "_$k"
		    $fname.list.types$var add checkbutton \
			    -variable $this-matvar_$var2 \
			    -label [lindex $matrix_types $k]
		}
	    } elseif {$type == "vector"} {
		$fname.list add cascade -label "Mat $j" \
			-menu $fname.list.types$var
		menu $fname.list.types$var
		$fname.list.types$var add command -label "Sel All" \
			-command "$this mat_sel_sub $this-matvar_$var \
			$num_v_type 1"
		$fname.list.types$var add command -label "Sel None" \
			-command "$this mat_sel_sub $this-matvar_$var \
			$num_v_type 0"
		for {set k 0} { $k < $num_v_type} {incr k} {
		    set var2 $var
		    append var2 "_$k"
		    $fname.list.types$var add checkbutton \
			    -variable $this-matvar_$var2 \
			    -label [lindex $vector_types $k]
		}
	    } elseif {$type == "scaler"} {
		$fname.list add checkbutton \
			-variable $this-matvar_$var \
			-label "Mat $j"
	    }
	}
    }
    method graphbutton {name var_index num_mat type} {
	set val_list {}
	set num_vals 0
	set var_root $this-matvar_$var_index
	append var_root "_o[set $this-var_orientation]"
#	puts "var_root = $var_root"
	# loop over all the materials	
	for {set j 0} { $j < $num_mat} {incr j} {
	    set mat_root $var_root
	    append mat_root "_$j"
	    switch $type {
		"matrix3" {
		    for {set k 0} { $k < $num_m_type} {incr k} {
			set tail "_$k"	
			if {[set $mat_root$tail] != 0} {
			    lappend val_list "$j" [lindex $matrix_types $k]
			    incr num_vals
			}
		    }
		}
		"vector" {
		    for {set k 0} { $k < $num_v_type} {incr k} {
			set tail "_$k"
			if {[set $mat_root$tail] != 0} {
			    lappend val_list "$j" [lindex $vector_types $k]
			    incr num_vals
			}
		    }
		}
		"scaler" {
		    if {[set $mat_root] != 0} {
			lappend val_list "$j" "invalid"
			incr num_vals
		    }
		}
	    }
	}
#	puts "Calling $this-c graph"
#	puts "name      = $name"
#	puts "num_mat   = $num_mat"
#	puts "var_index = $var_index"
#	puts "num_vals  = $num_vals"
#	puts "val_list  = $val_list"
	if {[llength $val_list] > 0} {
	    set call "$this-c graph $name $var_index $num_vals"
	    for {set i 0} { $i < [llength $val_list] } { incr i } {
		set insert [lindex $val_list $i]
		append call " $insert"
	    }
#	    puts "call =  $call"
	    eval $call
	}
    }
    method addVar {w name mat type i} {
	set fname "$w.var$i"
	frame $fname
	pack $fname -side top -fill x -padx 2 -pady 2

	label $fname.label -text "$name"
	pack $fname.label -side left -padx 2 -pady 2

	button $fname.button -text "Graph" -command "$this graphbutton $name $i $mat $type"
	pack $fname.button -side right -padx 2 -pady 2

	if {$mat > $num_materials} {
	    set num_materials $mat
#	    puts "num_materials is now $num_materials"
	}

	make_mat_menu $fname $mat $i $type
    }
    method setVar_list { args } {
	set var_list $args
	puts "var_list is now $var_list"
    }
    method setMat_list { args } {
	set mat_list $args
	puts "mat_list is now $mat_list"
    }
    method setType_list { args } {
	set type_list $args
	puts "type_list is now $type_list"
    }
    method buildVarFrame {w} {
	if {[llength $var_list] > 0} {
	    frame $w.vars -borderwidth 3 -relief ridge
	    pack $w.vars -side top -fill x -padx 2 -pady 2
	    
#	    puts "var_list length [llength $var_list]"
	    for {set i 0} { $i < [llength $var_list] } { incr i } {
		set newvar [lindex $var_list $i]
		set newmat [lindex $mat_list $i]
		set newtype [lindex $type_list $i]
		addVar $w.vars $newvar $newmat $newtype $i
	    }
	}
    }
    method buildColorMenus {w} {
	set n "$this-c needexecute "
	set i 0
	set b 1
	if {[set $this-nl] > 0} {
	    # color menu stuff
	    frame $w.colormenus -borderwidth 3 -relief ridge
	    pack $w.colormenus -side top -fill x -padx 2 -pady 2
	    
	    # set up the stuff for the grid colors
	    frame $w.colormenus.gridcolor
	    pack $w.colormenus.gridcolor -side left -fill y -padx 2 -pady 2
	    
	    setup_color $w.colormenus.gridcolor "grid" [set $this-nl] $n
	    
	    # set up the stuff for the node colors
	    frame $w.colormenus.nodecolor
	    pack $w.colormenus.nodecolor -side right -fill y -padx 2 -pady 2
	    
	    setup_color $w.colormenus.nodecolor "node" [set $this-nl] $n
	}
    }
    method make_entry {w text v c} {
	frame $w
	label $w.l -text "$text"
	pack $w.l -side left
	entry $w.e -textvariable $v -state disabled
	bind $w.e <Return> $c
	pack $w.e -side right
    }
    method ui {} {
	set w .ui[modname]
	if {[winfo exists $w]} {
	    wm deiconify $w
	    raise $w
	    return;
	}
	toplevel $w
	wm minsize $w 300 190
	set n "$this-c needexecute "

	#selection stuff
	frame $w.o
	pack $w.o -side top -fill x -padx 2 -pady 2

	frame $w.o.select
	pack $w.o.select -side left -fill x -padx 2 -pady 2
	
	checkbutton $w.o.select.plane_on -text "Selection plane" -variable \
		$this-plane_on -command $n
	pack $w.o.select.plane_on -side top -anchor w -pady 2 -ipadx 3

	checkbutton $w.o.select.nselect_on -text "Node Select" -variable \
		$this-node_select_on -command $n
	pack $w.o.select.nselect_on -side top -anchor w -pady 2 -ipadx 3

	radiobutton $w.o.select.node -variable $this-var_orientation \
		-command $n -text "Node Centered" -value 0
	pack $w.o.select.node -side top -anchor w -pady 2 -ipadx 3

	radiobutton $w.o.select.cell -variable $this-var_orientation \
		-command $n -text "Cell Centered" -value 1
	pack $w.o.select.cell -side top -anchor w -pady 2 -ipadx 3
	
	radiobutton $w.o.select.face -variable $this-var_orientation \
		-command $n -text "Face Centered" -value 2
	pack $w.o.select.face -side top -anchor w -pady 2 -ipadx 3
	
	button $w.o.select.findxy -text "Find XY" -command "$this-c findxy"
	pack $w.o.select.findxy -pady 2 -side top -ipadx 3 -anchor w

	button $w.o.select.findyz -text "Find YZ" -command "$this-c findyz"
	pack $w.o.select.findyz -pady 2 -side top -ipadx 3 -anchor w

	button $w.o.select.findxz -text "Find XZ" -command "$this-c findxz"
	pack $w.o.select.findxz -pady 2 -side top -ipadx 3 -anchor w

	# right side
	frame $w.o.r
	pack $w.o.r -side right -fill x -padx 2 -pady 2

	# sphere stuff
	frame $w.o.r.sphere
	pack $w.o.r.sphere -side top -anchor w -fill x -padx 2 -pady 2

	set r [expscale $w.o.r.sphere.radius -label "Radius:" \
		-orient horizontal -variable $this-radius -command $n ]
	pack $w.o.r.sphere.radius -side top -fill x

	scale $w.o.r.sphere.polygons -label "Polygons:" -orient horizontal \
	    -variable $this-polygons -command $n \
	    -from 8 -to 400 -tickinterval 392
	pack $w.o.r.sphere.polygons -side top -fill x

	# node ID
	make_entry $w.o.r.nodel "level index:" $this-index_l $n
	pack $w.o.r.nodel -side top -fill x -padx 2 -pady 2
	make_entry $w.o.r.nodex "x index:" $this-index_x $n
	pack $w.o.r.nodex -side top -fill x -padx 2 -pady 2
	make_entry $w.o.r.nodey "y index:" $this-index_y $n
	pack $w.o.r.nodey -side top -fill x -padx 2 -pady 2
	make_entry $w.o.r.nodez "z index:" $this-index_z $n
	pack $w.o.r.nodez -side top -fill x -padx 2 -pady 2

	makeFrames $w

	# close button
	button $w.close -text "Close" -command "wm withdraw $w"
	pack $w.close -side bottom -expand yes -fill x
#	button $w.gtest -text "Graph test" -command "$this graph_test"
#	pack $w.gtest -side bottom -expand yes -fill x
    }
    method reset_var_val {} {
	set var_val_list {}
    }
    method set_time { args } {
	set time_list $args
	puts "time_list =  $time_list"
    }
    method set_var_val { args } {
	set val_list $args
	lappend var_val_list $val_list
	puts "Args to set_var_val were $args"
#	puts "New var_val_list: $var_val_list"
    }
    method get_color { index } {
	set color_scheme {
	    { 255 0 0}  { 255 102 0}
	    { 255 204 0}  { 255 234 0}
	    { 204 255 0}  { 102 255 0}
	    { 0 255 0}    { 0 255 102}
	    { 0 255 204}  { 0 204 255}
	    { 0 102 255}  { 0 0 255}}
	#set color_scheme { {255 0 0} {0 255 0} {0 0 255} }
	set incr {}
	set upper_bounds [expr [llength $color_scheme] -1]
	for {set j 0} { $j < $upper_bounds} {incr j} {
	    set c1 [lindex $color_scheme $j]
	    set c2 [lindex $color_scheme [expr $j + 1]]
	    set incr_a {}
	    lappend incr_a [expr [lindex $c2 0] - [lindex $c1 0]]
	    lappend incr_a [expr [lindex $c2 1] - [lindex $c1 1]]
	    lappend incr_a [expr [lindex $c2 2] - [lindex $c1 2]]
	    lappend incr $incr_a
	}
	lappend incr {0 0 0}
#	puts "incr = $incr"
	set step [expr $num_colors / [llength $color_scheme]]
	set ind [expr $index % $num_colors] 
	set i [expr $ind / $step]
	set im [expr double($ind % $step)/$step]
#	puts "i = $i  im = $im"
	set curr_color [lindex $color_scheme $i]
	set curr_incr [lindex $incr $i]
#	puts "curr_color = $curr_color, curr_incr = $curr_incr"
	set r [expr [lindex $curr_color 0]+round([lindex $curr_incr 0] * $im)] 
	set g [expr [lindex $curr_color 1]+round([lindex $curr_incr 1] * $im)] 
	set b [expr [lindex $curr_color 2]+round([lindex $curr_incr 2] * $im)] 
	set c [format "#%02x%02x%02x" $r $g $b]
#	puts "r=$r, g=$g, b=$b, c=$c"
	return $c
    }
    method graph_data { id var args } {
	set w .graph[modname]$graph_id
        if {[winfo exists $w]} { 
            destroy $w 
	}
	toplevel $w
	incr graph_id
#	wm minsize $w 300 300

#	puts "id = $id"
#	puts "var = $var"
#	puts "args = $args"

	button $w.close -text "Close" -command "destroy $w"
	pack $w.close -side bottom -anchor s -expand yes -fill x
	
	blt::graph $w.graph -title "Plot of $var" \
		-height 250 -plotbackground gray99

	set max 1e-10
	set min 1e+10

	#seperate the materials from the types
	set args_mat {}
	set args_type {}
	for {set i 0} { $i < [llength $args] } {incr i} {
	    lappend args_mat [lindex $args $i]
	    incr i
	    lappend args_type [lindex $args $i]
	}
#	puts "args_mat = $args_mat"
#	puts "args_type = $args_type"
	
#	puts "length of var_val_list = [llength $var_val_list]"
#	puts "length of args_mat     = $args_mat"
	for { set i 0 } { $i < [llength $var_val_list] } {incr i} {
	    set mat_vals [lindex $var_val_list $i]
#	    puts "mat_vals = $mat_vals"
	    for { set j 0 } { $j < [llength $mat_vals] } {incr j} {
		set val [lindex $mat_vals $j]
		if { $max < $val } { set max $val }
		if { $min > $val } { set min $val }
	    }
	}
	
	if { ($max - $min) > 1000 || ($max - $min) < 1e-3 } {
	    $w.graph yaxis configure -logscale true -title $var
	} else {
	    $w.graph yaxis configure -title $var
	}
	
	$w.graph xaxis configure -title "Timestep" -loose true
	
	set vvlist_length [llength $var_val_list]
#	puts "length of var_val_list = [llength $var_val_list]"
	for { set i 0 } { $i < $vvlist_length } {incr i} {
#	    puts "adding line"
	    set mat_index  [lindex $args_mat $i]
	    set mat_type [lindex $args_type $i]
	    set mat_vals [lindex $var_val_list $i]
	    if {$i == 0} {
		set color_ind 0
	    } else {
		set color_ind [expr round(double($i) / \
			($vvlist_length-1) * ($num_colors -1))]
	    }
	    #	    set mat_vals [lindex $var_val_list $mat_index]
	    #	    set color_ind [expr round(double($mat_index) / ($num_materials-1) * ($num_colors - 1))]
	    set line_name "Material_$mat_index"
	    if {$mat_type != "invalid"} {
		append line_name "_$mat_type"
	    }
	    $w.graph element create $line_name -linewidth 2 -pixels 3 \
		    -color [$this get_color $color_ind] \
		    -xdata $time_list -ydata $mat_vals
	}
	
	pack $w.graph
    }
    method graph_test {} {
	set w .graph[modname]test
        if {[winfo exists $w]} { 
            destroy $w 
	}
	toplevel $w
#	wm minsize $w 300 300

	button $w.close -text "Close" -command "destroy $w"
	pack $w.close -side bottom -anchor s -expand yes -fill x
	
	blt::graph $w.graph -title "Test graph" -height 250 \
		-plotbackground gray99

	set x1 {0 1 2 3 4 5}
	set y1 {0 1 2 3 4}
	set x2 {0 1 2 3 4 5}
	set y2 {4 3 2 1 0}
	set x3 {0 1 2 3 4 5}
	set y3 {1 4 2 0 3}

	$w.graph yaxis configure -title "y-values"
	$w.graph xaxis configure -title "x-values" -loose true

	$w.graph element create line1 -linewidth 2 -pixels 3 \
		-color [$this get_color 0]  -xdata $x1 -ydata $y1
	$w.graph element create line2 -linewidth 2 -pixels 3 \
		-color [$this get_color 20] -xdata $x2 -ydata $y2
	$w.graph element create line3 -linewidth 2 -pixels 3 \
		-color [$this get_color 40] -xdata $x3 -ydata $y3

	$w.graph element create line4 -linewidth 2 -pixels 3 \
		-color [$this get_color 60]  -xdata $x1 -ydata $y1
	$w.graph element create line5 -linewidth 2 -pixels 3 \
		-color [$this get_color 80] -xdata $x2 -ydata $y2
	$w.graph element create line6 -linewidth 2 -pixels 3 \
		-color [$this get_color 100] -xdata $x3 -ydata $y3

	$w.graph element create line7 -linewidth 2 -pixels 3 \
		-color [$this get_color 120]  -xdata $x1 -ydata $y1
	$w.graph element create line8 -linewidth 2 -pixels 3 \
		-color [$this get_color 140] -xdata $x2 -ydata $y2
	$w.graph element create line9 -linewidth 2 -pixels 3 \
		-color [$this get_color 160] -xdata $x3 -ydata $y3

	$w.graph element create line10 -linewidth 2 -pixels 3 \
		-color [$this get_color 180]  -xdata $x1 -ydata $y1
	$w.graph element create line11 -linewidth 2 -pixels 3 \
		-color [$this get_color 200] -xdata $x2 -ydata $y2
	$w.graph element create line12 -linewidth 2 -pixels 3 \
		-color [$this get_color 220] -xdata $x3 -ydata $y3

	$w.graph element create line13 -linewidth 2 -pixels 3 \
		-color [$this get_color 240]  -xdata $x1 -ydata $y1
	$w.graph element create line14 -linewidth 2 -pixels 3 \
		-color [$this get_color 280] -xdata $x2 -ydata $y2
	$w.graph element create line15 -linewidth 2 -pixels 3 \
		-color [$this get_color 1100] -xdata $x3 -ydata $y3

	pack $w.graph
    }
}	
	
    

    









