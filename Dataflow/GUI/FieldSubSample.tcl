#
#  For more information, please see: http://software.sci.utah.edu
# 
#  The MIT License
# 
#  Copyright (c) 2004 Scientific Computing and Imaging Institute,
#  University of Utah.
# 
#  License for the specific language governing rights and limitations under
#  Permission is hereby granted, free of charge, to any person obtaining a
#  copy of this software and associated documentation files (the "Software"),
#  to deal in the Software without restriction, including without limitation
#  the rights to use, copy, modify, merge, publish, distribute, sublicense,
#  and/or sell copies of the Software, and to permit persons to whom the
#  Software is furnished to do so, subject to the following conditions:
# 
#  The above copyright notice and this permission notice shall be included
#  in all copies or substantial portions of the Software.
# 
#  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
#  OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
#  THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
#  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
#  DEALINGS IN THE SOFTWARE.
#


# GUI for FieldSubSample module
# by Allen R. Sanderson
# March 2002

# This GUI interface is for sub sampling a topologically structured field.

itcl_class SCIRun_FieldsCreate_FieldSubSample {
    inherit Module
    constructor {config} {
        set name FieldSubSample
        set_defaults
    }

    method set_defaults {} {
	global power_app_command
	set    power_app_command ""

	global $this-wrap
	global $this-dims

	set $this-wrap 0
	set $this-dims 3

	for {set i 0} {$i < 3} {incr i 1} {
	    if { $i == 0 } {
		set index i
	    } elseif { $i == 1 } {
		set index j
	    } elseif { $i == 2 } {
		set index k
	    }

	    global $this-$index-dim
	    global $this-$index-start
	    global $this-$index-start2
	    global $this-$index-stop
	    global $this-$index-stop2
	    global $this-$index-stride
	    global $this-$index-stride2
	    global $this-$index-wrap

	    set $this-$index-dim      2
	    set $this-$index-start    0
	    set $this-$index-start2  "0"
	    set $this-$index-stop     1
	    set $this-$index-stop2   "1"
	    set $this-$index-stride   1
	    set $this-$index-stride2 "1"
	    set $this-$index-wrap     0

	    trace variable $this-$index-dim w "$this update_set_size_callback"
	}
	
	trace variable $this-dims w "$this update_set_size_callback"
    }

    method set_power_app_cmd { cmd } {
	global power_app_command
	set power_app_command $cmd
    }

    method ui {} {

	global $this-wrap
	global $this-dims

	set tmp 0.0

        set w .ui[modname]
        if {[winfo exists $w]} {
            raise $w
            return
        }

	if { [set $this-wrap] } {
	    set wrap normal
	} else {
	    set wrap disable
	}

        toplevel $w

	frame $w.main

	frame $w.main.l
	label $w.main.l.direction -text "Index"  -width 5 -anchor w -just left
	label $w.main.l.start     -text "Start"  -width 5 -anchor w -just left
	label $w.main.l.stop      -text "Stop"   -width 5 -anchor w -just left
	label $w.main.l.stride    -text "Stride" -width 6 -anchor w -just left
	label $w.main.l.wrap      -text "Wrap"   -width 4 -anchor w -just left

	pack $w.main.l.direction -side left -padx  20
	pack $w.main.l.start     -side left -padx  70
	pack $w.main.l.stop      -side left -padx 110
	pack $w.main.l.stride    -side left -padx  40
	pack $w.main.l.wrap      -side left

#	grid $w.main.l.direction $w.main.l.start $w.main.l.stop $w.main.l.stride $w.main.l.wrap

	for {set i 0} {$i < 3} {incr i 1} {
	    if { $i == 0 } {
		set index i
	    } elseif { $i == 1 } {
		set index j
	    } elseif { $i == 2 } {
		set index k
	    }

	    global $this-$index-dim
	    global $this-$index-start
	    global $this-$index-start2
	    global $this-$index-stop
	    global $this-$index-stop2
	    global $this-$index-stride
	    global $this-$index-stride2
	    global $this-$index-wrap

	    # Update the sliders to have the new end values.
	    if { [set $this-wrap] == 0 } {    
		set $this-$index-wrap 0
	    }

	    if [set $this-$index-wrap] {
		set start_val 0
		set stop_val [expr [set $this-$index-dim] - 1]
	    } else {
		set start_val 1
		set stop_val [expr [set $this-$index-dim] - 2]
	    }

	    frame $w.main.$index

	    label $w.main.$index.l -text " $index :" -width 3 -anchor w -just left

	    pack $w.main.$index.l -side left

	    scaleEntry4 $w.main.$index.start \
		0 $stop_val 200 \
		$this-$index-start $this-$index-start2 $index

	    scaleEntry2 $w.main.$index.stop \
		$start_val [expr [set $this-$index-dim] - 1] 200 \
		$this-$index-stop $this-$index-stop2

	    scaleEntry2 $w.main.$index.stride \
		1 [expr [set $this-$index-dim] - 1] 100 $this-$index-stride $this-$index-stride2

	    checkbutton $w.main.$index.wrap -variable $this-$index-wrap \
		    -state $wrap -disabledforeground "" \
		    -command "$this wrap $index"

	    pack $w.main.$index.l $w.main.$index.start $w.main.$index.stop \
		    $w.main.$index.stride $w.main.$index.wrap -side left
#	    grid $w.main.$index.l $w.main.$index.start $w.main.$index.stop 
#		    $w.main.$index.stride $w.main.$index.wrap
	}

	if { [set $this-dims] == 3 } {
	    pack $w.main.l $w.main.i $w.main.j $w.main.k -side top -padx 10 -pady 5
	} elseif { [set $this-dims] == 2 } {
	    pack $w.main.l $w.main.i $w.main.j -side top -padx 10 -pady 5
	} elseif { [set $this-dims] == 1 } {
	    pack $w.main.l $w.main.i -side top -padx 10 -pady 5	    
	}

	pack $w.main -side top -fill x -expand 1
	
	global power_app_command

	if { [in_power_app] } {
	    makeSciButtonPanel $w $w $this -no_execute -no_close -no_find \
		"\"Close\" \"wm withdraw $w; $power_app_command\" \"Hides this GUI\""
	} else {
	    makeSciButtonPanel $w $w $this
	}
	 
	moveToCursor $w
    }

    method scaleEntry2 { win start stop length var1 var2 } {
	frame $win 
	pack $win -side top -padx 5

	scale $win.s -from $start -to $stop -length $length \
	    -variable $var1 -orient horizontal -showvalue false \
	    -command "$this updateSliderEntry $var1 $var2"

	entry $win.e -width 4 -text $var2

	bind $win.e <KeyRelease> "$this manualSliderEntry $start $stop $var1 $var2"

	pack $win.s -side left
	pack $win.e -side bottom -padx 5
    }

    method updateSliderEntry {var1 var2 someUknownVar} {
	set $var2 [set $var1]
    }

    method manualSliderEntry { start stop var1 var2 } {

	if { ![string is integer [set $var2]] } {
	    set $var2 [set $var1] }

	if { [set $var2] < $start } {
	    set $var2 $start }
	
	if { [set $var2] > $stop } {
	    set $var2 $stop }
	
	set $var1 [set $var2]
    }


    method wrap { index } {

	global $this-$index-start
	global $this-$index-start2
	global $this-$index-stop
	global $this-$index-stop2
	global $this-$index-dim
	global $this-$index-wrap

	set w .ui[modname]

	if [ expr [winfo exists $w] ] {

	    # Update the sliders to have the new end values.

	    if [ set $this-$index-wrap ] {
		set start_val 0
		set stop_val  [expr [set $this-$index-dim] - 1]
	    } else {
		set start_val [expr [set $this-$index-start] + 1]
		set stop_val  [expr [set $this-$index-dim] - 2]
	    }

	    $w.main.$index.start.s configure -from 0 -to $stop_val
	    $w.main.$index.stop.s configure -from $start_val -to [expr [set $this-$index-dim] - 1]

	    bind $w.main.$index.start.e <KeyRelease> \
		"$this manualSliderEntry4 0 $stop_val $this-$index-start $this-$index-start2 $index"
	    bind $w.main.$index.stop.e  <KeyRelease> \
		"$this manualSliderEntry $start_val  [expr [set $this-$index-dim] - 1] $this-$index-stop $this-$index-stop2"
	}
    }

    method scaleEntry4 { win start stop length var1 var2 index } {
	frame $win 
	pack $win -side top -padx 5

	scale $win.s -from $start -to $stop -length $length \
	    -variable $var1 -orient horizontal -showvalue false \
	    -command "$this updateSliderEntry4 $index"

	entry $win.e -width 4 -text $var2

	bind $win.e <KeyRelease> \
	    "$this manualSliderEntry4 $start $stop $var1 $var2 $index"

	pack $win.s -side left
	pack $win.e -side bottom -padx 5
    }

    method updateSliderEntry4 { index someUknownVar } {

	global $this-$index-start
	global $this-$index-start2
	global $this-$index-stop
	global $this-$index-stop2

	wrap $index

	set $this-$index-start2 [set $this-$index-start]
	set $this-$index-stop2  [set $this-$index-stop]
    }

    method manualSliderEntry4 { start stop var1 var2 index } {

	if { ![string is integer [set $var2]] } {
	    set $var2 [set $var1] }

	if { [set $var2] < $start } {
	    set $var2 $start }
	
	if { [set $var2] > $stop } {
	    set $var2 $stop }
	
	set $var1 [set $var2]

	updateSliderEntry4 $index 0
    }

    method update_index { } {
	global $this-i-index
	global $this-j-index
	global $this-k-index

	global $this-i-index2
	global $this-j-index2
	global $this-k-index2

	set $this-i-index2 [set $this-i-index]
	set $this-j-index2 [set $this-j-index]
	set $this-k-index2 [set $this-k-index]
    }

    method update_set_size_callback { name1 name2 op } {
	set_size
    }

    method set_size { } {
	global $this-dims
	global $this-i-dim
	global $this-j-dim
	global $this-k-dim
	global $this-wrap

	set w .ui[modname]

	if [ expr [winfo exists $w] ] {
	    pack forget $w.main.i
	    pack forget $w.main.k
	    pack forget $w.main.j
	    
	    if { [set $this-dims] == 3 } {
		pack $w.main.l $w.main.i $w.main.j $w.main.k -side top -padx 10 -pady 5
	    } elseif { [set $this-dims] == 2 } {
		pack $w.main.l $w.main.i $w.main.j -side top -padx 10 -pady 5
	    } elseif { [set $this-dims] == 1 } {
		pack $w.main.l $w.main.i -side top -padx 10 -pady 5
	    }
	}

	for {set i 0} {$i < 3} {incr i 1} {
	    if { $i == 0 } {
		set index i
	    } elseif { $i == 1 } {
		set index j
	    } elseif { $i == 2 } {
		set index k
	    }

	    global $this-$index-start
	    global $this-$index-start2
	    global $this-$index-stop
	    global $this-$index-stop2
	    global $this-$index-stride
	    global $this-$index-stride2

	    set $this-$index-wrap 0

	    set stop_val1 [expr [set $this-$index-dim] - 1]
	    set stop_val2 [expr [set $this-$index-dim] - 2]

	    if [ expr [winfo exists $w] ] {

		if { [set $this-wrap ] } {
		    $w.main.$index.wrap configure -state normal
		} else {
		    $w.main.$index.wrap configure -state disabled 
		}

		# Update the sliders to the new bounds.
		$w.main.$index.start.s configure -from 0 -to $stop_val2
		$w.main.$index.stop.s  configure -from 0 -to $stop_val1
		$w.main.$index.stride.s  configure -from 1 -to $stop_val1

		bind $w.main.$index.start.e <KeyRelease> \
		    "$this manualSliderEntry4 0 $stop_val2 $this-$index-start $this-$index-start2 $index"
		bind $w.main.$index.stop.e  <KeyRelease> \
		    "$this manualSliderEntry  1 $stop_val1 $this-$index-stop $this-$index-stop2"
		bind $w.main.$index.stride.e  <KeyRelease> \
		    "$this manualSliderEntry  1 $stop_val1 $this-$index-stride $this-$index-stride2"
	    }

	    # Update the stop values to be at the initials values.
	    set $this-$index-start 0	    
	    set $this-$index-stop  $stop_val1
	    set $this-$index-stride  1

	    # Update the text values.
	    set $this-$index-start2 [set $this-$index-start]
	    set $this-$index-stop2  [set $this-$index-stop]
	    set $this-$index-stride2  [set $this-$index-stride]
	}
    }
}
