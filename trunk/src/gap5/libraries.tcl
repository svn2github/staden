package require Plotchart

# Registers the xcombobox as a suitable editing tool for tablelist
;proc tablelist_add_xcombobox {} {
    set name xcombobox
    array set ::tablelist::editWin [list \
	$name-creationCmd   "xcombobox %W" \
	$name-putValueCmd   "%W set %T" \
	$name-getValueCmd   "%W get" \
	$name-putTextCmd    "%W set %T" \
	$name-getTextCmd    "%W get" \
	$name-putListCmd    "" \
	$name-getListCmd    "" \
	$name-selectCmd     "" \
	$name-invokeCmd     {[%W subwidget button] invoke} \
	$name-fontOpt       -font \
	$name-useFormat     1 \
	$name-useReqWidth   0 \
	$name-usePadX       1 \
	$name-isEntryLike   1 \
	$name-focusWin      {[%W subwidget entry]} \
	$name-reservedKeys  {Left Right Up Down}
       ]
}

tablelist_add_xcombobox

;proc ListLibrariesPopulate {io w} {
    # Clear any existing data
    $w selection clear 0 end
    $w delete 0 end

    set db [$io get_database]
    set nc [$db get_num_libraries]
    for {set i 0} {$i < $nc} {incr i} {
	set lib [$io get_library [$db get_library_rec $i]]
	$lib update_stats
	set name    [$lib get_name]
	set orient  [$lib get_orient]
	set type    [$lib get_machine]
	set mean    [$lib get_insert_size]
	set sd      [$lib get_insert_sd]
	set count   [$lib get_count]
	set type [lindex [list unknown Sanger Illumina SOLiD 454 Helicos IonTorrent PacBio ONT] $type]
	
	# Find most likely library orientation
	set max 0
	set o 0
	set tot 0
	for {set j 0} {$j < 3} {incr j} {
	    if {$max < [lindex $count $j]} {
		set max [lindex $count $j]
		set o $j
	    }
	    incr tot [lindex $count $j]
	}
	set mean [lindex $mean $o]
	set sd   [expr {int([lindex $sd   $o]+.5)}]
	set or   [lindex [list "-> <-" "<- ->" "-> -> / <- <-"] $o]

	$w insert end [list $i $name $tot $type $mean $sd $or]
	$lib delete
    }

    if {[$w sortcolumn] != -1} {
        $w sortbycolumn [$w sortcolumn] -[$w sortorder]
    }

    #$db delete
}

proc ListLibraries_editStartCmd {tbl row col text} {
    set w [$tbl editwinpath]

    puts [info level [info level]]

    switch $col {
	3 {
	    # Machine type
            $w configure -values [list unknown Sanger Illumina SOLiD 454 Helicos IonTorrent PacBio ONT]
	    $w configure -fixed_list 1
	}
    }

    return $text
}

proc ListLibraries_editEndCmd {io tbl row col text} {
    set w [$tbl editwinpath]

    puts [info level [info level]]

    set db [$io get_database]
    set rec [$db get_library_rec $row]
    set lib [$io get_library $rec]
    puts "Processing [$lib get_name]"

    switch $col {
       1 {
	   # Entry name
	   puts [$w get]
	   $lib set_name [$w get]
       }

       3 {
	   # Machine type
	   set id [lsearch [list unknown Sanger Illumina SOLiD 454 Helicos IonTorrent PacBio ONT] [$w get]]
	   puts [$w get]:$id
	   $lib set_machine_type $id
       }
    }

    $lib delete

    $io flush

    return $text
}

proc ListLibraries {io} {
    set t [xtoplevel .list_libraries]
    if {$t == ""} return
    wm title $t "List Libraries"
    wm geometry $t 550x600

    # Menus
    global list_libraries_menu
    menu $t.menubar
    $t configure -menu $t.menubar
    create_menus $list_libraries_menu $t.menubar

    # Create our tablelist
    tablelist $t.list \
	-columns {7 Index 15 "Name" 10 "Pair count" 8 "Type" 10 "Insert size" 6 "s.d." 15 "Orientation"} \
        -labelcommand tablelist::sortByColumn \
	-exportselection 0 \
	-stretch 0 \
        -yscrollcommand [list $t.yscroll set] \
	-editstartcommand ListLibraries_editStartCmd \
	-editendcommand "ListLibraries_editEndCmd $io"

    $t.list columnconfigure 0 -sortmode integer
    $t.list columnconfigure 1 -editable yes -editwindow entry
    $t.list columnconfigure 2 -sortmode integer
    $t.list columnconfigure 3 -editable yes -editwindow xcombobox
    $t.list columnconfigure 4 -sortmode integer
    $t.list columnconfigure 5 -sortmode real
    
    # Add scrollbars
    scrollbar $t.yscroll -command "$t.list yview"

    # Canvas
    canvas $t.c
    scrollbar $t.x -orient horiz -command "$t.c xview"
    scrollbar $t.y -orient vert  -command "$t.c yview"
    $t.c configure -xscrollcommand "$t.x set"
    $t.c configure -yscrollcommand "$t.y set"
    

    # Layout
    grid columnconfigure $t 0 -weight 1
    grid rowconfigure $t 0 -weight 1
    grid rowconfigure $t 1 -weight 1
    grid $t.list $t.yscroll -sticky nsew
    grid $t.c $t.y -sticky nsew
    grid $t.x -sticky nsew

    # Add some bindings
#    bind [$t.list bodypath] <<use>> "+ListLibrariesDetailed $io %W $t.list %x %y"
    bind [$t.list bodypath] <<select>> "+ListLibrariesDetailed $io %W $t %x %y"

    # Populate the list
    ListLibrariesPopulate $io $t.list
}


# Double-click callback from the library list.
;proc ListLibrariesDetailed {io w t x y} {
    # Move x,y to parent tablelist coordinates
    incr y [lindex [split [winfo geometry $w] +] 2]
    incr x [lindex [split [winfo geometry $w] +] 1]

    set row [$t.list containing $y]
    if {$row == -1} return

    incr row
    set idx [lindex [$w get $row.0 $row.end] 0]
    set db [$io get_database]
    set rec [$db get_library_rec $idx]

    draw_lib $io $t.c $rec
}

# Detailed data about a library including the distribution graph
;proc draw_lib {io c rec} {
    if {![winfo exists $c]} {
	#if {[xtoplevel $w] == ""} return
	scrollbar $w.xs -orient horiz -command "$c xview"
	scrollbar $w.ys -orient vert  -command "$c yview"
	canvas $c \
	    -xscrollcommand "$w.xs set" \
	    -yscrollcommand "$w.ys set" \
	    -height 400 \
	    -width 600
	grid rowconfigure    $w 0 -weight 1 -minsize 0
	grid columnconfigure $w 0 -weight 1 -minsize 0
	grid $c    -row 0 -column 0 -sticky nsew
	grid $w.xs -row 1 -column 0 -sticky nsew
	grid $w.ys -row 0 -column 1 -sticky nsew
    }

    # Get coordinate list
    set lib [$io get_library $rec]
    set dist [$lib get_dist]

    set l_0 {}
    set l_1 {}
    set l_2 {}
    set maxy 1
    set maxx 1
    set x_at_maxy 1
    for {set j 0} {$j < 3} {incr j} {
	set last_x 0
	set last_y 0
	set col [lindex {black orange3 blue} $j]
	foreach {x y} [lindex $dist $j] {
	    if {$last_x == 0} {set last_x $x; continue}
	    #set y [expr {log($y+1)}]
	    lappend l_$j $x $y
	    #$c create line $last_x $last_y $x $y -tags type_$j -fill $col -width 3
	    set last_x $x
	    set last_y $y
	    if {$maxy < $y} {set maxy $y; set x_at_maxy $x}
	    if {$maxx < $x} {set maxx $x}
	}
    }

    set bsize [expr {int(log($x_at_maxy)/log(2)-7)}]
    if {$bsize < 1} {set bsize 1}
    set maxy [expr {$maxy * $bsize*1.05}]

    # Fine tune X range
    set minx $maxx
    set maxx 1
    set y100 [expr {$maxy / 100}]
    foreach {x y} "$l_0 $l_1 $l_2" {
	if {$y >= $y100} {
	    if {$x < $minx} {set minx $x}
	    if {$x > $maxx} {set maxx $x}
	}
    }

    # Create the plots
    incr minx -20
    if {$minx < 1} {set minx 1}
    incr maxx 20

    foreach {xn xi} [nice_num $minx [expr {$maxx-$minx+1}] 0] {break}
    foreach {yn yi} [nice_num 0 $maxy 0] {break}

    $c delete all
    set s [::Plotchart::createXYPlot $c \
	       [list $xn $maxx $xi] \
	       [list 0 $maxy $yi]]
    $s xtext "Insert size"
    $s ytext "Count"

    # Plot it
    for {set j 0} {$j < 3} {incr j} {
	$s dataconfig lib_${rec}_$j \
	    -colour [lindex {red blue magenta} $j]
	$s legend lib_${rec}_$j \
	    [lindex {"-->  <--" "<--  -->" "-> -> / <- <-"} $j]
	foreach {x y} [set l_$j] {
	    $s plot lib_${rec}_$j $x [expr {$y*$bsize}]
	}
    }

    $c itemconfigure data -width 2
    $c configure -scrollregion [$c bbox all]

    $lib delete
}
