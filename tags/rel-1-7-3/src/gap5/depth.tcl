# A horribly hacky example of a tk canvas view of contig 0 from a tgap
# database. It's inefficient in caching data too, but despite that is a
# proof of concept that we're fast enough right now.
#
# Warning, this is full of hard coded Tk pathnames, Tcl global variables,
# assumptions of contig length, etc. All this ensures we'll not feel tempted
# to steal any of this hacked up demo code :-)

package require Tk

# Maximum width of template-dislay before we omit drawing it
set max_template_display_width 500000

#
# Constructs a container window for drawing horizontal tracks of data.
# It consists of standard controls like menus, horizontal scrollbars and
# a scale control along with zero or more children plots, named "tracks".
#
# Each track has it's own ID value and function associated with it.
# The 1.5Dplot passes high level events (like scroll in X) down to all
# suitable child tracks.
#
proc 1.5Dplot {w io {cnum {}}} {
    global $w

    # Create the window
    if {[toplevel $w] == ""} return

    set wid 800
    set hei 600

    wm geometry $w ${wid}x${hei}

    if {$cnum == ""} {
	set cnum [$io contig_order 0]
    }

    set ${w}(io) $io
    set ${w}(cnum) $cnum
    set ${w}(contig) [$io get_contig $cnum]
    set ${w}(length) [[set ${w}(contig)] get_length]
    set ${w}(tracks) {}
    set ${w}(move_status) 0
    set ${w}(yzoom) 1
    set ${w}(xorigin) 0
    set ${w}(width) $wid
    set ${w}(height) $hei
    set ${w}(border) 500
    set ${w}(x1)    0
    set ${w}(x2)    100
    set ${w}(last_x1) -9999999
    set ${w}(last_x2) -9999999
    set ${w}(ntracks) 0

    set ${w}(x1) 0
    set ${w}(x2) [expr {[set ${w}(x1)]+1000}]

    set xstart 0
    set ${w}(xorigin) [set ${w}(x1)]
    set ${w}(xzoom) [expr {([set ${w}(x2)]-[set ${w}(x1)]+1)/
			   double([set ${w}(width)])}]

    # X Zoom control
    scale $w.xzoom -from 1 -to 1000 -resolution 0.1 -orient horiz -label Zoom \
	-command "set_zoom $w"
    $w.xzoom set 16
    grid columnconfigure $w 0 -weight 1
    grid $w.xzoom -row 998 -sticky nsew
    grid rowconfigure $w 998 -weight 0

    # X scrollbar
    scrollbar $w.xscroll \
	-command "scrollx1.5 $w" \
	-repeatinterval 5 \
	-orient horiz
    $w.xscroll set \
	[expr {[set ${w}(x1)]/double([set ${w}(length)])}] \
	[expr {[set ${w}(x2)]/double([set ${w}(length)])}]

    grid $w.xscroll -row 999 -sticky nsew
    grid rowconfigure $w 999 -weight 0

    # An information label
    frame $w.label
    grid $w.label -row 1000 -sticky nsew
    label $w.label.l -textvariable ${w}(info)
    pack $w.label.l

    bind $w <5> "zoom1.5 $w 0 1 1.1"
    bind $w <4> "zoom1.5 $w 0 1 [expr {1/1.1}]"

    bind $w <Any-Configure> "resize1.5 $w"

    redraw_plot $w
}

# Resize events
proc resize1.5 {w} {
    global $w

    # We get rogue events when highlighting. Skip these
    if {[winfo width  $w] == [set ${w}(width)] && \
	[winfo height $w] == [set ${w}(height)]} return

    puts "\n==================="
    puts "resize to [winfo width $w]x[winfo height $w]"

    set ${w}(width) [winfo width $w]
    set ${w}(height) [winfo height $w]

    foreach id [set ${w}(tracks)] {
	set t $w.track$id
	global $t
	set width  [winfo width  [set ${t}(canvas)]]
	set height [winfo height [set ${t}(canvas)]]
	set ${t}(width)  $width
	set ${t}(height) $height
    }

    set ${w}(x2) [expr {[set ${w}(xzoom)]*[set ${w}(width)] +
			+[set ${w}(x1)] - 1}]
    redraw_plot $w
}

# Scroll in X
proc scrollx1.5 {w cmd args} {
    global $w

    set sbar $w.xscroll
    set clen [set ${w}(length)]

    puts ""
    #parray $w

    if {$cmd == "moveto"} {
	set xpos [expr {int([lindex $args 0]*$clen)}]
    } elseif {$cmd == "set_xpos"} {
	set xpos [lindex $args 0]
    } elseif {$cmd == "scroll"} {
	set xpos [expr {[lindex [$sbar get] 0]*$clen}]
	if {[lindex $args 1] == "pages"} {
	    set wid [expr {[c2x $w [set ${w}(width)]] - [c2x $w 0]}]
	    set xpos [expr {$xpos + $wid/2*[lindex $args 0]}]
	} else {
	    set wid [expr {[c2x $w [lindex $args 0]] - [c2x $w 0]}]
	    set xpos [expr {$xpos+5*$wid}]
	}
    } else {
	set xpos $cmd
    }

    set wid [expr {[set ${w}(x2)] - [set ${w}(x1)]}]
    set ${w}(x1) $xpos
    set ${w}(x2) [expr {$xpos+$wid}]

    # Update scrollbar
    $sbar set \
	[expr {[set ${w}(x1)]/double([set ${w}(length)])}] \
	[expr {[set ${w}(x2)]/double([set ${w}(length)])}]

    redraw_plot $w
}

# X zoom only via scalebox
proc set_zoom {w val} {
    global $w

    # Change the scale
    set wid [set ${w}(width)]
    set scale [set ${w}(xzoom)]
    set mid [c2x $w [expr {$wid/2.0}]]
    set ${w}(xzoom) [expr {pow(($val+3)/10,3)}]

    # Adjust xorigin to ensure $mid stays the same
    set ${w}(xorigin) [expr {$mid - ($wid/2)*[set ${w}(xzoom)]}]

    # Reset x1/x2
    set ${w}(x1) [c2x $w 0]
    set ${w}(x2) [c2x $w $wid]

    # Simulate a scrollbar movement to update size, position and redraw.
    set sbar $w.xscroll

    $sbar set \
	[expr {[set ${w}(x1)]/double([set ${w}(length)])}] \
	[expr {[set ${w}(x2)]/double([set ${w}(length)])}]

    redraw_plot $w
}

# Zoom the plot. z>1 => magnify/zoom in, z<1 => zoom out
# x and y are booleans to govern whether we want to zoom in x, y or both.
proc zoom1.5 {w x y z} {
    global $w

    if {$y} {
	set ${w}(yzoom) [expr {[set ${w}(yzoom)]/$z}]
    }

    if {$x} {
	#set ${w}(xzoom) [expr {$z*[set ${w}(xzoom)]}]
	set x1 [set ${w}(x1)]
	set x2 [set ${w}(x2)]
	set mid [expr {($x1+$x2)/2.0}]
	set wid [expr {$x2-$x1+1}]
	set wid [expr {$wid * $z}];
    
	puts "old x=$x1..$x2"

	set x1 [expr {int($mid-$wid/2)}]
	set x2 [expr {int($mid+$wid/2)}]
	set ${w}(x1) $x1
	set ${w}(x2) $x2

	puts "new x=$x1..$x2"

	$w.xscroll set \
	    [expr {$x1/double([set ${w}(length)])}] \
	    [expr {$x2/double([set ${w}(length)])}]

	set ${w}(xorigin) [expr $x1]
	set ${w}(xzoom)   [expr {($x2-$x1+1)/double([set ${w}(width)])}]
    }

    redraw_plot $w seq_depth
    return
}

#
# Functions to convert x/y base-level coordinates into canvas coordinates (c)
# and back again.
#
proc x2c {w pos} {
    global $w
    return [expr {($pos - [set ${w}(xorigin)]) / [set ${w}(xzoom)]}]
}

proc c2x {w pos} {
    global $w
    return [expr {$pos * [set ${w}(xzoom)] + [set ${w}(xorigin)]}]
}

# m* is the size of the margins.
# +ve = from that edge.
# -ve = from opposite edge
proc add_plot {w func height args} {
    puts [info level [info level]]

    global $w
    incr ${w}(ntracks)
    set tnum [set ${w}(ntracks)]
    lappend ${w}(tracks) $tnum
    set t $w.track$tnum
    global $t

    set ${t}(w)      $w
    set ${t}(num)    $tnum
    set ${t}(func)   $func
    set ${t}(height) $height

    # Heights
    if {$height < 0} {
	set height [expr {-$height}]
	set weight 1
	set ${t}(ys) $w.yscroll$tnum
	scrollbar $w.yscroll$tnum -orient vert -command "yscroll_plot $w $t"
	grid $w.yscroll$tnum -row $tnum -column 1 -sticky nsew
    } else {
	set weight 0
    }
    set ${t}(y1) 0
    set ${t}(y2) $height

    set ${t}(canvas) [eval canvas $t -height $height $args \
			  [list -yscrollcommand [list yscroll_plot_set $w $t]]]
    grid $t -row $tnum -sticky nsew
    grid rowconfigure $w $tnum -weight $weight

    #redraw_plot $w
}

proc yscroll_plot {w t cmd {opt1 {}} {opt2 {}}} {
    global $w $t
    puts [info level [info level]]

    set y1 [set ${t}(y1)]
    set h  [set ${t}(height)]
    switch $cmd {
	scroll {
	    set y1 [set ${t}(y1)]
	    if {$opt2 == "pages"} {
		incr y1 [expr {int($opt1*$h/2.0)}]
	    } else {
		incr y1 [expr {int($opt1*1)}]
	    }
	}

	moveto {
	    set y1 [expr {int($opt1 * [set ${t}(scroll_height)])}]
	}
    }

    set ${t}(y1) $y1
    set y2 [expr {$y1+$h}]
    set wid [set ${t}(width)]

    [set ${t}(canvas)] configure -scrollregion [list 0 $y1 $wid $y2]
    [set ${t}(ys)] set \
	[expr {$y1/[set ${t}(scroll_height)]}] \
	[expr {$y2/[set ${t}(scroll_height)]}]

    #redraw_plot $w [set ${t}(num)]
}

proc yscroll_plot_set {w t args} {
    puts [info level [info level]]
}

# The top-level redraw function. This works out what region is visible and
# fires off the individual redraw functions registered.
proc redraw_plot {w {track_types {}}} {
    global $w
    parray $w
    

    set bd [set ${w}(border)]

    puts "Redrawing"

    set x1 [set ${w}(x1)]
    set x2 [set ${w}(x2)]
    set ${w}(xorigin) [expr $x1]

    puts ""
    puts ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
    puts "xorigin=[set ${w}(xorigin)] xzoom=[set ${w}(xzoom)] yzoom=[set ${w}(yzoom)]"
    puts "base coord=$x1..$x2   canvas coord=[x2c $w $x1]..[x2c $w $x2]"
    puts ""

    if {$track_types != {}} {
	set tracks {}
	foreach id [set ${w}(tracks)] {
	    set t $w.track$id
	    global $t
	    if {[lsearch -exact $track_types [set ${t}(func)]] != -1} {
		lappend tracks $id
	    }
	}
    } else {
	set tracks [set ${w}(tracks)]
    }
    puts tracks=$tracks

    foreach id $tracks {
	set t $w.track$id
	global $t

	parray $t

	set d [set ${t}(canvas)]
	$d delete all
	set y1 [set ${t}(y1)]
	set y2 [expr {[set ${t}(y1)]+[set ${t}(height)]}]
	puts [time {[set ${t}(func)] $w $t $x1 $x2 $y1 $y2}]
    }

    set ${w}(last_x1) $x1
    set ${w}(last_x1) $x2
}


#
# Plots the sequence read depth.
#
proc seq_depth {w t x1 x2 y1 y2} {
    global $w $t

    puts [info level [info level]]

    set c    [set ${w}(contig)]
    set d    [set ${t}(canvas)]
    set wid  [set ${w}(width)]
    set yz   [set ${w}(yzoom)]

    set inc [expr {($x2-$x1+1)/double($wid)}]
    if {$inc < 1} {set inc 1}

    set l [$c read_depth [expr {int($x1)}] [expr {int($x2)}] $inc]
    set len [llength $l]
    puts ">>> Region $x1..$x2 y=$y1..$y2 in steps of $inc => $len items"

    set h [expr {$y2-$y1+1}]

    set wx1 [x2c $w $x1]
    set wx2 [x2c $w $x2]
    set wid [expr {$wx2-$wx1+1}]
    if {$wid < 0} {set wid 1}
    set bpv [expr {$len/double($wid)}]

    set line {}
    set max $h
    if {$bpv < 1} {
	for {set i 0} {$i < $len} {incr i} {
	    set avg [expr {[lindex $l $i]*$yz}]
	    if {$avg > $max} {set avg $max}
	    lappend line [expr {$wx1+$i/$bpv}] [expr {$h-$avg}]
	}
    } else {
	for {set i 0} {$i < $len} {incr i} {
	    set avg [expr {[lindex $l $i]*$yz}]
	    if {$avg > $max} {set avg $max}
	    lappend line [expr {$wx1+$i}] [expr {$h-$avg}]
	}
    }

    $d create line $line -fill black
    puts BPV=$bpv/[lrange $line 0 1]..[lrange $line [expr {[llength $line]-4}] end]
}

proc piarray {n} {
    upvar 1 $n a

    foreach name [lsort -integer [array names a]] {
	puts [format "${n}(%2d) %f" $name $a($name)]
    }
}

# Given a list of horizontal lines this distributes the elements around
# so they do not clash in Y
proc compute_max_y {w lines} {
    global $w

    set ymax 5
    set ymap [lindex [lindex $lines 0] 1]
    for {set i 1} {$i < $ymax} {incr i} {
	lappend ymap -999999
    }
    set ysize 1

    set all_y {}

    set n 1
    foreach l $lines {
	# First and last coords represent the range of lines
	set x1 [lindex $l 0]
	set x2 [lindex $l end-2]

	while {$x1 < [lindex $ymap [expr {$n%$ymax}]]} {
	    # Compute new and old modulus
	    set old_mod [expr {$n%$ymax}]
	    set new_mod [expr {$n%($ymax+1)}]
	    set mod_delta [expr {$new_mod - $old_mod}]

	    # Rotate if needed
	    if {$mod_delta > 0} {
		set list [concat \
			    [lrange $ymap [expr {$ymax-$mod_delta}] end] \
			    [lrange $ymap 0 [expr {$ymax-$mod_delta-1}]]]
	    } elseif {$mod_delta < 0} {
		set list [concat \
			    [lrange $ymap [expr {-$mod_delta}] end] \
			    [lrange $ymap 0 [expr {-$mod_delta-1}]]]
	    } else {
		set list $ymap
	    }


	    # Insert previous item
	    if {[llength $all_y] > $ymax} {
		set app [lindex $all_y end-$ymax]
		set all_y [lreplace $all_y end-$ymax end-$ymax]
	    } else {
		set app -99999999
	    }
	    set ymap [linsert $list $new_mod $app]
	    incr ymax
	}

	lappend all_y [expr {$x2+5}]
	lset ymap [expr {$n%$ymax}] [expr {$x2+5}]
	incr n
    }

    return $ymax
}

#
# Plots the sequence read depth.
#
proc seq_seqs {w t x1 x2 y1 y2} {
    global $w $t
    global max_template_display_width

    puts SEQ_SEQS_START,$y1..$y2

    set c    [set ${w}(contig)]
    set d    [set ${t}(canvas)]
    set wid  [set ${w}(width)]
    set yz   [set ${w}(yzoom)]
    set io   [set ${w}(io)]
    
    if {![info exists ${w}(init)]} {
	set ${w}(init) 1
	$d bind all <Any-Enter> "seq_seqs_bind $t $w $d"
    }

    set wx1 [x2c $w $x1]
    set wx2 [x2c $w $x2]
    set wid [expr {$wx2-$wx1+1}]

    if {[expr {$x2-$x1}] > $max_template_display_width} return

    set reads [$c seqs_in_range [expr {int($x1)}] [expr {int($x2+0.999)}]]
    puts "nreads between $x1..$x2=[llength $reads]"
    puts SEQ_SEQS_DRAW_START
    #set reads [lrange $reads 0 10000]
    catch {unset done}

    # Obtain the list of reads to plot and generate horizontal lines.
    # For now we have no X coordinates though
    set cn1 [set ${w}(cnum)]
    foreach r $reads {
	foreach {st1 en1 rec comp} $r {break}
	if {[info exists done($rec)]} continue

	set unordered 0
	set r1 [$io get_sequence $rec]
	set done($rec) 1

	set xst1 [x2c $w $st1]
	set xen1 [x2c $w $en1]

	set pair [$r1 get_pair]
	set dir1 [expr {[$r1 get_length] >= 0 ? "f" : "r"}]
	$r1 delete

	# Paired end
	if {$pair} {
	    set line {}
	    set r2 [$io get_sequence $pair]
	    set st2 [$r2 get_position]
	    set cn2 [$r2 get_contig]
	    set en2 [expr {$st2+abs([$r2 get_length])-1}]
	    set dir2 [expr {[$r2 get_length] >= 0 ? "f" : "r"}]
	    $r2 delete

	    # Temporary hack: only show long read-pairs
	    #if {$en2 - $st1 < 10000} continue

	    lappend line $xst1 $xen1 $dir1 $rec
	    set done($pair) 1

	    #if {$en2 >= $x1-100 && $st2 <= $x2+100 && $cn2 == $cn1} 
	    if {$cn2 == $cn1} {
		# Both semi-visible, so draw
		if {$st2 >= $st1} {
		    lappend line $xen1 [x2c $w $st2] t $rec
		    lappend line [x2c $w $st2] [x2c $w $en2] $dir2 $rec
		} else {
		    set line [concat [x2c $w $en2] $xst1 t $rec $line]
		    # The read at the other end won't be visible anyway as
		    # if it was then it'd have already processed it as
		    # an 'r1' abiove.
		    #lappend line [x2c $w $st2] [x2c $w $en2] $dir2 $rec
		    
		    # However in this case we won't be sorted in left to
		    # right order either, so set a flag indicating as
		    # such.
		    set unordered 1
		}
	    } else {
		# Other end is invisible, so draw dotted line
		if {$cn2 != $cn1} {
		    set mode D
		} else {
		    set mode d
		}
		if {$dir1 == "f"} {
		    lappend line $xen1 [expr {$xen1+20}] $mode $rec
		} else {
		    set line [concat [expr {$xst1-20}] $xst1 $mode $rec $line]
		}
	    }

	    # NB: line consists of multiple segments, but must be in
	    # left to right order
	} else {
	    set line [list $xst1 $xen1 s $rec]
	}
	set sz [expr {int(log([lindex $line end-2]-[lindex $line 0]))}]
	if {$unordered && [info exists lines_sz($sz)]} {
	    set p1 [lindex $line 0]
	    set c 0
	    foreach l $lines_sz($sz) {
		if {[lindex $l 0] > $p1} break;
		incr c
	    }
	    set lines_sz($sz) [linsert $lines_sz($sz) $c $line]
	} else {
	    lappend lines_sz($sz) $line
	}
    }

    # Plot
    set y 3
    # Approx .3sec here
    set yoff [set ${t}(y1)]

    set ysep [expr {10-log($x2-$x1)}]
    set ysep [expr {2+($ysep > 0 ? int($ysep) : 0)}]

    set ystart 0
    set ymax 0
    foreach sz [lsort -integer -decreasing [array names lines_sz]] {
	set ymax [compute_max_y $w $lines_sz($sz)]
	foreach l $lines_sz($sz) {
	    set yp [expr {($y%$ymax+$ystart)*$ysep+3}]
	    foreach {x1 x2 st rec} $l {
		switch $st {
		    s {$d create line $x1 $yp $x2 $yp -capstyle round \
			   -activewidth 4 -activefill white \
			   -width 2 -fill black -tags rec_$rec}
		    f {$d create line $x1 $yp $x2 $yp -capstyle round \
			   -activewidth 4 -activefill white \
			   -width 2 -fill blue -tags rec_$rec}
		    r {$d create line $x1 $yp $x2 $yp -capstyle round \
			   -activewidth 4 -activefill white \
			   -width 2 -fill red -tags rec_$rec}
		    t {$d create line $x1 $yp $x2 $yp -capstyle round \
			   -activewidth 4 -activefill white \
			   -width 1 -fill grey60 -tags rec_$rec}
		    d {$d create line $x1 $yp $x2 $yp -capstyle round \
			   -activewidth 4 -activefill white \
			   -width 1 -fill grey60 -dash . -tags rec_$rec}
		    D {$d create line $x1 $yp $x2 $yp -capstyle round \
			   -activewidth 4 -activefill white \
			   -width 1 -fill orange -dash . -tags rec_$rec}
		}
	    }
	    incr y
	}
	incr ystart $ymax
    }

    # Sync Y scrollbar
    set hei [expr {($ymax+$ystart)*$ysep+3.0}]
    set ${t}(scroll_height) $hei
    puts "[set ${t}(ys)] set [expr {$y1/$hei}] [expr {$y2/$hei}]"
    [set ${t}(ys)] set [expr {$y1/$hei}] [expr {$y2/$hei}]

    puts SEQ_SEQS_DRAW_END
}


# Mouse-over event callback
proc seq_seqs_bind {t w canvas} {
    global $w $t

    # Find the current reading record number
    set tags [$canvas gettags current]
    set rec [lindex [regexp -inline {(^| )rec_([0-9]+)} $tags] 2]
    if {$rec == ""} return

    # Produce some information about it
    set io [set ${w}(io)]
    set r1 [$io get_sequence $rec]
    set c1 [$r1 get_contig]
    set p1 [$r1 get_position]

    set info "[$r1 get_name] length [$r1 get_len]"

    set pair [$r1 get_pair]
    if {$pair} {
	set r2 [$io get_sequence $pair]
	set c2 [$r2 get_contig]
	set p2 [$r2 get_position]
	set l2 [$r2 get_length]
	append info ", pair [$r2 get_name] length $l2"
	if {$c1 != $c2} {
	    append info " (contig\#$c2)"
	} else {
	    set size [expr {$p2+abs($l2)-$p1+1}]
	    append info " (insert size $size)"
	}
    } else {
	append info " single ended (rec \#$rec)"
    }

    set ${w}(info) $info
}

proc lsort_callback {a b} {
    set base 4
    set l1 [expr {int(log([lindex $a end-2]-[lindex $a 0])/log($base))}]
    set l2 [expr {int(log([lindex $b end-2]-[lindex $b 0])/log($base))}]

    if {$l2 > $l1} {
	return 1;
    } elseif {$l1 > $l2} {
	return -1;
    }

    # Same approximate length, so sort by position still
    set p1 [lindex $a 0]
    set p2 [lindex $b 0]
    
    return [expr {$p2>$p1 ? -1 : ($p2<$p1 ? 1 : 0)}]
}

#
# Plots the sequence ruler
#
proc seq_ruler {w t x1 x2 y1 y2} {
    global $w $t

    set c    [set ${w}(contig)]
    set d    [set ${t}(canvas)]
    set wid  [set ${w}(width)]

    set wx1 [x2c $w $x1]
    set wx2 [x2c $w $x2]

    set h  [expr {$y2-$y1}]
    set h1 [expr {$y1+0.9*($y2-$y1)}]
    set h2 [expr {$y1+0.1*($y2-$y1)}]
    set h3 [expr {$y1+0.7*($y2-$y1)}]
    set h4 [expr {$y1+0.3*($y2-$y1)}]
    set h5 [expr {$y1+0.5*($y2-$y1)}]

    $d create rectangle $wx1 0 $wx2 $h -fill \#ffcccc

    $d create line $wx1 $h3 $wx2 $h3 -fill black

    foreach {x1 step p1 pstep fmt} [nice_num $x1 [expr {$x2-$x1}] 0] break
    set dist [expr {$step/[set ${w}(xzoom)]}]
    for {set i $x1; set j $p1} {$i < $x2} {set i [expr {$i+$step}]; set j [expr {$j+$pstep}]} {
	$d create line [x2c $w $i] $h2 [x2c $w $i] $h3
	if {$dist < 200} {
	    $d create text [x2c $w $i] $h3 -anchor nw \
		-text [format $fmt $j]
	}
    }

    set step2 [expr {$step/2.0}]
    set pstep2 [expr {$pstep/2.0}]
    for {set i $x1; set j $p1} {$i < $x2} {set i [expr {$i+$step2}]; set j [expr {$j+$pstep2}]} {
	$d create line [x2c $w $i] $h4 [x2c $w $i] $h3
	if {$dist >= 200} {
	    $d create text [x2c $w $i] $h3 -anchor nw \
		-text [format $fmt $j]
	}
    }

    set step3 [expr {$step/10.0}]
    for {set i $x1} {$i < $x2} {set i [expr {$i+$step3}]} {
	$d create line [x2c $w $i] $h5 [x2c $w $i] $h3
    }
}

# Picks a nice number for val. Ie value to only a few decimal places
proc nice_num {v r p} {
    set v [expr {int($v)}]
    set e [expr {int(log($r)/log(10)+1.0e-10)-$p}]
    set i [expr {int(pow(10,$e))}]
    set n [expr {($v/$i)*$i}]

    # Construct sensible ways of printing the numbers too
    set n2 $n
    set i2 $i
    set level 0
    while {$i2 >= 1000} {
	incr level
	set i2 [expr {$i2/1000.0}]
	set n2 [expr {$n2/1000.0}]
    }
    set suffix [lindex {{} K M G T} $level]
    set len [string length [regsub {(.*\.(0$)?)|(^[^.]*$)} [expr {$i2/2}] {}]]
    set format "%-.${len}f$suffix"

    return [list $n $i $n2 $i2 $format]
}


##############################################################################
#user interface dialog box for reading coverage histogram
proc ReadingCoverage { io } {
    global gap5_defs

    set f [keylget gap5_defs READING_COVERAGE.WIN]
    if {[xtoplevel $f -resizable 0] == ""} return
    wm title $f "reading coverage"

    contig_id $f.id \
	    -io $io \
	    -range 0
    
    #OK and Cancel buttons
    okcancelhelp $f.ok_cancel \
	    -ok_command "ReadingCoverage2 $io $f $f.id" \
	    -cancel_command "destroy $f" \
	    -help_command "show_help gap4 {Consistency-ReadingCov}" \
	    -bd 2 \
	    -relief groove

    pack $f.id $f.ok_cancel -side top -fill both
}

##############################################################################
#stand alone quality display
proc ReadingCoverage2 { io f id} {
    
    if {[set contign [contig_id_gel $id]] == ""} {bell; return}
    SetContigGlobals $io $contign

    # stop windows from hiding the plot
    destroy $f


    set pwin .read_depth
    set cname [lindex $contign 0]
    1.5Dplot $pwin $io [db_info get_contig_num $io $cname]
    add_plot $pwin seq_depth 50  -bd 2 -relief raised
    add_plot $pwin seq_seqs -200 -bd 2 -relief raised
    add_plot $pwin seq_ruler 50  -bd 2 -relief raised
}

# Test version when running as a script in its own right
if {[string match "*depth.tcl" $argv0]} {
    source $env(STADTABL)/shlib.conf
    load $env(STADLIB)/$env(MACHINE)-binaries/${lib_prefix}tk_utils${lib_suffix}
    load_package tk_utils
    tk_utils_init
    load_package gap5
    load libtgap.so g5

    #set dbname "/lustre/seq/scratch1/jkb/pe_rmdup.0"
    #set dbname "/lustre/seq/scratch1/jkb/1112_s_1"
    foreach dbname $argv {break;}
    set io [g5::open_database -name $dbname -access rw]

    wm withdraw .
    bind . <1> {puts "[llength [info command]] commands"}

    #add_plot $w bg_grid
    set pwin .plot
    1.5Dplot $pwin $io
    add_plot $pwin seq_depth  50 -bd 2 -relief raised
    add_plot $pwin seq_seqs -200 -bd 2 -relief raised
    add_plot $pwin seq_ruler  50 -bd 2 -relief raised
    #add_plot $pwin bg_grid
}