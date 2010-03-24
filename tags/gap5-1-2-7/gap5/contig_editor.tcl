#
# The contig editor consists of multiple components, following the standard
# model, view, controller pattern.
#
# 1) Controller: A top-level window, w, containing one or two editors plus
#    all the usual GUI buttons and menus.
#    Display settings are local to this and are held within the global tcl
#    array named after $w. In the code I typically upvar this to a
#    local array named opt().
#    Button/checkbutton GUI textvariables are stored within opt
#    starting with capital letters.
#        opt(PackSequences)
#        opt(HideAnno)
#
# 2) Model: The contig itself (with an associated io). The gapio ($io) is
#    typically a child io so we can edit (copy on write).
#    The contig has a global array named contigIO_$crec, usually
#    upvared to cio. This is the level at which we handle undo and
#    contig registration events.
#        cio(Undo) = <commands>
#        cio(Redo) = <commands>
#        cio(io)   = io=0x69ec60
#        cio(base) = io=0x68afa0
#        cio(crec) = 298459
#        cio(ref)  = 2
#
# 3) View: A names and editor tk widget (ednames and editor), attached to the
#    io/contig. This will be something like $w.seqs.sheet.
#    It doesn't have much Tcl data locally as the settings are mainly
#    held in the C widgets.
#        .e1.ed1.pane.seq.sheet(displayPos) = 658122
#        .e1.ed1.pane.seq.sheet(parent)     = .e1.ed1.pane
#        .e1.ed1.pane.seq.sheet(reg)        = 2
#        .e1.ed1.pane.seq.sheet(top)        = .e1

# The join editor is a stack of two contig editors. As such we have
# toplevel window with settings that govern both visible editors.
# $opt(all_editors) is a list of editors visible for window.

#-----------------------------------------------------------------------------
# IO/contig specific components. (child IOs, undo history)
#-----------------------------------------------------------------------------

# Conceptually each edit had an equal and opposite edit, although it may
# require some compound edits to achieve these.
#
# When we make a change we store on the undo list the command(s) required
# to reverse that change. (We actually just store an opcode & operands as
# this is more space efficient - see io_undo_exec.) 
#
# Note that in order for undo to restore the editing cursor position we
# typically store a cursor move event plus the edit event. These are grouped
# as a tcl list so they can be undone with as a single event. Cursor
# positioning is not a data-change, rather a view change. Hence the editor
# that has the cursor moved is the editor the user clicked Undo in (while
# the edit itself may have been made in another editor window on the same
# contig).
#
# Redo, for now, has been disabled until we get a realiable undo
# implementation fully tested.

proc io_child {io crec} {
    upvar \#0 contigIO_$crec cio

    if {![info exists cio(io)]} {
	set cio(base) $io
	set child     [$io child]
	set cio(io)   $child
	set cio(crec) $crec
	set cio(ref)  0
    }

    incr cio(ref)

    return $cio(io)
}

proc io_detach {crec} {
    upvar \#0 contigIO_$crec cio

    if {![info exists cio]} return

    incr cio(ref) -1

    if {[set cio(ref)] == 0} {
	$cio(io) close
	unset cio
    }
}

proc io_undo_exec {w crec cmdu} {
    upvar \#0 contigIO_$crec cio
    global $w

    set io [$w io]

    foreach cmd $cmdu {
	foreach {code op1 op2 op3 op4} $cmd break
	switch -- $code {
	    C_SET {
		$w set_cursor $op1 $op2 $op3 1
		# This may change the Cutoffs status, so update button
		set top [set ${w}(top)]
		global $top
		set ${top}(Cutoffs) [lindex [$w configure -display_cutoffs] 4]
	    }

	    B_REP {
		set seq [$io get_sequence $op1]
		$seq replace_base $op2 $op3 $op4
		$seq delete
	    }

	    B_INS {
		set seq [$io get_sequence $op1]
		$seq insert_base $op2 $op3 $op4
		$seq delete
	    }

	    B_DEL {
		set seq [$io get_sequence $op1]
		$seq delete_base $op2
		$seq delete
	    }

	    B_MOVE {
		$w decr_contig
		set c [$io get_contig $cio(crec)]
		foreach {p f} [$c remove_sequence $op1] break;
		$c add_sequence $op1 $op2 $p $f
		$c delete
		$w incr_contig

		eval $w set_cursor [$w get_cursor relative]
	    }

	    B_CUT {
		set seq [$io get_sequence $op1]
		$seq set_clips $op2 $op3
		$seq delete
	    }

	    C_INS {
		set contig [$io get_contig $op1]
		$contig insert_base $op2
		foreach seq $op3 {
		    foreach {rec pos base val cut} $seq break;
		    set seq [$io get_sequence $rec]
		    $seq replace_base $pos $base $val
		    foreach {b q c} [$seq get_base $pos] break
		    if {$c != $cut} {
			# Base was deleted from cutoff, but now in
			# used portion. Adjust clips to compensate
			set orient [$seq get_orient]
			set len [$seq get_length]
			foreach {l r} [$seq get_clips] break;
			if {$orient != 0} {
			    set pos [expr {abs($len)-$pos-1}]
			}

			if {$pos == [expr {$l-1}]} {
			    incr l
			} elseif {$pos == [expr {$r-1}]} {
			    incr r -1
			}

			$seq set_clips $l $r
		    }
		    $seq delete
		}
	    }
	    
	    C_DEL {
		set contig [$io get_contig $op1]
		$contig delete_base $op2
		$contig delete
	    }

	    T_DEL {
		set tag [$io get_anno_ele $op1]
		$tag remove
	    }

	    T_NEW {
		array set d $op1
		set rec [$io new_anno_ele $d(otype) $d(orec) $d(start) $d(end)]
		set t [$io get_anno_ele $rec]
		$t set_comment $d(anno)
		$t set_type $d(type)
		$t delete
	    }

	    T_MOVE {
		set s [$io get_sequence $op1]
		$s move_annos $op2
		$s delete
	    }

	    T_MOD {
		array set d $op2
		set tag [$io get_anno_ele $op1]
		if {[$tag get_comment] != $d(anno)} {
		    $tag set_comment $d(anno)
		}
		if {[$tag get_type] != $d(type)} {
		    $tag set_type $d(type)
		}
		$tag delete
	    }
	    
	    default {
		puts stderr "Unknown undo command: $cmd"
	    }
	}
    }
}

proc io_store_undo {crec cmdu cmdr} {
    upvar \#0 contigIO_$crec cio

    lappend cio(Undo) [list $cmdu $cmdr]
    set cio(Redo) ""

    contig_notify -io $cio(base) -cnum $crec -type GENERIC \
	-args [list TASK_GENERIC "" data {undo normal redo disable}]
    contig_notify -io $cio(base) -cnum $crec -type CHILD_EDIT -args ""
}

proc io_undo {ed crec} {
    upvar \#0 contigIO_$crec cio
    
    foreach {cmdu cmdr} [lindex [set cio(Undo)] end] break
    io_undo_exec $ed $crec $cmdu
    #eval $cmdu
    #lappend cio(Redo) [list $cmdu $cmdr]
    set cio(Undo) [lrange $cio(Undo) 0 end-1]

    if {[llength $cio(Undo)] == 0} {
	contig_notify -io $cio(base) -cnum $crec -type GENERIC \
	    -args [list TASK_GENERIC "" data {undo disable}]
    }
    contig_notify -io $cio(base) -cnum $crec -type GENERIC \
	-args [list TASK_GENERIC "" data {redo normal}]
    contig_notify -io $cio(base) -cnum $crec -type CHILD_EDIT -args ""
}

proc io_undo_state {crec} {
    upvar \#0 contigIO_$crec cio
    if {[info exists cio(Undo)] && [llength $cio(Undo)] != 0} {
	return normal
    }
    return disabled
}

proc io_redo {crec} {
    upvar \#0 contigIO_$crec cio
    
    foreach {cmdu cmdr} [lindex [set cio(Redo)] end] break
    eval $cmdr
    lappend cio(Undo) [list $cmdu $cmdr]
    set cio(Redo) [lrange $cio(Redo) 0 end-1]

    if {[llength $cio(Redo)] == 0} {
	contig_notify -io $cio(base) -cnum $crec -type GENERIC \
	    -args [list TASK_GENERIC "" data {redo disable}]
    }
    contig_notify -io $cio(base) -cnum $crec -type GENERIC \
	-args [list TASK_GENERIC "" data {undo normal}]
    contig_notify -io $cio(base) -cnum $crec -type CHILD_EDIT -args ""
}

proc io_redo_state {crec} {
    upvar \#0 contigIO_$crec cio
    if {[info exists cio(Redo)] && [llength $cio(Redo)] != 0} {
	return normal
    }
    return disabled
}

#-----------------------------------------------------------------------------
# Contig registration hookups. This data is obviously held per view rather
# than per contig or per io.
#-----------------------------------------------------------------------------
proc contig_register_callback {ed type id args} {
    global $ed
    set w [set ${ed}(top)]
    global $w

    #puts [info level [info level]]

    switch $type {
	QUERY_NAME {
	    return "Contig Editor"
	}

	CHILD_EDIT -
	LENGTH {
	    $ed redraw
	}

	GENERIC {
	    if {$ed != [set ${w}(curr_editor)]} return
	    foreach {component state} [lindex $args 2] {
		switch $component {
		    "undo" {
			$w.toolbar.undo configure -state $state
		    }
		    "redo" {
			#$w.toolbar.redo configure -state $state
		    }
		}
	    }
	}
	
	CURSOR_NOTIFY {
	    foreach a $args {
		foreach {k v} $a break;
		set arg($k) $v
	    }

	    # Only move the cursor if it's not sent by oursleves and
	    # it's our primary cursor.
	    # EDIT: Removed " || $arg(id) == 0"
	    if {[set ${ed}(reg)] != $arg(sent_by) && \
		    ($arg(id) == [$ed cursor_id])} {
		$ed set_cursor 17 [set ${w}(contig)] $arg(abspos)
	    }
	}

	JOIN_TO {
	    foreach a $args {
		foreach {k v} $a break;
		set arg($k) $v
	    }

	    # NB: What happens when we have contig 1 and contig 2 open, plus
	    # a join editor on 1+2 which we then join. We end up with two
	    # $io objects both for the same contig.
	    #
	    # We can ignore this whole problem as we can only join
	    # after saving.
	    #io_detach $arg(contig_num)
	    set ${w}(io) [io_child [set ${w}(io_base)] $arg(contig)]
	    $ed io [set ${w}(io)]

	    # Update the editor cached contig record details
	    set ${w}(-contig) $arg(contig)
	    $ed incr_contig $arg(contig)
	    upvar \#0 contigIO_$arg(contig) cio
	    set cio(Undo) ""
	    set cio(Redo) ""
	    set cio(io) [$ed io]
	    set cio(crec) $arg(contig)
	    incr cio(ref)
	    
	    # If cursor was on the old contig record, move it to the
	    # new contig instead.
	    foreach {type rec pos} [$ed get_cursor relative] break
	    if {$rec == $arg(contig_num)} {
		set rec $arg(contig)
		incr pos $arg(offset)
		$ed set_cursor $type $rec $pos 0
	    }

	    # Finally adjust the display.
	    set pos [$ed xview]
	    incr pos $arg(offset)
	    $ed xview $pos
	}

	GET_LOCK {
	    foreach a $args {
		foreach {k v} $a break;
		set arg($k) $v
	    }

	    # Check for lock 2 => WRITE
	    if {[expr {$arg(lock)&2}] == 2} {
		# Attempt to exit => save dialogue box
		if {[info exists ${w}(-contig2)] || ![editor_exit $w 1]} {
		    return [expr {$arg(lock) & ~2}]  ;# Disallow lock
		}
		update idletask                      ;# Ensure it's flushed
	    }

	    return $arg(lock)                        ;# Permit lock
	}

	REGISTER -
	DEREGISTER {
	    # nothing to do
	}
	
	default {
	    puts "Event '$type $id $args' not handled"
	}
    }

    return
}

#-----------------------------------------------------------------------------
# User dialogue for starting up the editors
#-----------------------------------------------------------------------------
proc EditContig2 {io t id} {
    if {[set reading [contig_id_gel $id]] == ""} return
    if {[set crec [contig_id_rec $id]] == ""} return

    destroy $t
    catch {edit_contig -io $io -contig $crec -reading $reading}
    SetContigGlobals $io $reading
}

proc EditContig {io} {
    set t .cedialog
    if {[xtoplevel $t -resizable 0] == ""} { return }
    wm title $t "Edit contig"

    contig_id $t.id -io $io -range 0 -command "EditContig2 $io $t $t.id"
    okcancelhelp $t.but \
	-ok_command "EditContig2 $io $t $t.id" \
	-cancel_command "destroy $t" \
	-help_command "show_help gap4 {Editor}"

    pack $t.id -side top -fill both
    pack $t.but -side bottom -fill both
}

proc JoinContig2 {io t id1 id2} {
    if {[set read1 [contig_id_gel $id1]] == ""} return
    if {[set read2 [contig_id_gel $id2]] == ""} return
    if {[set crec1 [contig_id_rec $id1]] == ""} return
    if {[set crec2 [contig_id_rec $id2]] == ""} return

    destroy $t
    join_contig -io $io \
	-contig  $crec1 -reading  $read1 -pos  1 \
	-contig2 $crec2 -reading2 $read2 -pos2 1
    SetContigGlobals $io $read1
}

proc JoinContig {io} {
    set t .jedialog
    if {[xtoplevel $t -resizable 0] == ""} { return }
    wm title $t "Join contigs"

    contig_id $t.id1 -io $io -range 0 -default "" -trace 2
    contig_id $t.id2 -io $io -range 0 -default "" -trace 0
    okcancelhelp $t.but \
	-ok_command "JoinContig2 $io $t $t.id1 $t.id2" \
	-cancel_command "destroy $t" \
	-help_command "show_help gap4 {Editor-Joining}"

    #initialise current frame
    SetCurFrame $t [list $t.id1 $t.id2]
    bind [entrybox_path $t.id1.ent] <<select>> "SetCurFrame $t {$t.id1 $t.id2}"
    bind [entrybox_path $t.id2.ent] <<select>> "SetCurFrame $t {$t.id2 $t.id1}"
    pack $t.id1 $t.id2 -side top -fill both
    pack $t.but -side bottom -fill both
}

#-----------------------------------------------------------------------------
# Externally usable functions.
#-----------------------------------------------------------------------------

# The top-level interface called from TCL
proc edit_contig {args} {
    eval contig_editor [next_editor] $args
}

# The top-level interface called from TCL
proc join_contig {args} {
    eval contig_editor [next_editor] $args
}

#-----------------------------------------------------------------------------
# Internally usable functions.
#-----------------------------------------------------------------------------

# Allocates a window pathname for the contig editor
set editor_num 0
proc next_editor {} {
    global editor_num
    incr editor_num
    return ".e$editor_num"
}

#
# Creates a contig editor mega-widget. The layout top to bottom is:
# Menu bar
# Tool/button bar
# Paned window (left/right split)
# Status line
#
# The paned window contains the actual editor itself with the left side being
# Name + scrollbar underneath and the right side being editor + 2 scrollbars
# below and right.
#
# Returns the pathname of the newly created object: '$path'
#
proc contig_editor {w args} {
    global gap5_defs
    upvar \#0 $w opt

    # Initialise the $path global array - the instance data for this widget
    foreach {arg val} $args {
	set opt($arg) $val
    }

    set opt(win) $w
    set opt(Disagreements) 0
    set opt(DisagreeMode)  1
    set opt(DisagreeCase)  1
    set opt(PackSequences) 1
    set opt(HideAnno)      0
    set opt(Status)        "--- Status info here ---"

    set opt(io_base) $opt(-io)
    set opt(io) [io_child $opt(-io) $opt(-contig)]

    set join [info exists opt(-contig2)]

    #set opt(contig) [contig_order_to_number -io $opt(-io) -order 0]
    set opt(contig) $opt(-contig)
    if {[info exists opt(-reading)]} {
	set opt(-reading) [$opt(io) seq_name2rec $opt(-reading)]
	if {$opt(-reading) == -1} { set opt(-reading) 0 }
    } else {
	set opt(-reading) 0
    }
    if {![info exists opt(-pos)]} { set opt(-pos) 0 }
    if {$join} {
	set opt(contig2) $opt(-contig2)
	set opt(io2) [io_child $opt(-io) $opt(-contig2)]
	if {[info exists opt(-reading2)]} {
	    set opt(-reading2) [$opt(io) seq_name2rec $opt(-reading2)]
	    if {$opt(-reading2) == -1} { set opt(-reading2) 0 }
	} else {
	    set opt(-reading2) 0
	}
	if {![info exists opt(-pos2)]}     { set opt(-pos2) 0 }
    }

    # Create the window layout
    if {![winfo exists $w]} {
	toplevel $w
	set c [$opt(io) get_contig $opt(contig)]
	#$c dump_ps /tmp/tree.ps
	if {[$opt(io) read_only]} {
	    set extra "   *** READ-ONLY ***"
	} else {
	    set extra ""
	}
	if {$join} {
	    set c2 [$opt(io2) get_contig $opt(contig2)]
	    wm title $w "Join: [$c get_name] / [$c2 get_name]$extra"
	} else {
	    wm title $w "Edit: [$c get_name]$extra"
	}
    }
    wm resizable $w 1 1

    # The toolbar 
    set tool [frame $w.toolbar -bd 0]
    checkbutton $tool.cutoffs \
	-variable ${w}(Cutoffs) \
	-text Cutoffs \
	-command "editor_cutoffs $w"
    checkbutton $tool.quality \
	-variable ${w}(Quality) \
	-text Quality \
	-command "editor_quality $w"
    button $tool.undo    -text Undo -command "editor_undo $w" \
	-state [io_undo_state $opt(contig)]
#    button $tool.redo    -text Redo -command "editor_redo $w" \
	-state [io_redo_state $opt(contig)]
    button $tool.search  -text Search \
	-command "create_search_win $w.search \"editor_search $w\" 0"
    button $tool.save -text Save -command "editor_save $w"
    wm protocol $w WM_DELETE_WINDOW "editor_exit $w"
    pack $tool.undo $tool.search $tool.cutoffs $tool.quality \
	-side left
    pack $tool.save -side right

    # Highlights of the current editor so we know what window the button
    # applies to.
    bind $tool.undo <Any-Enter> "editor_hl \[set ${w}(curr_editor)\] red"
    bind $tool.undo <Any-Leave> "editor_hl \[set ${w}(curr_editor)\] \#d9d9d9"
#    bind $tool.redo <Any-Enter> "editor_hl \[set ${w}(curr_editor)\] red"
#    bind $tool.redo <Any-Leave> "editor_hl \[set ${w}(curr_editor)\] \#d9d9d9"

    if {$join} {
	set opt(Lock) 1; # See default in tkEditor.c link_to command
	checkbutton $tool.lock \
	    -variable ${w}(Lock) \
	    -text Lock \
	    -command "editor_lock $w"
	pack $tool.lock -side left

	button $tool.align \
	    -text Align -padx 2 \
	    -command "\[set ${w}(curr_editor)\] join_align"
	button $tool.alignL \
	    -text "<" -padx 2 \
	    -command "\[set ${w}(curr_editor)\] join_align 0 1"
	button $tool.alignR \
	    -text ">" -padx 2 \
	    -command "\[set ${w}(curr_editor)\] join_align 1 0"
	pack $tool.alignL $tool.align $tool.alignR \
	    -side left -fill both -padx 0

	button $tool.join -text Join -command "editor_join $w"
	pack $tool.join -side right

	if {[$c get_rec] == [$c2 get_rec]} {
	    $tool.join configure -state disabled
	}
    }

    if {[$opt(io) read_only]} {
	$tool.save configure -state disabled
	catch {$tool.join configure -state disabled}
    }

    # The editor(s) itself
    if {$join} {
	set pane0 $w.ed0
	set e [editor_pane $w $pane0 1 2 opt]
	set opt(editor2) $e
	lappend opt(all_editors) $e
    }
    set pane1 $w.ed1
    set e [editor_pane $w $pane1 0 "" opt]
    set opt(editor1) $e
    lappend opt(all_editors) $e

    # Difference bar for the join editor
    if {$join} {
	set diffs [diff_pane $w.diffs]
	$e link_to $opt(editor2) $diffs.pane.seq.sheet
    }

    # The bottom status line
    set status $w.status
    frame $status -bd 2 -relief groove
    label $status.dummy
    label $status.l -textvariable ${w}(Status)
    pack  $status.dummy -fill both
    place $status.l -relx 0

    # Menu bar
    global contig_editor_main_menu
    $w configure -menu $w.menubar
    menu $w.menubar 
    create_menus $contig_editor_main_menu $w.menubar
    bind  $w.menubar <ButtonPress> "tag_repopulate_menu \[set ${w}(curr_editor)\]"

    # Packing
    grid rowconfigure $w 3 -weight 1
    grid columnconfigure $w 0 -weight 1
    grid $tool   -sticky nsew -row 0
    if {$join} {
	grid rowconfigure $w 1 -weight 1
	grid $pane0  -sticky nsew -row 1
	grid $diffs  -sticky nsew -row 2
    }
    grid $pane1  -sticky nsew -row 3
    grid $status -sticky nsew -row 4

    # Synchronised pane movement
    set opt(panes) $pane1.pane
    if {$join} {
	lappend opt(panes) $diffs.pane $pane0.pane
    }

    foreach p $opt(panes) {
	bind $p <ButtonRelease-1> "+sync_panes %W $w 1 1"
	bind $p <ButtonRelease-2> "+sync_panes %W $w 0 1"

	bind $p <Any-B2-Motion> "+sync_panes %W $w 0 0"
    }
    sync_panes $pane1.pane opt 0 1

    # Grid control
    # In theory this should work, but in practice getting the panedwindow
    # widget to correctly honour the side of the internal components is
    # a massive exercise in frustration, let alone making the pane "sash"
    # only move in increments of one font element.

    update idletasks
    array set font [font metrics sheet_font]
#    wm grid $w [winfo width .e1] [winfo height .e1] \
#	[font measure sheet_font A] $font(-linespace)


    $e redraw
}

proc editor_search {w args} {
    global $w
    set ed [set ${w}(curr_editor)]
    eval $ed search $args
}

proc editor_lock {w} {
    global $w
    set ed [set ${w}(curr_editor)]
    $ed lock [set ${w}(Lock)]
}

proc editor_save {w} {
    upvar \#0 $w opt

    foreach ed $opt(all_editors) {
	if {[$ed save] != 0} {
	    bell
	}
    }
}

proc editor_join {w} {
    upvar \#0 $w opt

    foreach ed $opt(all_editors) {
	if {[$ed save] != 0} {
	    bell
	    return
	}
    }

    $ed join
    editor_exit $w
}

# Returns true if we really want to exit
#         false if we cannot exit (user hit cancel, or failed to save).
proc editor_exit {w {get_lock 0}} {
    global $w
    set ed [set ${w}(curr_editor)]

    if {[winfo exists $w.save_dialog]} return

    # Two styles of exit dialog depending on whether we arrived here with
    # $get_lock true, indicating this wasn't a user controlled exit but rather
    # a request originating in another window due to the requirement of
    # taking write-access to this contig.
    if {![[set ${w}(io)] read_only] && [$ed edits_made]} {
	if {$get_lock} {
	    set ret [tk_dialog \
			 $w.save_dialog \
			 "Quit editor?" \
			 "Another window wishes to modify this contig, shutting down the editor in the process.\n\nTo deny this hit Cancel.\n\nOtherwise the editor will exit. You can choose whether to save changes when this happens." \
			 "" \
			 2 \
			 Save {Don't Save} Cancel]

	    set ret [lindex {yes no cancel} $ret]
	} else {
	    set ret [tk_messageBox \
			 -icon question \
			 -title "Save changes" \
			 -message "Edits have been made. Save changes?" \
			 -default yes \
			 -type yesnocancel \
			 -parent $w]
	}
	
	if {$ret == "cancel"} {
	    return 0
	} elseif {$ret == "yes"} {
	    if {[$ed save] != 0} {
		bell
		return 0
	    }
	}
    }

    set detach ""
    foreach ed [set ${w}(all_editors)] {
	global $ed
	set id [set ${ed}(reg)]
	contig_deregister -io [set ${w}(io_base)] -id $id
	lappend detach [$ed contig_rec]
    }

    destroy $w

    foreach crec $detach {
	io_detach $crec
    }

    return 1
}

proc display_pos_set {ed pos} {
    global $ed
    set ${ed}(displayPos) $pos
}

proc jog_editor {cmd dist} {
    $cmd xview scroll $dist units
}

# Lays out the editor pane, with options of the scrollbars being above or
# below the editor (used in the join editor).
# Above indicates whether the scrollbars are above the text panels or below.
proc editor_pane {top w above ind arg_array} {
    upvar $arg_array opt
    global gap5_defs

    if {$above != 0} {
	set above 1
	set jogrow 0
	set scrollrow 1
	set textrow 2
	set cattop 0
    } else {
	set textrow   0
	set scrollrow 1
	set jogrow 2
	set cattop 1
    }

    set f $w
    frame $f -bd 3

    set w $f.pane
    panedwindow $w -orient horiz -bd 1 -relief sunken \
	-showhandle 1 -sashrelief raised

    frame $w.name -bd 0 -highlightthickness 0
    frame $w.seq -bd 0 -highlightthickness 0

    $w add $w.name $w.seq

    # Seqs panel
    set ed $w.seq.sheet
    set ed [editor $ed \
		-width 80 \
		-height 16 \
		-xscrollcommand "display_pos_set $ed \[$ed xview\];
                                 $w.seq.x set" \
		-yscrollcommand "$w.seq.y set" \
		-qual_fg     [keylget gap5_defs CONTIG_EDITOR.QUAL_IGNORE] \
		-diff1_bg    [keylget gap5_defs CONTIG_EDITOR.DIFF1_BG] \
		-diff2_bg    [keylget gap5_defs CONTIG_EDITOR.DIFF2_BG] \
		-diff1_fg    [keylget gap5_defs CONTIG_EDITOR.DIFF1_FG] \
		-diff2_fg    [keylget gap5_defs CONTIG_EDITOR.DIFF2_FG] \
		-stripe_bg   [keylget gap5_defs CONTIG_EDITOR.STRIPE_BG] \
		-qualcolour0 [keylget gap5_defs CONTIG_EDITOR.QUAL0_COLOUR] \
		-qualcolour1 [keylget gap5_defs CONTIG_EDITOR.QUAL1_COLOUR] \
		-qualcolour2 [keylget gap5_defs CONTIG_EDITOR.QUAL2_COLOUR] \
		-qualcolour3 [keylget gap5_defs CONTIG_EDITOR.QUAL3_COLOUR] \
		-qualcolour4 [keylget gap5_defs CONTIG_EDITOR.QUAL4_COLOUR] \
		-qualcolour5 [keylget gap5_defs CONTIG_EDITOR.QUAL5_COLOUR] \
		-qualcolour6 [keylget gap5_defs CONTIG_EDITOR.QUAL6_COLOUR] \
		-qualcolour7 [keylget gap5_defs CONTIG_EDITOR.QUAL7_COLOUR] \
		-qualcolour8 [keylget gap5_defs CONTIG_EDITOR.QUAL8_COLOUR] \
		-qualcolour9 [keylget gap5_defs CONTIG_EDITOR.QUAL9_COLOUR] \
		-bd 0 \
	        -consensus_at_top $cattop \
	        -stack_mode  $opt(PackSequences) \
	        -hide_anno   $opt(HideAnno)]
    set opt(curr_editor) $ed

    # X and y scrollbars
    scrollbar $w.seq.x -orient horiz -repeatinterval 30
    scrollbar $w.seq.y -orient vert 

    # Names panel
    set edname $w.name.sheet
    ednames $edname \
	-width 15 \
	-height 16 \
	-xscrollcommand "$w.name.x set" \
	-bd 0

    scrollbar $w.name.x -orient horiz

    entry $w.name.pos \
	-textvariable ${ed}(displayPos) \
	-font {Helvetica -12} \
	-width 15
    bind $w.name.pos <Return> "editor_goto $ed $w.seq.sheet"

    # The jog control for scrolling the editor
    set posh [font metrics [$w.name.pos cget -font] -linespace]
    incr posh -2
    jog $w.seq.jog \
	-orient horiz \
	-command "jog_editor $w.seq.sheet" \
	-repeatinterval 50 -width $posh

    # Pack it all in the paned window
    grid rowconfigure $w.name $textrow -weight 1
    grid columnconfigure $w.name 0 -weight 1
    grid $w.name.sheet -row $textrow   -sticky nsew
    grid $w.name.x     -row $scrollrow -sticky nsew
    grid $w.name.pos   -row $jogrow    -sticky nsew

    $w.name configure -width [expr {11*[font measure sheet_font A]}]

    grid rowconfigure $w.seq $textrow -weight 1
    grid columnconfigure $w.seq 0 -weight 1
    grid $w.seq.sheet $w.seq.y -row $textrow   -sticky nsew
    grid $w.seq.x              -row $scrollrow -sticky nsew
    grid $w.seq.jog            -row $jogrow -sticky nsew

    focus $w.seq.sheet

    # Initialise with an IO and link name/seq panel together
    global $ed $edname
    set ${ed}(parent) $w
    set ${ed}(top) $top
    set ${ed}(reg) [contig_register \
			-io $opt(io_base) \
			-contig $opt(contig$ind) \
			-command "contig_register_callback $ed" \
			-flags [list ALL GENERIC CHILD_EDIT]]
    $ed init $opt(io$ind) $opt(contig$ind) $opt(-reading$ind) $opt(-pos$ind) $w.name.sheet

    if {$ind == 2} {
	set ${ed}(side) top
    } else {
	set ${ed}(side) bottom
    }
    set ${edname}(ed) $ed

    # Force new style mode
    $w.name.x set 0.0 0.1
    $w.seq.x set 0.0 0.1
    $w.seq.y set 0.0 0.1

    $w.seq.x configure -command "$w.seq.sheet xview"
    $w.seq.y configure -command "$w.seq.sheet yview"
    $w.name.x configure -command "$w.name.sheet xview"

    grid rowconfigure $f 0 -weight 1
    grid columnconfigure $f 1 -weight 1
    grid $f.pane           -row 0 -columnspan 2 -sticky nsew

    # Force scrollbar to be set to the correct size.
    $w.name.sheet xview moveto 0.0

    return $ed
}

# The "differences" bar that separates a pair of join editors
proc diff_pane {w} {
    frame $w -bd 0 -padx 3
    set p [panedwindow $w.pane -orient horiz -bd 2 -relief sunken \
	      -showhandle 1 -sashrelief raised]
    pack $p -fill both -expand 1

    frame $p.name -bd 0 -highlightthickness 0
    label $p.name.diff -text "Differences"
    pack $p.name.diff

    frame $p.seq -bd 0 -highlightthickness 0
    sheet $p.seq.sheet
    pack $p.seq.sheet -fill both -expand 1

    $p add $p.name $p.seq
    
    return $w
}

# Synchronises the pane position between a list of panes
proc sync_panes {w arg_array proxy round} {
    upvar $arg_array opt

    if {[winfo class $w] != "Panedwindow"} {
	set w $w.pane
    }

    if {$proxy} {
	foreach {px y} [$w proxy coord] {}
    } else {
	foreach {px y} [$w sash coord 0] {}
    }

    if {$round} {
	set fs [font measure sheet_font A]
	set x [expr {int($px / $fs)*$fs+[$w cget -bd]+[$w cget -sashpad]+3}]
    } else {
	set x $px
    }

    foreach pane $opt(panes) {
	if {[winfo class $pane] != "Panedwindow"} {
	    set pane $pane.pane
	}

	if {!$round && $pane == $w} continue

	after idle "$pane sash place 0 $x $y"
    }
}

# Highlights the curr_editor window
proc editor_hl {ed col} {
    set w [winfo parent [winfo parent [winfo parent $ed]]]
    $w configure -bg $col
}

# A coordinate jump
proc editor_goto {ed w} {
    upvar \#0 $ed eopt

    set pos $eopt(displayPos)

    if {[regexp {^[-+]?[0-9]+$} $pos] == 1} {
	$w xview $pos
    }
}

# Callback for cutoffs button
proc editor_cutoffs {w} {
    upvar \#0 $w opt

    foreach ed $opt(all_editors) {
	$ed configure -display_cutoffs $opt(Cutoffs)
	# Move cursor to the consensus if disabling the cutoffs has now
	# hidden it.
	if {$opt(Cutoffs) == 0} {
	    set curr [$ed get_cursor relative]
	    eval $ed set_cursor [$ed get_cursor absolute] 0
	    eval $ed set_cursor $curr 0
	}
	$ed redraw
    }
}

# Callback for quality button
proc editor_quality {w} {
    upvar \#0 $w opt

    foreach ed $opt(all_editors) {
	$ed configure -display_quality $opt(Quality)
	$ed configure -display_mapping_quality $opt(Quality)
	$ed redraw
    }
}

proc editor_disagreements {w} {
    upvar \#0 $w opt

    foreach ed $opt(all_editors) {
	if {$opt(Disagreements)} {
	    $ed configure -display_differences $opt(DisagreeMode)
	} else {
	    $ed configure -display_differences 0
	}
	$ed configure -differences_case_sensitive $opt(DisagreeCase)
	$ed redraw
    }
}

proc editor_toggle_annos {w} {
    global $w

    if {[info exists ${w}(all_editors)]} {
	# Called either via the menu, after setting ther new value
	upvar \#0 $w opt
    } else {
	# Or via control-Q on an editor, to toggle the value
	upvar \#0 [set ${w}(top)] opt
	set opt(HideAnno) [expr {1-$opt(HideAnno)}]
    }

    foreach ed $opt(all_editors) {
	$ed configure -hide_anno $opt(HideAnno)
	$ed redraw
    }
}

proc ed2name {w} {
    return [regsub {\.seq\.sheet} $w .name.sheet]
}

proc name2ed {w} {
    return [regsub {\.name\.sheet} $w .seq.sheet]
}

proc set_editor_pack_sequences {w} {
    upvar \#0 $w opt

    foreach ed $opt(all_editors) {
	$ed configure -stack_mode $opt(PackSequences)
	$ed xview scroll 1 units
	$ed xview scroll -1 units

	# Force redraw of name scrollbar
	[ed2name $ed] xview scroll 0 units
    }
}

proc set_differences_quality_callback {w val} {
    upvar \#0 $w opt

    foreach ed $opt(all_editors) {
	$ed configure -display_differences_quality $val
	$ed redraw
    }
}

proc set_differences_quality {w} {
    upvar \#0 $w opt
    global gap5_defs

    set ed $opt(curr_editor)

    set t $ed.qual_diff
    if {[xtoplevel $t -resizable 0] == ""} {return}
    wm title $t "Set differences quality"

    set start [lindex [$ed configure -display_differences_quality] 3]

    scalebox $t.qual \
	-title "Quality" \
	-orient horizontal \
	-from 0 \
	-to 100 \
	-width 5 \
	-default $start \
	-type CheckInt \
	-command "set_differences_quality_callback $w"
    $t.qual.scale configure -length 150

    okcancelhelp $t.ok \
	-ok_command "keylset gap5_defs CONTIG_EDITOR.DISAGREE_QUAL \[scalebox_get $t.qual\]; destroy $t" \
	-cancel_command "set_differences_quality_callback $w $start; destroy $t" \
        -help_command "show_help gap4 {Editor-Differences Quality}"

    pack $t.qual $t.ok -side top -fill both
}

#-----------------------------------------------------------------------------
# Undo support
proc store_undo {w cmdu cmdr} {
    io_store_undo [$w contig_rec] $cmdu $cmdr
}

proc editor_undo {top} {
    upvar \#0 $top opt
    
    set w $opt(curr_editor)
    io_undo $w [$w contig_rec]
}

# proc editor_redo {top} {
#     upvar \#0 $top opt
#     set w $opt(curr_editor)
#     io_redo [$w contig_rec]
# } 


#-----------------------------------------------------------------------------
# Basic sequence editing function
proc editor_edit_base {w call where} {
    upvar $w opt

    set io [$w io]
    upvar $opt(top) top
    
    if {$where == "" || [$io read_only]} {
	bell
	return
    }

    foreach {type rec pos} $where break;
    if {$type == 18} {
	set seq [$io get_sequence $rec]
	foreach {old_call old_conf} [$seq get_base $pos] break
	$seq replace_base $pos $call 30
	$seq delete

	foreach {type rec _pos} [$w get_cursor relative] break
	$w cursor_right
	store_undo $w \
	    [list \
		 [list C_SET $type $rec $pos] \
		 [list B_REP $rec $pos $old_call $old_conf] ] {}
    }

    $w redraw
}

proc editor_insert_gap {w where} {
    upvar $w opt

    set io [$w io]

    if {$where == "" || [$io read_only]} {
	bell
	return
    }

    foreach {type rec pos} $where break;
    if {$type == 18} {
	set seq [$io get_sequence $rec]
	$seq insert_base $pos * 20
	$seq delete

	store_undo $w \
	    [list \
		 [list C_SET $type $rec $pos] \
		 [list B_DEL $rec $pos] ] {}
    } else {
	set contig [$io get_contig $rec]
	$contig insert_base $pos
	$contig delete

	store_undo $w \
	    [list \
		 [list C_SET $type $rec $pos] \
		 [list C_DEL $rec $pos] ] {}
    }
    $w cursor_right
    
    $w redraw
}

proc editor_delete_base {w where} {
    upvar $w opt

    set io [$w io]

    if {$where == "" || [$io read_only]} {
	bell
	return
    }

    foreach {type rec pos} $where break;

    $w cursor_left
    incr pos -1

    if {$type == 18} {
	set seq [$io get_sequence $rec]
	foreach {old_base old_conf} [$seq get_base $pos] break
	foreach {l0 r0} [$seq get_clips] break;
	$seq delete_base $pos

	# Identify if we're in a situation where we need to undo the clip
	# points too. This occurs when, for example, we delete the last base
	# of the left-cutoff data. Undoing this would be an insertion to the
	# start of the used portion, making that base now visible instead.
	# We use a belt and braces method by temporarily adding a base back
	# to check our clip points are consistent
	$seq insert_base $pos A 0
	foreach {l1 r1} [$seq get_clips] break;
	$seq delete_base $pos

	$seq delete

	if {$l0 != $l1 || $r0 != $r1} {
	    store_undo $w \
		[list \
		     [list C_SET $type $rec [expr {$pos+1}]] \
		     [list B_INS $rec $pos $old_base $old_conf] \
		     [list B_CUT $rec $l0 $r0]] {}
	} else {
	    store_undo $w \
		[list \
		     [list C_SET $type $rec [expr {$pos+1}]] \
		     [list B_INS $rec $pos $old_base $old_conf]] {}
	}
    } else {
	set contig [$io get_contig $rec]

	# get pileup at consensus position so we can restore it
	set pileup [$contig get_pileup $pos]

	$contig delete_base $pos
	$contig delete

	store_undo $w \
	    [list \
		 [list C_SET $type $rec [expr {$pos+1}]] \
		 [list C_INS $rec $pos $pileup]] {}
    }
    
    $w redraw
}

proc editor_move_seq {w where direction} {
    upvar $w opt

    set io [$w io]

    if {$where == "" || [$io read_only]} {
	bell
	return
    }

    foreach {type rec _pos} $where break;
    if {$type != 18} {
	# sequences only
	bell
	return
    }

    set seq  [$io get_sequence $rec]
    set cnum [$seq get_contig]
    set pos  [$seq get_position]
    $seq delete

    store_undo $w \
	[list \
	     [list C_SET $type $rec $_pos] \
	     [list T_MOVE $rec [expr {-$direction}]] \
	     [list B_MOVE $rec $pos] ] {}

    $w decr_contig

    set c [$io get_contig $cnum]
    foreach {pair_rec flags} [$c remove_sequence $rec] break;
    incr pos $direction
    $c add_sequence $rec $pos $pair_rec $flags
    $c delete

    set seq [$io get_sequence $rec]
    $seq move_annos $direction
    $seq delete

    $w incr_contig

    # A bit obscure, but it ensures edSetApos() is called in C, keeping
    # cached absolute and relative positions in sync after the move.
    eval $w set_cursor [$w get_cursor relative]

    $w redraw
}

proc editor_clip_seq {w where end} {
    upvar $w opt

    set io [$w io]

    if {$where == "" || [$io read_only]} {
	bell
	return
    }

    foreach {type rec pos} $where break;
    if {$type != 18} {
	# sequences only
	bell
	return
    }

    set seq  [$io get_sequence $rec]
    set orient [$seq get_orient]

    foreach {l r} [$seq get_clips] break;
    store_undo $w \
	[list \
	     [list B_CUT $rec $l $r] \
	     [list C_SET $type $rec $pos]] {}

    switch $end {
	l {
	    if {$orient} {
		set r [expr {abs([$seq get_length])-$pos}]
	    } else {
		set l [expr {$pos+1}]
	    }
	}
	r {
	    if {$orient} {
		set l [expr {abs([$seq get_length])-$pos+1}]
	    } else {
		set r $pos
	    }
	}
    }

    if {$l <= $r} {
	$seq set_clips $l $r
    } else {
	bell
    }
    $seq delete

    $w redraw
}

# Updates the editor status line for editor $w.
# x and y are optional, but if set specify the location to display (eg
# for mouse motion events).
proc update_brief {w {name 0} {x {}} {y {}}} {
    global gap5_defs
    global $w

    if {$x != "" && $y != ""} {
	foreach {type rec pos} [$w get_number $x $y] break
    } else {
	foreach {type rec pos} [$w get_number] break
    }

    if {$name} {
	set w [name2ed $w]
	global $w
    }

    if {![info exists type]} {
	set w [set ${w}(top)]
	global $w
	set ${w}(Status) ""
	return
    }

    if {$name} {
	switch $type {
	    18 {
		set msg [$w get_seq_status $type $rec $pos \
			     [keylget gap5_defs READ_BRIEF_FORMAT]]
	    }
	    17 {
		set msg [$w get_seq_status $type $rec $pos \
			     [keylget gap5_defs CONTIG_BRIEF_FORMAT]]
	    }
	    21 {
		set msg "tag in name?"
	    }
	    default {
		set msg "Unknown data type $type"
	    }
	}
    } else {
	switch $type {
	    18 {
		set msg [$w get_seq_status $type $rec $pos \
			     [keylget gap5_defs BASE_BRIEF_FORMAT1]]
	    }
	    17 {
		set msg [$w get_seq_status $type $rec $pos \
			     [keylget gap5_defs BASE_BRIEF_FORMAT2]]
	    }
	    21 {
		set msg [$w get_seq_status $type $rec $pos \
			     [keylget gap5_defs TAG_BRIEF_FORMAT]]
		regsub -all "\n" $msg { \n } msg
	    }
	    default {
		set msg "Unknown data type $type"
	    }
	}
    }

    set w [set ${w}(top)]
    global $w
    set ${w}(Status) $msg
}

proc editor_name_select {w where} {
    global $w

    if {$where == ""} return

    foreach {type rec pos} $where break;
    if {$type != 18} return

    # Add name to XA_PRIMARY selection
    set seq [[$w io] get_sequence $rec]
    set name [$seq get_name]

    # FIXME: underline name too, via $ed call?
#    [set ${w}(ed)] ...

    # Ignore parameters for now. Assume reading name length <= maxbytes.
    catch {rename editor_selection_handler {}}
    proc editor_selection_handler {offset maxbytes} [list return $name]

    selection own $w
    selection handle $w editor_selection_handler

    # For Windows...
    clipboard clear
    clipboard append $name
}

#-----------------------------------------------------------------------------
# Tag editor windows

# Functions to make tag edits or to be called by the undo/redo stack.
proc U_tag_change {w rec new_a} {
    set io [$w io]

    #-- Get existing tag
    set old_a ""
    if {$rec != -1} {
	set tag [$io get_anno_ele $rec]
	set d(strand)  0; # fixme
	set d(type)    [$tag get_type]
	foreach {d(start) d(end)} [$tag get_position] break;
	set d(otype)   [$tag get_obj_type]
	set d(orec)    [$tag get_obj_rec]
	set d(anno)    [$tag get_comment]
	set d(default) "?"

	set old_a [array get d]
	unset d
    }

    #-- Create, Modify or Delete existing tag
    array set d $new_a

    if {$new_a == ""} {
	# Delete
	$tag remove

	store_undo $w \
	    [list \
		 [list T_NEW $old_a]] {}

#	    [list U_tag_change $w -1 $old_a] \
	    [list U_tag_change $w $rec ""]

    } elseif {$rec == -1} {
	# Create
	set rec [$io new_anno_ele $d(otype) $d(orec) $d(start) $d(end)]
	set t [$io get_anno_ele $rec]
	$t set_comment $d(anno)
	$t set_type $d(type)
	$t delete
	
	store_undo $w \
	    [list \
		 [list T_DEL $rec]] {}

    } else {
	# Modify
	if {[$tag get_comment] != $d(anno)} {
	    $tag set_comment $d(anno)
	}
	if {[$tag get_type] != $d(type)} {
	    $tag set_type $d(type)
	}
	$tag delete

	store_undo $w \
	    [list \
		 [list T_MOD $rec $old_a]] {}
    }

    unset d
}

proc tag_repopulate_menu {w} {
    foreach {rtype rrec rpos} [$w get_cursor relative] break
    foreach {atype arec apos} [$w get_cursor absolute] break

    upvar contigIO_$arec cio
    upvar $w wa

    set me $wa(top).menubar.commands.edit_tag
    $me delete 0 end
    $me configure -tearoff 0

    set md $wa(top).menubar.commands.delete_tag
    $md delete 0 end
    $md configure -tearoff 0

    # Perhaps not the most efficient approach, but it works
    set c [$cio(io) get_contig $cio(crec)]
    foreach anno [$c anno_in_range $apos $apos] {
	if {$rrec == [lindex $anno 8]} {
	    foreach {start end rec itype} $anno break
	    set type ""
	    while {$itype > 0} {
		set type "[format "%c" [expr {$itype & 0xff}]]$type"
		set itype [expr {$itype >> 8}]
	    }
	    $me add command -label "$type $start..$end \#$rec" \
		-command "tag_editor_launch $w $rec"
	    $md add command -label "$type $start..$end \#$rec" \
		-command "tag_editor_delete $w $rec"
	}
    }
    $c delete
}

proc tag_editor_launch {w where} {
    if {[llength $where] != 1} {
	foreach {type rec pos} $where break;
	if {$type != 21} return
    } else {
	set rec $where
    }

    set tag [[$w io] get_anno_ele $rec]
    global .Tag.$rec
    upvar \#0 .Tag.$rec d

    set d(strand)  0; # fixme
    set d(type)    [$tag get_type]
    set d(otype)   [$tag get_obj_type]
    set d(orec)    [$tag get_obj_rec]
    set d(anno)    [$tag get_comment]
    set d(default) "?"
    
    $tag delete

    create_tag_editor $w.tag_$rec "tag_editor_callback $w $rec" .Tag.$rec
}

proc tag_editor_delete {w where} {
    if {[llength $where] != 1} {
	foreach {type rec pos} $where break;
	if {$type != 21} return
    } else {
	set rec $where
    }

    U_tag_change $w $rec ""
#    set tag [[$w io] get_anno_ele $rec]
#    $tag remove

    $w redraw
}

proc tag_editor_create {w} {
    global $w
    set rec -1
    
    puts [$w select get]
    foreach {otype orec start end} [$w select get] break;

    global .Tag.$rec
    upvar \#0 .Tag.$rec d

    set d(strand)  0; # fixme
    set d(type)    "COMM"
    set d(otype)   $otype
    set d(orec)    $orec
    set d(start)   $start
    set d(end)     $end
    set d(anno)    "default"
    set d(default) "?"

    create_tag_editor $w.tag_$rec "tag_editor_callback $w $rec" .Tag.$rec
}

proc tag_editor_callback {w rec cmd args} {
    upvar \#0 .Tag.$rec d
    set f $w.tag_$rec
    set io [$w io]

    switch $cmd {
	"save" {
#	    if {$rec == -1} {
#		# Allocate a new item
#		set rec [$io new_anno_ele $d(otype) $d(orec) $d(start) $d(end)]
#		puts "New tag with rec $rec"
#		set t [$io get_anno_ele $rec]
#		$t set_comment $d(anno)
#		$t set_type $d(type)
#		$t delete
#	    } else {
#		set t [$io get_anno_ele $rec]
#
#		if {[$t get_comment] != $d(anno)} {
#		    $t set_comment $d(anno)
#		}
#		if {[$t get_type] != $d(type)} {
#		    $t set_type $d(type)
#		}
#		$t delete
#	    }
	    U_tag_change $w $rec [array get d]

	    $w redraw
	    destroy $f
	    unset d
	}

	"quit" {
	    destroy $f
	    unset d
	}

	"move" {
	    puts "move"
	}

	"copy" {
	    puts "copy"
	}

	default {
	    puts "ERROR: Unknown tag callback command: $cmd"
	}
    }
}

proc editor_menu {w x y} {
    upvar $w opt

    puts $opt(top).menubar.commands
    tk_popup $opt(top).menubar.commands \
	[expr $x+[winfo rootx $w]] [expr $y+[winfo rooty $w]]
}

#-----------------------------------------------------------------------------
# Trace display
proc show_trace {w loc} {
    set io [$w io]
    
    if {$loc == ""} return
    foreach {type rec pos} $loc break;

    if {$type == 18} {
	set s [$io get_seq $rec]
	set name [$s get_name]
	set t [trace_add $w.traces $name $w $name]
	after 100 "$t xview 0"

	global $w
	lappend ${w}(Traces) $t
    }
}

#-----------------------------------------------------------------------------
# Generic bindings
bind Editor <Any-Enter> {
    focus %W
}

bind EdNames <Any-Motion> {update_brief %W 1 @%x @%y}

bind Editor <Any-Motion> {update_brief %W 0 @%x @%y}

# Jump to read-pair
bind EdNames <3> {
    global %W

    set ed [set %W(ed)]
    foreach {type rec pos} [$ed get_number @%x @%y] break
    if {![info exists type]} {
	return
    }

    if {$type == 18} {
	set other_end [$ed get_template_seqs $rec]
	if {[llength $other_end] != 1} return
	set s [[$ed io] get_seq $other_end]
	$ed set_cursor 18 $other_end 1
    }
}

bind Editor <<select>> {
    focus %W
    set w [winfo toplevel %W]
    if {![string match [set ${w}(curr_editor)] %W]} {
	set ${w}(curr_editor) %W
	$w.toolbar.undo configure -state [io_undo_state [%W contig_rec]]
#	$w.toolbar.redo configure -state [io_redo_state [%W contig_rec]]
    }
    set _sel [%W get_number @%x @%y 0 1]
    if {$_sel == ""} {
	unset _sel
	return
    } else {
	eval %W set_cursor $_sel 0
    }
    unset _sel
    %W select clear
}

bind Editor <Key-Left>		{%W cursor_left; update_brief %W}
bind Editor <Control-Key-b>	{%W cursor_left; update_brief %W}

bind Editor <Key-Right>		{%W cursor_right; update_brief %W}
bind Editor <Control-Key-f>	{%W cursor_right; update_brief %W}

bind Editor <Key-Up>		{%W cursor_up;    update_brief %W}
bind Editor <Control-Key-p>	{%W cursor_up;    update_brief %W}

bind Editor <Key-Down>		{%W cursor_down;  update_brief %W}
bind Editor <Control-Key-n>	{%W cursor_down;  update_brief %W}

# Not all X11 servers have these keysyms
catch {
    bind Editor <Key-KP_Left>	{%W cursor_left;  update_brief %W}
    bind Editor <Key-KP_Right>	{%W cursor_right; update_brief %W}
    bind Editor <Key-KP_Up>	{%W cursor_up;    update_brief %W}
    bind Editor <Key-KP_Down>	{%W cursor_down;  update_brief %W}
}

bind Editor <Control-Key-a>	{%W read_start;   update_brief %W}
bind Editor <Alt-Key-a>		{%W read_start2;  update_brief %W}
bind Editor <Meta-Key-a>	{%W read_start2;  update_brief %W}
bind Editor <Escape><Key-a>	{%W read_start2;  update_brief %W}

bind Editor <Control-Key-e>	{%W read_end;     update_brief %W}
bind Editor <Alt-Key-e>		{%W read_end2;    update_brief %W}
bind Editor <Meta-Key-e>	{%W read_end2;    update_brief %W}
bind Editor <Escape><Key-e>	{%W read_end2;    update_brief %W}

bind Editor <Alt-Key-comma>	{%W contig_start; update_brief %W}
bind Editor <Meta-Key-comma>	{%W contig_start; update_brief %W}
bind Editor <Escape><Key-comma>	{%W contig_start; update_brief %W}

bind Editor <Alt-Key-period>	{%W contig_end;   update_brief %W}
bind Editor <Meta-Key-period>	{%W contig_end;   update_brief %W}
bind Editor <Escape><Key-period> {%W contig_end;  update_brief %W}

bind Editor <Double-1> {%W display_trace}
bind Editor <Control-Key-t> {%W display_trace}

# Editing commands
bind Editor <Key-a> {editor_edit_base %W a [%W get_number]}
bind Editor <Key-c> {editor_edit_base %W c [%W get_number]}
bind Editor <Key-g> {editor_edit_base %W g [%W get_number]}
bind Editor <Key-t> {editor_edit_base %W t [%W get_number]}
bind Editor <Key-i> {editor_insert_gap %W [%W get_number]}
bind Editor <Key-asterisk> {editor_insert_gap %W [%W get_number]}
bind Editor <Key-Delete> {editor_delete_base %W [%W get_number]}
bind Editor <Key-BackSpace> {editor_delete_base %W [%W get_number]}

bind Editor <Control-Key-Left>  {editor_move_seq %W [%W get_number] -1}
bind Editor <Control-Key-Right> {editor_move_seq %W [%W get_number]  1}

bind Editor <Control-Key-q>     {editor_toggle_annos %W}

bind Editor <Key-less>          {editor_clip_seq %W [%W get_number] l}
bind Editor <Key-greater>       {editor_clip_seq %W [%W get_number] r}

# MouseWheel scrolling
bind Editor  <MouseWheel> {%W yview scroll [expr {-(%D)}] units}
bind EdNames <MouseWheel> {%W yview scroll [expr {-(%D)}] units}
if {[tk windowingsystem] eq "x11"} {
    bind Editor <4>               {%W yview scroll  -1 units}
    bind Editor <5>               {%W yview scroll  +1 units}
    bind Editor <Control-4>       {%W yview scroll -10 units}
    bind Editor <Control-5>       {%W yview scroll +10 units}

    bind Editor <Shift-4>         {%W xview scroll  -1 units}
    bind Editor <Shift-5>         {%W xview scroll  +1 units}
    bind Editor <Shift-Control-4> {%W xview scroll -10 units}
    bind Editor <Shift-Control-5> {%W xview scroll +10 units}

    bind EdNames <4>               {%W yview scroll  -1 units}
    bind EdNames <5>               {%W yview scroll  +1 units}
    bind EdNames <Control-4>       {%W yview scroll -10 units}
    bind EdNames <Control-5>       {%W yview scroll +10 units}

    bind EdNames <Shift-4>         {%W xview scroll  -1 units}
    bind EdNames <Shift-5>         {%W xview scroll  +1 units}
}

bind Editor <Key-Page_Down> {%W xview scroll  +1000 units}
bind Editor <Key-Page_Up>   {%W xview scroll  -1000 units}

# Selection control for adding tags
bind Editor <<select-drag>> {%W select to @%x}
#bind Editor <<select-release>> {puts "select release"}

bind EdNames <<select-drag>> {editor_name_select %W [%W get_number @%x @%y]}

# Tag editing
bind Editor <Key-F11> {tag_editor_launch %W [%W get_number @%x @%y]}
bind Editor <Key-F12> {tag_editor_delete %W [%W get_number @%x @%y]}

bind Editor <<menu>> {
    set _sel [%W get_number @%x @%y 0 1]
    if {$_sel == ""} {
	unset _sel
	return
    } else {
	eval %W set_cursor $_sel
    }
    unset _sel

    tag_repopulate_menu %W
    editor_menu %W %x %y
}