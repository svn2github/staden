#
# Copyright (c) Medical Research Council, Laboratory of Molecular Biology,
# 1998. All rights reserved.
#
# This file is part of the Staden Package. See the Staden Package copyright
# notice for information on the restrictions for usage and distribution, and
# for a disclaimer of all warranties.
#

#
# This handles creating a config file for fetching data from a simple
# 'database' stored in a UNIX file system. The database consists of a fixed
# number of columns with each column representing one piece of information
# The search key is the leftmost column (numbered 0) which will always be
# the sequence identifier.
#
# Eg:
# # ID            ST      PR      SP      TN        SI
# xb54a3.s1       1       1       41      xb54a3    1400..2000
# xb54b12.s1      2       1       41      xb54b12   1400..2000
# xb54b12.r1      2       2       -24     xb54b12   1400..2000
# xb54b12.r1L     2       2       -24     xb54b12   1300..2000

proc unix_database {} {
    global line_types conf_file

    set w .db

    if {[modal $w] == ""} {
	return
    }
    wm title $w "Simple text database"

    # Create windows
    get_fname $w.file \
	-text "Database filename" \
	-type load
    pack $w.file -side top -fill both

    frame $w.adframe
    button $w.adframe.add -text "Add column" \
	-command "unix_database_add $w"
    button $w.adframe.del -text "Delete column" \
	-command "unix_database_del $w"
    pack $w.adframe.add $w.adframe.del -side left
    pack $w.adframe -side top -fill both

    frame $w.lines -bd 2 -relief sunken
    pack $w.lines -side top -fill both -expand 1

    okcancelhelp $w.buttons \
	-ok_command "unix_database_generate $w" \
	-cancel_command "unix_database_cancel $w" \
	-help_command "show_help pregap4 {Pregap4-Database-Simple}"
    pack $w.buttons -side top -fill both

    # Initialise data structures
    upvar #0 $w data
    if {![info exists data]} {
        set lines ""
        foreach i [lsort -ascii $line_types] {
	    if {$i != "ID" || $i != "EN" || $i != "SQ"} {
	        lappend lines $i
	    }
        }
       set data(line_types) $lines
    }

    # Find existing config
    set data(num_lines) 0
    set conf [get_conf_data_section unix_database]
    set cur_lines ""
    foreach l [split $conf \n] {
	if {[regexp "^catch {load_db .*}$" $l]} {
	    $w.file insert end [lindex [lindex $l 1] 1]
	}
	if {[regsub "^proc (..)_com {}.*$" $l {\1} var]} {
	    lappend cur_lines $var
	    # Delete the procedure, incase we delete it from the GUI
	    # and then hit OK.
	    # We reinitialise the memory config when OK or cancel is hit
	    # The catch is because we may have defined the same command twice.
	    catch {rename ${var}_com {}}
	}
    }
    foreach l $cur_lines {
	unix_database_add $w $l
    }

}

proc unix_database_add {w {def {}}} {
    upvar #0 $w data

    set i [incr data(num_lines)]
    label $w.lines.l$i -text "Line type"
    eval tk_optionMenu $w.lines.m$i ${w}(line$i) $data(line_types)
    if {"$def" != ""} {
	set data(line$i) $def
    }
    label $w.lines.f$i -text "Column number $i"
    grid $w.lines.l$i $w.lines.m$i $w.lines.f$i -row $i -padx 5 -pady 5
}

proc unix_database_del {w} {
    upvar #0 $w data
    set i $data(num_lines)
    if {$i > 0} {
	incr data(num_lines) -1
	destroy $w.lines.l$i $w.lines.m$i $w.lines.f$i
    }
}

# Generates the Tcl code for the unix-style database. Note that this code
# contains strings containing tcl code, so it may be confusing. Search for
# "End of" to see the extents of the function.
proc unix_database_generate {w} {
    global pregap4_defs
    upvar #0 $w data

    vfuncheader "Generate database configuration."

    # Consistency check
    for {set lineno 1} {$lineno <= $data(num_lines)} {incr lineno} {
	set l $data(line$lineno)
	if {[info exists check($l)]} {
	    tk_messageBox \
		-icon error \
		-title Error \
		-message "Line type $l appears in more than one column." \
		-type ok \
		-parent $w
	    return
	}
	set check($l) 1
    }
    catch {unset check}

    # Clear if the number of database columns is zero
    if {$data(num_lines) == 0} {
	vmessage "Cleared configuration in config file."
        set_conf_data_section unix_database ""
        write_conf_file
        destroy $w

	if {[winfo exists [keylget pregap4_defs CONFIG.WIN]]} {
	    config_panel_restart
	}

	return
    }

    # More consistency checks...
    if {[$w.file get2] == ""} {
	tk_messageBox \
		-icon error \
		-title Error \
		-message "No filename has been entered." \
		-type ok \
		-parent $w
	return
    }

    # Create the configuration...
    set new_conf \
{# ------------------------------------------------------------------
# This has been automatically generated by Pregap4. The following
# procedure is used to load up the database file into memory when
# pregap4 is started. The subsequent *_com procedures query this
# database during the augment stage.

proc load_db {db_file} {
    global db_arr
    set fd [open $db_file r]
    while {[gets $fd line] != -1} {
	set db_arr([lindex $line 0]) [lrange $line 1 end]
    }
    close $fd
}
proc find_item {type index} {
    global lines db_arr
    set id $lines(ID)
    if {[info exists db_arr($id)]} {
	set l $db_arr($id)
    } else {
	set l ""
	foreach i [array names db_arr] {
	    if {[string match $i $id]} {
		set l $db_arr($i)
		break
	    }
	}
	set db_arr($id) $l
    }
    return [set lines($type) [lindex $l $index]]
}
}

    set db_file [$w.file get2]
    append new_conf "catch {load_db [list $db_file]}\n\n"

    for {set i 1} {$i <= $data(num_lines)} {incr i} {
	set n $data(line$i)
	append new_conf \
	    "proc ${n}_com {} {return \[find_item $n [expr $i-1]\]}\n"
	vmessage "Wrote ${n}_com to config file."
    }

    set_conf_data_section unix_database $new_conf
    write_conf_file

    destroy $w
    uplevel #0 {eval [get_conf_data_section unix_database]}

    if {[winfo exists [keylget pregap4_defs CONFIG.WIN]]} {
	config_panel_restart
    }
}
# End of unix_database_generate function

proc unix_database_cancel {w} {
    destroy $w
    uplevel #0 {eval [get_conf_data_section unix_database]}
}