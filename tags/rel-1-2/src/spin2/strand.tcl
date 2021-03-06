proc strand {path args} {

    set b1 top
    set b2 bottom
    radiolist $path \
	-title   Strand \
	-default 1 \
	-orient horizontal \
	-buttons [format {{%s} {%s}} \
		      [list $b1] [list $b2]]

    eval radiolist_configure $path $args
    return $path
}

proc strand_both {path args} {

    set b1 top
    set b2 bottom
    set b3 both
    radiolist $path \
	-title   Strand \
	-default 1\
	-orient horizontal \
	-buttons [format {{%s} {%s} {%s}} \
		      [list $b1] [list $b2] [list $b3]]
    eval radiolist_configure $path $args
    return $path
}

proc strand_get {path} {
    global TOP_S BOTTOM_S

    if {[radiolist_get $path] == 1} {
	return $TOP_S
    } elseif {[radiolist_get $path] == 2} {
	return $BOTTOM_S
    } else {
	return [expr $TOP_S|$BOTTOM_S]
    }
}