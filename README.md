# sk_s17V6
finding valid 17 clues sudokus in a solution grid
this project is an attempt to have a code faster than the code in the project sk_s17.

A new project is open due to the significant changes done.

The code loaded remains focused on the distribution 665 of the clues in bands and stacks.

The main features changed toward the previous project are the followings

creation of an external loop taking advantage of small unavoidable sets loaded for bands 1+2
Split of the filter on UAs bands 1+2 in a lot of the 64 first and others, immediate switch to vectors at the start

Use of the stack constraint (valid band in stack) through the list of relevant uas

Re use as much as possible of the uas seen in the final control.

A better performance is expected, especially when the number of bands 3 attached to a pair bands1 ,band2 is small.
