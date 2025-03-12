Generally try to use best fits for previous concentration as starting guess for next concentration

If dense phase or fits = False (i.e., too many iterations) then use a file called initial_guess_currentconc rather than best_guess_previousconc

Only for D262V do I run the check for if slope drops drastically only fit through where slope drops

Started trying to add additional days but Day 3 is very noisy for the mutants so currently am focusing on Day 1 and Day 2, have now added Day 4

Tried to start each D262N Day 2 with fibril1 and fibril2 equal but still always made fibril2 go towards zero

Redid WT with same starting as the mutants so all besides day 2 have same initial conditions

Day 4 starts at time 5 but still using same fitting script
