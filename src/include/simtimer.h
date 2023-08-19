/*
Multistrand nucleic acid kinetic simulator
Copyright (c) 2008-2023 California Institute of Technology. All rights reserved.
The Multistrand Team (help@multistrand.org)
*/

#ifndef __SIMTIMER_H__
#define __SIMTIMER_H__

class SimOptions;

class SimTimer {

public:
	SimTimer(SimOptions& myOptions);

	void advanceTime(void);
	bool wouldBeHit(const double);
	bool checkHit(const double);
	int checkHitBi(const double collisionRate);

	bool checkForNewNucleotide(void);
	friend std::ostream& operator<<(std::ostream&, SimTimer&);

	double rate = 0.0;
	double stime = 0.0;
	double maxsimtime = 0.0;
	double last_trajectory_time = 0.0;
	long stopcount = 0;
	long stopoptions = 0;

	int nuclAdded = 0;

	// inspection needs to set this to a non-random variable
	double rchoice = 0.0;


	SimOptions* simOptions = NULL;

};

#endif
