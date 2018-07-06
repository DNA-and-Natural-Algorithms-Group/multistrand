// Some simulation datapoints
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
	int numActiveNucl = 0; // for co-transcriptional folding.

	// inspection needs to set this to a non-random variable
	double rchoice = 0.0;

private:

	SimOptions* simOptions = NULL;

};


#endif
